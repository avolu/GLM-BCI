clear all;


malexflag = 0; % user flag
if malexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\GLM_BCI_PAPER\RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\GLM_BCI_PAPER\PROCESSED_DATA\'; % save directory
    
    %Meryem Laptop
    %     path.code = 'C:\Users\m\Documents\GitHub\GLM-BCI'; addpath(genpath(path.code)); % code directory
    %     path.dir = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    %     path.save = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\GLM_BCI_PAPER\RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\GLM_BCI_PAPER\PROCESSED_DATA\'; % save directory
end

% #####
% 50 or 100% hrf?
hrf_amp = 100;

%% simulated data file names %% load ground truth hrf
switch hrf_amp
    case 50
        filename = 'resting_sim_50_shorterHRF';
        hrf = load([path.code '\sim HRF\hrf_simdat_50_shorterHRF.mat']);
    case 100
        filename = 'resting_sim_100_shorterHRF';
        hrf = load([path.code '\sim HRF\hrf_simdat_100_shorterHRF.mat']);
end

% single trial linear detrending for no-GLM approach?
stlindet = true;

%% save folder name
sfoldername = '\CV_results_data';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj86', 'Subj91', 'Subj92', 'Subj94', 'Subj95', 'Subj96', 'Subj97', 'Subj98', 'Subj99', 'Subj100', 'Subj101', 'Subj102','Subj103', 'Subj104'}; % potentially exclude Subj94 (only 4 channels HRF added because of motion)


%% Options/Parameter Settings
rhoSD_ssThresh = 15;  % mm
flag_save = 1;
flag_conc = 1; % if 1 CCA inputs are in conc, if 0 CCA inputs are in intensity
% results eval parameters
eval_param.HRFmin = -2;
eval_param.HRFmax = 15; % used only for block design runs
eval_param.Hb = 1; % 1 HbO / 0 HbR (for block only)
eval_param.pre = 5;  % HRF range in sec to calculate ttest
eval_param.post = 8;
flag_detrend = 0; % input paramater to load_nirs function: performing linear detrend if 1, no detrending if 0 during "pre-processing"
drift_term = 1; % input parameter to hmrDeconvHRF_DriftSS function: performing linear detrend for GLM_SS and GLM_CCA during single trial estimation
polyOrder_drift_hrfestimate = 1; % input parameter to hmrDeconvHRF_DriftSS function: polynomial order, performs linear/polynomial detrending during estimation of HRF from training data
flag_hrf_resid = 0; % 0: hrf only; 1: hrf+yresid

% Features/structs for feature extraction function
fparam.swdw=[0,4;8,11]; % need to discuss this selection!
ival = [eval_param.HRFmin eval_param.HRFmax];

% motion artifact detection
motionflag = false;
motioncorr_flag = false;

%% Eval plot flag (developing/debugging purposes only)
evalplotflag = 0; % compares dc, hrf_ss, hrf_tcca, true hrf for hrf added channels
evalplotflag_glm = 0; % displays raw signal, model fit, yresid, hrf, ss, drift etc for sanity check (now for glm_ss only)

%omit first omsec seconds in training data (to get rid of initial
%motion artifacts)
omsec = 20;

%% dimensions
% beta [# of basis] X [HbO/R] X [# of channels]

%% GO
tic;
for sbj = 1:numel(sbjfolder) % loop across subjects
    disp(['subject #' num2str(sbj)]);
    
    % change to subject directory
    cd([path.dir filesep sbjfolder{sbj} filesep]);
    
    %% load data
    [fq, t, AUX, d_long, d_short, d0_long, d0_short, d, d0, SD, s, lstLongAct{sbj},lstShortAct{sbj},lstHrfAdd{sbj}] = load_nirs(filename,flag_conc,flag_detrend);
    
    %% check if the number of time points is odd/even, if odd make it even... (number of embedded should be the same)
    if mod(size(AUX,1),2) == 1
        AUX(end,:)=[];
        d(end,:)=[];
        d_long(end,:)=[];
        d_short(end,:)=[];
        d0(end,:)=[];
        d0_long(end,:)=[];
        d0_short(end,:)=[];
        t(end,:)=[];
        s(end,:)=[];
    end
    
    % create data split indices
    len = size(AUX,1);
    cvIDX = omsec*fq:len;
    
    %% convert testing fNIRS data to concentration and detect motion artifacts
    dod = hmrIntensity2OD(d(cvIDX,:));
    
    
    if motioncorr_flag
    dod = hmrMotionCorrectSplineSG(dod, d(cvIDX,:), t(cvIDX,:), SD, 0.99, 10, 1);
    end
    
    s = s(cvIDX,:);
    t = t(cvIDX,:);
    % motion artifact removal?
    if motionflag
        [tIncAuto] = hmrMotionArtifact(dod,fq,SD,ones(size(d,1),1),0.5,1,30,5);
        [s,tRangeStimReject] = enStimRejection(t,s,tIncAuto,ones(size(d,1),1),[-2  10]);
    end
    % conversion
    dod = hmrBandpassFilt(dod, fq, 0, 0.5);
    dc = hmrOD2Conc( dod, SD, [6 6]);
    dc_linear_detrend = linear_detrend(dc, t); % keeping this for no GLM
    %% run test and train CV splits
    onset_stim = find(s(:,1)==1); % condition 1: stimulus
    onset_stim_rest = find(s(:,2)==1); % condition 2: rest
    
    % in case of a really early stim not allowing prestim period
    if onset_stim(1) < abs(eval_param.HRFmin*fq)% if first stim is too early
        onset_stim(1) = [];
        onset_stim_rest(1) = [];
    elseif onset_stim_rest(1) < abs(eval_param.HRFmin*fq) % if first rest is too early
        onset_stim_rest(1) = [];
    end
    
    % if after cutting the data into cvIDX, rest mark is earlier than the
    % first stim mark, remove rest (this is hypothetical now)
    if onset_stim_rest(1) < onset_stim(1)
        onset_stim_rest(1) = [];
    end
    
    % in case there is no rest block after the last stim block
    if onset_stim(end)>onset_stim_rest(end)
        onset_stim(end) = [];
    end
    
    % all pre/post stimulus timings for both conditions
    %       cc = 1 (HRF),% cc = 2 (rest)
    pre_stim_t{1} = onset_stim+eval_param.HRFmin*fq;
    post_stim_t{1} = onset_stim+eval_param.HRFmax*fq;
    pre_stim_t{2} = onset_stim_rest+eval_param.HRFmin*fq;
    post_stim_t{2} = onset_stim_rest+eval_param.HRFmax*fq;
    
    % onset_stimulus for both conditions
    onset_stim_all{1} = onset_stim;
    onset_stim_all{2} = onset_stim_rest;
    
    for os = 1:size(onset_stim,1)%% loop around each stimulus
        % write train/test flag matrix
        TTM{sbj}.tstidx(os)=os;
        TTM{sbj}.tnridx(os,:) = setdiff(1:size(onset_stim,1), os);
        
        for cc=1:2
%             %% Save normal raw data in single trials and remove baseline
%             y_raw(:,:,:,os,cc)= dc_linear_detrend(pre_stim_t{cc}(os):post_stim_t{cc}(os),:,:)- ...
%                 repmat(mean(dc_linear_detrend(pre_stim_t{cc}(os):onset_stim_all{cc}(os),:,:),1),numel(pre_stim_t{cc}(os):post_stim_t{cc}(os)),1);
%             %% Linear detrend single trial
%             if stlindet
%                 tdet= t(pre_stim_t{cc}(os):post_stim_t{cc}(os));
%                 y_raw(:,:,:,os,cc) = linear_detrend(y_raw(:,:,:,os,cc), tdet);
%             end
            
            %% Save normal raw data in single trials 
            y_raw(:,:,:,os,cc)= dc_linear_detrend(pre_stim_t{cc}(os):post_stim_t{cc}(os),:,:);
            % Linear detrend single trial
            if stlindet
                tdet= t(pre_stim_t{cc}(os):post_stim_t{cc}(os));
                y_raw(:,:,:,os,cc) = linear_detrend(y_raw(:,:,:,os,cc), tdet);
            end
            % remove baseline 
            bl = repmat(mean(y_raw(1:abs(eval_param.HRFmin*fq),:,:,os,cc),1),numel(1:size(y_raw,1)),1);  
            y_raw(:,:,:,os,cc)= y_raw(:,:,:,os,cc)- bl;
            
            %% save info for t-test based channel selection (do only for HRF condition)
            if cc==1
                chselInfo{sbj}.BL_RAW(os,:,:) = squeeze(mean(y_raw(1:abs(eval_param.HRFmin*fq),:,:,os,cc),1));
                chselInfo{sbj}.PEAK_RAW(os,:,:) = squeeze(mean(y_raw(abs(eval_param.HRFmin*fq)+(4*fq:8*fq),:,:,os,cc),1));
            end
        end
        
        % *****************************************************
        %% estimate HRF regressor with GLM with SS from all training trials
        % *****************************************************
        % zero elements from pre-stimulus to the end of the following rest block
        dod_new = dod; dod_new(pre_stim_t{1}(os):post_stim_t{2}(os),:) = 0;
        dc_new = hmrOD2Conc(dod_new, SD, [6 6]);
        s_new = s; s_new(pre_stim_t{1}(os):post_stim_t{2}(os),:) = 0;
        % GLM with SS: generate HRF regressor
        [HRF_regressor_SS, yavgstd_ss, tHRF, nTrialsSS, d_ss, yresid_ss, ysum2_ss, beta, yR_ss] = ...
            hmrDeconvHRF_DriftSS(dc_new, s_new, t, SD, [], [], [eval_param.HRFmin eval_param.HRFmax], ...
            1, 1, [0.5 0.5], rhoSD_ssThresh, 1, polyOrder_drift_hrfestimate, 0,hrf,lstHrfAdd{sbj},0, ...
            [pre_stim_t{1}(os)  post_stim_t{2}(os)]);
        
        % *****************************************************
        %% calculate single trial HRF estimates using HRF regressor and pre-used single training trials, as well as unseen test block
        % *****************************************************
        % results are used for feature-extraction and classifier
        % training + testing, while statistics of unseen test block are
        % only used for generating the test feature vector for the
        % classifier
        for k = 1:size(onset_stim,1)%% loop around each stimulus
            
            %% Perform GLM on single trials
            % for both conditions cc = 1 (stim) and cc = 2 (rest)
            for cc=1:2
                % for both hrf regressors rr=1 (hrf stim) and rr=2 (hrf rest)
                for rr=1:2
                    %% GLM with SS
                    [yavg_ss(:,:,:,k,cc,rr), yavgstd_ss, tHRF, nTrialsSS, ynew_ss(:,:,:,k,cc,rr), yresid_ss, ysum2_ss, beta_ss(:,:,:,k,cc,rr), yR_ss] = ...
                        hmrDeconvHRF_DriftSS(dc(pre_stim_t{cc}(k):post_stim_t{cc}(k),:,:), ...
                        s(pre_stim_t{cc}(k):post_stim_t{cc}(k),cc), ...
                        t(pre_stim_t{cc}(k):post_stim_t{cc}(k),:), ...
                        SD, [], [], [eval_param.HRFmin eval_param.HRFmax], 1, 5, ...
                        squeeze(HRF_regressor_SS(:,:,:,rr)), rhoSD_ssThresh, 1, drift_term, ...
                        0, hrf,lstHrfAdd{sbj},evalplotflag_glm,[] );
                end
            end
            
            % use HRF + residual?
            if flag_hrf_resid
                yavg_ss = ynew_ss;
%               yavg_cca = ynew_cca;
            end
            
            % display current state:
            disp(['sbj ' num2str(sbj) ', epoch ' num2str(os) ])
        end
        
        %% get features/markers
        % stim on blocks
        % for both conditions
        for cc=1:2
            % for both regressors
            for rr=1:2
                % short separation GLM
                [FMss{sbj,os}(:,:,:,:,cc,rr), FMclab] = getFeaturesAndMetrics(yavg_ss(:,:,:,:,cc,rr), fparam, ival, hrf);
            end
        end
        
        %% save info for t-test based channel selection (do only for HRF condition)
        cc=1;
        rr= 1;
        chselInfo{sbj}.BL_SS(os,:,:,:) = squeeze(mean(yavg_ss(1:abs(eval_param.HRFmin*fq),:,:,:,cc,rr),1));
        chselInfo{sbj}.PEAK_SS(os,:,:,:) = squeeze(mean(yavg_ss(abs(eval_param.HRFmin*fq)+(4*fq:8*fq),:,:,:,cc,rr),1));

        %% save weights
        FWss{sbj,os}=beta_ss;
    end
    %% no GLM: feature extraction
    for cc=1:2
        [FMdc{sbj}(:,:,:,:,cc), FMclab] = getFeaturesAndMetrics(y_raw(:,:,:,:,cc), fparam, ival, hrf);
    end
end
% clear vars
clear vars AUX d d0 d_long d0_long d_short d0_short t s REG_trn ADD_trn


%% save data
if flag_save
    disp('saving data...')
    save([path.save '\FV_results_SSvsNo_ldrift' num2str(drift_term) '_resid' num2str(flag_hrf_resid) 'stlindriftSWAPPED_' num2str(stlindet) '_hrf_amp' num2str(hrf_amp) '_mot_corr' num2str(motioncorr_flag) '_20soffs_ttstchsel.mat'], 'FMdc', 'FMss', 'FWss', 'TTM', 'lstHrfAdd', 'lstLongAct', 'lstShortAct', 'FMclab', 'chselInfo');
end

toc;







