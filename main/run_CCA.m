clear all;

malexflag = 1; % user flag
if malexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = path.code; % save directory
    
    %Meryem Laptop
    %     path.code = 'C:\Users\m\Documents\GitHub\GLM-BCI'; addpath(genpath(path.code)); % code directory
    %     path.dir = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    %     path.save = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = path.code; % save directory
end

% #####
%% simulated data file names
filename = 'resting_sim_100_shorterHRF';
%% load ground truth hrf
hrf = load([path.code '\sim HRF\hrf_simdat_100_shorterHRF.mat']);
%% save folder name
sfoldername = '\CV_results_data';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};

%% Options/Parameter Settings
rhoSD_ssThresh = 15;  % mm
flag_save = 0;
flag_conc = 1; % if 1 CCA inputs are in conc, if 0 CCA inputs are in intensity
% results eval parameters
eval_param.HRFmin = -2;
eval_param.HRFmax = 15; % used only for block design runs
eval_param.Hb = 1; % 1 HbO / 0 HbR (for block only)
eval_param.pre = 5;  % HRF range in sec to calculate ttest
eval_param.post = 8;
flag_detrend = 0; % input paramater to load_nirs function: performing linear detrend if 1, no detrending if 0 during "pre-processing"
drift_term = 1; % input parameter to hmrDeconvHRF_DriftSS function: performing linear detrend for GLM_SS and GLM_CCA during single trial estimation
polyOrder_drift_hrfestimate = 3; % input parameter to hmrDeconvHRF_DriftSS function: polynomial order, performs linear/polynomial detrending during estimation of HRF from training data
flag_hrf_resid = 0; % 0: hrf only; 1: hrf+yresid
% CCA parameters
flags.pcaf =  [0 0]; % no pca of X or AUX
flags.shrink = true;
% perform regularized (rtcca) (alternatively old approach)
rtccaflag = true;

% Features/structs for feature extraction function
fparam.swdw=[0,4;8,15]; % need to discuss this selection!
ival = [eval_param.HRFmin eval_param.HRFmax];

% get features from ground truth
hrfdat.x = hrf.hrf_conc;
hrfdat.fs = 25;
hrfdat.t = hrf.t_hrf';
[FVgt] = featureExtract(hrfdat, fparam);

% motion artifact detection
motionflag = false;
% plot flag
flag_plot = true;
% include tcca results or not in plots?
flag_plotCCA = true;


% Validation parameters
% tlags = 0:1:10;
% stpsize = 2:2:24;
% cthresh = 0:0.1:0.9;
tlags = 3;
stpsize = 2;
cthresh = 0.3;
%% CCA with optimum parameters
tl = tlags;
sts = stpsize;
ctr = cthresh;

% set stepsize for CCA
param.tau = sts; %stepwidth for embedding in samples (tune to sample frequency!)
% set correlation trheshold for CCA to 0 so we dont lose anything within
% the tCCA function
param.ct = 0;   % correlation threshold
% set number of tCCA regressors used for GLM
tcca_nReg = 2;


%% Eval plot flag (developing/debugging purposes only)
evalplotflag = 0; % compares dc, hrf_ss, hrf_tcca, true hrf for hrf added channels
evalplotflag_glm = 0; % displays raw signal, model fit, yresid, hrf, ss, drift etc for sanity check (now for glm_ss only)

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
    
    %% lowpass filter AUX signals
    AUX = hmrBandpassFilt(AUX, fq, 0, 0.5);
    %% AUX signals
    AUX = [AUX, d0_short]; % full AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
    %% zscore AUX signals
    AUX = zscore(AUX);
    
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
    spltIDX = {1:floor(len/5),floor(len/5)+1:len};
    trntst = {[2,1], [1,2]};
    
    %% run test and train CV splits
    for tt = 1%1:2
        tstIDX = spltIDX{trntst{tt}(1)};
        trnIDX = spltIDX{trntst{tt}(2)};
        
        %% convert testing fNIRS data to concentration and detect motion artifacts
        dod = hmrIntensity2OD(d(tstIDX,:));
        s = s(tstIDX,:);
        t = t(tstIDX,:);
        % motion artifact removal?
        if motionflag
            [tIncAuto] = hmrMotionArtifact(dod,fq,SD,ones(size(d,1),1),0.5,1,30,5);
            [s,tRangeStimReject] = enStimRejection(t,s,tIncAuto,ones(size(d,1),1),[-2  10]);
        end
        % conversion
        dod = hmrBandpassFilt(dod, fq, 0, 0.5);
        dc{tt} = hmrOD2Conc( dod, SD, [6 6]);        
        dc_linear_detrend{tt} = linear_detrend(dc{tt}, t); % keeping this for no GLM
        %% run test and train CV splits
        onset_stim = find(s(:,1)==1); % condition 1: stimulus
        onset_stim_rest = find(s(:,2)==1); % condition 2: rest
        
        % in case of a really early stim not allowing prestim period
        if onset_stim(1) < abs(eval_param.HRFmin*fq)
            onset_stim(1) = [];
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
        
        for os = 1:size(onset_stim,1)%% loop around each stimulus
            % write train/test flag matrix
            TTM{sbj}.tstidx(os)=os;
            TTM{sbj}.tnridx(os,:) = setdiff(1:size(onset_stim,1), os);
            
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
            %% estimate HRF regressor with GLM with tCCA from all training trials
            % *****************************************************
            %% learn tCCA filters
            param.NumOfEmb = ceil(tl*fq / sts);
            % Temporal embedding of auxiliary data from training split to
            % learn the filter matrix
            % Perform CCA on training data % AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
            % use test data of LD channels without synth HRF
            X = d0_long(trnIDX,:);
            % new tCCA with shrinkage
            [REG_trn,  ADD_trn] = rtcca(X,AUX(trnIDX,:),param,flags);
            % save trained tCCA filter matrix (first tcca_nReg regressors)
            Atcca = ADD_trn.Av(:,1:tcca_nReg);
            
            %% generate tCCA regressor and zero elements from test trial pre-stimulus to the end of the following rest block
            aux_emb = tembz(AUX(tstIDX,:), param);
            % calculate tcca regressors
            REG_tcca = aux_emb*Atcca;     
                    
            % zero out test blocks
            REG_tcca(pre_stim_t{1}(os):post_stim_t{2}(os),:) = 0;
            %% Generate HRF regressor with GLM + tCCA and the tCCA regressor 
            [HRF_regressor_CCA, yavgstd_ss, tHRF, nTrialsSS, d_ss, yresid_ss, ysum2_ss, beta, yR_ss] = ...
                hmrDeconvHRF_DriftSS(dc_new, s_new, t, SD, REG_tcca, [], [eval_param.HRFmin eval_param.HRFmax], ...
                1, 1, [0.5 0.5], rhoSD_ssThresh, 1, polyOrder_drift_hrfestimate, 0,hrf,lstHrfAdd{sbj},0, ...
                [pre_stim_t{1}(os)  post_stim_t{2}(os)]); %% MERYEM
            
            
            % *****************************************************
            %% calculate single trial HRF estimates using HRF regressor and pre-used single training trials, as well as unseen test block
            % *****************************************************
            % results are used for feature-extraction and classifier
            % training + testing, while statistics of unseen test block are
            % only used for generating the test feature vector for the
            % classifier
            for k = 1:size(onset_stim,1)%% loop around each stimulus
                
                %% Save normal raw data (single trials)
                for cc=1:2
                    y_raw(:,:,:,k,cc)= dc_linear_detrend{tt}(pre_stim_t{cc}(os):post_stim_t{cc}(os),:,:);
                end
                
                %% Perform GLM on single trials
                % for both conditions cc = 1 (HRF) and cc = 2 (rest)
                for cc=1:2
                    %% GLM with SS
                    [yavg_ss(:,:,:,k,cc), yavgstd_ss, tHRF, nTrialsSS, ynew_ss(:,:,:,k,cc), yresid_ss, ysum2_ss, beta_ss(:,:,:,k,cc), yR_ss] = ...
                        hmrDeconvHRF_DriftSS(dc{tt}(pre_stim_t{cc}(k):post_stim_t{cc}(k),:,:), ...
                        s(pre_stim_t{cc}(k):post_stim_t{cc}(k),cc), ...
                        t(pre_stim_t{cc}(k):post_stim_t{cc}(k),:), ...
                        SD, [], [], [eval_param.HRFmin eval_param.HRFmax], 1, 5, ...
                        squeeze(HRF_regressor_SS(:,:,:,cc)), rhoSD_ssThresh, 1, drift_term, ...
                        0, hrf,lstHrfAdd{sbj},evalplotflag_glm,[] );
                    
                    %% Perform GLM with CCA
                    [yavg_cca(:,:,:,k,cc), yavgstd_cca, tHRF, nTrials(tt), ynew_cca(:,:,:,k,cc), yresid_cca, ysum2_cca, beta_cca(:,:,:,k,cc), yR_cca] = ...
                        hmrDeconvHRF_DriftSS(dc{tt}(pre_stim_t{cc}(k):post_stim_t{cc}(k),:,:), ...
                        s(pre_stim_t{cc}(k):post_stim_t{cc}(k),cc), ...
                        t(pre_stim_t{cc}(k):post_stim_t{cc}(k),:), ...
                        SD, REG_tcca(pre_stim_t{cc}(k):post_stim_t{cc}(k),:), [], [eval_param.HRFmin eval_param.HRFmax], 1, 5, ...
                        squeeze(HRF_regressor_CCA(:,:,:,cc)), 0, 0, drift_term, 0, hrf,lstHrfAdd{sbj},0,[]);
                end
                
                % use HRF + residual?
                if flag_hrf_resid
                    yavg_ss = ynew_ss;
                    yavg_cca = ynew_cca;
                end
                
                if evalplotflag  % plotting all hrf added channels for a single subject
                    a=dc{tt}([onset_stim(os)+eval_param.HRFmin*fq]:[onset_stim(os)+eval_param.HRFmax*fq],:,:)-repmat(mean(dc{tt}([onset_stim(os)+eval_param.HRFmin*fq]:onset_stim(os),:,:),1),numel([onset_stim(os)+eval_param.HRFmin*fq]:[onset_stim(os)+eval_param.HRFmax*fq]),1);
                    figure;subplot(1,3,1);plot(tHRF,squeeze(a(:,1,lstHrfAdd{sbj}(:,1))));ylim([-1e-6 1.5e-6]);title('none'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,1),'k','LineWidth',2);ylabel('HbO (hrf only)');
                    subplot(1,3,2);plot(tHRF,squeeze(yavg_ss(:,1,lstHrfAdd{sbj}(:,1),os))); ylim([-1e-6 1.5e-6]);title('ss'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,1),'k','LineWidth',2);
                    subplot(1,3,3);plot(tHRF,squeeze(yavg_cca(:,1,lstHrfAdd{sbj}(:,1),os)));ylim([-1e-6 1.5e-6]);title('cca'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,1),'k','LineWidth',2);
                    figure;subplot(1,3,1);plot(tHRF,squeeze(a(:,2,lstHrfAdd{sbj}(:,1))));ylim([-1e-6 0.5e-6]);title('none'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,2),'k','LineWidth',2);ylabel('HbR (hrf only)');
                    subplot(1,3,2);plot(tHRF,squeeze(yavg_ss(:,2,lstHrfAdd{sbj}(:,1),os)));ylim([-1e-6 0.5e-6]); title('ss'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,2),'k','LineWidth',2);
                    subplot(1,3,3);plot(tHRF,squeeze(yavg_cca(:,2,lstHrfAdd{sbj}(:,1),os)));ylim([-1e-6 0.5e-6]);title('cca'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,2),'k','LineWidth',2);
                    
                    figure;subplot(1,3,1);plot(tHRF,squeeze(a(:,1,lstHrfAdd{sbj}(:,1))));ylim([-1e-6 1.5e-6]);title('none'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,1),'k','LineWidth',2);ylabel('HbO (hrf+resid)');
                    subplot(1,3,2);plot(tHRF,squeeze(ynew_ss(:,1,lstHrfAdd{sbj}(:,1),os))); ylim([-1e-6 1.5e-6]);title('ss'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,1),'k','LineWidth',2);
                    subplot(1,3,3);plot(tHRF,squeeze(ynew_cca(:,1,lstHrfAdd{sbj}(:,1),os)));ylim([-1e-6 1.5e-6]);title('cca'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,1),'k','LineWidth',2);
                    figure;subplot(1,3,1);plot(tHRF,squeeze(a(:,2,lstHrfAdd{sbj}(:,1))));ylim([-1e-6 0.5e-6]);title('none'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,2),'k','LineWidth',2);ylabel('HbR (hrf+resid)');
                    subplot(1,3,2);plot(tHRF,squeeze(ynew_ss(:,2,lstHrfAdd{sbj}(:,1),os)));ylim([-1e-6 0.5e-6]); title('ss'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,2),'k','LineWidth',2);
                    subplot(1,3,3);plot(tHRF,squeeze(ynew_cca(:,2,lstHrfAdd{sbj}(:,1),os)));ylim([-1e-6 0.5e-6]);title('cca'); hold on; plot(hrf.t_hrf,hrf.hrf_conc(:,2),'k','LineWidth',2);
                end
                % display current state:
                disp(['sbj ' num2str(sbj) ', epoch ' num2str(os) ])
            end
            
            %% get features/markers
            % stim on blocks
            % for both conditions
            for cc=1:2
                % no GLM
                [FMdc{sbj,os}(:,:,:,:,cc), clab] = getFeaturesAndMetrics(y_raw(:,:,:,:,cc), fparam, ival, hrf);
                % short separation GLM
                [FMss{sbj,os}(:,:,:,:,cc), clab] = getFeaturesAndMetrics(yavg_ss(:,:,:,:,cc), fparam, ival, hrf);
                % tCCA GLM
                [FMcca{sbj,os}(:,:,:,:,cc), clab] = getFeaturesAndMetrics(yavg_cca(:,:,:,:,cc), fparam, ival, hrf);
            end
            %% save weights
            FWss{sbj,os}=beta_ss;
            FWcca{sbj,os}=beta_cca;
        end
    end
    % clear vars
    clear vars AUX d d0 d_long d0_long d_short d0_short t s REG_trn ADD_trn
end

%% save data
if flag_save
    disp('saving data...')
    save([path.save '\FV_results.mat'], 'FMdc', 'FMss', 'FMcca', 'FWss', 'FWcca', 'TTM', 'lstHrfAdd', 'lstLongAct', 'lstShortAct');
end

toc;







