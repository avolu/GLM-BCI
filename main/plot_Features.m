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
%% load data
load([path.save '\FV_results.mat'])
%% load ground truth hrf
hrf = load([path.code '\sim HRF\hrf_simdat_100_shorterHRF.mat']);

% include tcca results in plots?
flag_plotCCA = true;

% Features/structs for feature extraction function
eval_param.HRFmin = -2;
eval_param.HRFmax = 15; % used only for block design runs
fparam.swdw=[0,4;8,11]; % need to discuss this selection!
ival = [eval_param.HRFmin eval_param.HRFmax];

% get features from ground truth
hrfdat.x = hrf.hrf_conc;
hrfdat.fs = 25;
hrfdat.t = hrf.t_hrf';
[FVgt] = featureExtract(hrfdat, fparam);

%% Get HRF features from all augmented channels to compare against  ground truth
%Sort through results and append
FV_Raw = cell(2);
FV_SS = cell(2,2);
FV_CCA = cell(2,2);
for sbj = 1:numel(TTM)
    % only look at the crossvalidated test results. These are stored where cell and trial index
    % coincide (== os)
    for os = 1:numel(TTM{sbj}.tstidx)
        for cc=1:2 % stim and resting condition
            % channel indices that have or dont have gt HRF
            idxChHrf = lstHrfAdd{sbj}(:,1);
            idxChNoHrf = setdiff(lstLongAct{sbj},squeeze(lstHrfAdd{sbj}(:,1)));
            if size(idxChHrf,1) > size(idxChNoHrf,1)
                idxChHrf = idxChHrf(1:size(idxChNoHrf,1));
            else
                idxChNoHrf = idxChNoHrf(1:size(idxChHrf,1));
            end
            % number of available channels
            nHrf = size(FMdc{sbj}(:,:,idxChHrf,:,cc));
            nNoHrf = size(FMdc{sbj}(:,:,idxChNoHrf,:,cc));
            % extract and append crossvalidated features (from testing trial), new dimension is F x C x I,
            % where F: # of Features, C: # Number of Chromophores, I: # of all
            % trials (epochs*channels)            
            FV_Raw{cc} = cat(3, FV_Raw{cc}, FMdc{sbj}(:,:,idxChHrf,os,cc));
            for rr=1:2 % stim and resting hrf regressor
                FV_SS{cc,rr} = cat(3, FV_SS{cc,rr}, FMss{sbj,os}(:,:,idxChHrf,os,cc,rr));
                FV_CCA{cc,rr} = cat(3, FV_CCA{cc,rr}, FMcca{sbj,os}(:,:,idxChHrf,os,cc,rr));
            end
        end
    end
end


%% Paired T-Tests for STIM trials with HRF added and HRF STIM regressor
cc = 1;
rr = 1;
for ff = 1:9
    for ch=1:3
        [h_co(ff,ch,1),p_co(ff,ch,1)]= ttest(squeeze(FV_Raw{cc,rr}(ff,ch,:)),squeeze(FV_SS{cc,rr}(ff,ch,:)));
        [h_co(ff,ch,2),p_co(ff,ch,2)]= ttest(squeeze(FV_Raw{cc,rr}(ff,ch,:)),squeeze(FV_CCA{cc,rr}(ff,ch,:)));
        [h_co(ff,ch,3),p_co(ff,ch,3)]= ttest(squeeze(FV_SS{cc,rr}(ff,ch,:)),squeeze(FV_CCA{cc,rr}(ff,ch,:)));
    end
end

%% (Box)Plot results (normal metrics) for STIM trials with HRF added and HRF STIM regressor
cc = 1;
rr = 1;
figure
labels = {'No GLM', 'GLM SS', 'GLM tCCA'};
chrom = {' HbO', ' HbR'};
% for all features
for ff=1:9
    % for both chromophores
    for ch=1:2
        subplot(2,9,(ch-1)*9+ff)
        xtickangle(35)
        hold on
        %% boxplots
        % with cca
        if flag_plotCCA
            boxplot([squeeze(FV_Raw{cc,rr}(ff,ch,:)), squeeze(FV_SS{cc,rr}(ff,ch,:)), squeeze(FV_CCA{cc,rr}(ff,ch,:))], 'labels', labels)
            H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,ch,1:3)));
        else
            % without cca
            boxplot([squeeze(FV_Raw{cc,rr}(ff,ch,:)), squeeze(FV_SS{cc,rr}(ff,ch,:))], 'labels', labels(1:2))
            H=sigstar({[1,2]},squeeze(p_co(ff,ch,1)));
        end
        
        % ground truth
        if ff<8
            hAx=gca;                                   % retrieve the axes handle
            %xtk=hAx.XTick;                             % and the xtick values to plot() at...
            xtk = [0.5, 1, 2, 3, 3.5];
            hold on
            if ff==5
                FVgt.x(ff,ch) = 6; % set time to peak to 6 seconds (due to gt hrf plateau)
            end
            hL=plot(xtk, ones(numel(xtk))*FVgt.x(ff,ch),'--g');
        end
        if ff<5
            ylabel('\muMol')
        end
        if ff==5
            ylabel('sec')
        end
        if ff==9
            ylabel('Mol')
        end
        
        title([FMclab{ff} chrom{ch}])
    end
end


%% (Box)Plot results (metric errors)
figure
labels = {'No GLM', 'GLM SS', 'GLM tCCA'};
chrom = {' HbO', ' HbR'};
% for all features
for ff=1:9
    % for both chromophores
    for ch=1:2
        subplot(2,9,(ch-1)*9+ff)
        xtickangle(35)
        hold on
        %% boxplots
        if ff==5
            FVgt.x(ff,ch) = 6; % set time to peak to 6 seconds (due to gt hrf plateau)
        end
        
        if ff<8
            % without GLM, with GLM+SS, with GLM+CCA
            if flag_plotCCA
                boxplot([abs(squeeze(FV_Raw{cc,rr}(ff,ch,:))-FVgt.x(ff,ch)), abs(squeeze(FV_SS{cc,rr}(ff,ch,:))-FVgt.x(ff,ch)), abs(squeeze(FV_CCA{cc,rr}(ff,ch,:))-FVgt.x(ff,ch))], 'labels', labels)
                H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,ch,1:3)));
            else
                boxplot([abs(squeeze(FV_Raw{cc,rr}(ff,ch,:))-FVgt.x(ff,ch)), abs(squeeze(FV_SS{cc,rr}(ff,ch,:))-FVgt.x(ff,ch))], 'labels', labels(1:2))
                H=sigstar({[1,2]},squeeze(p_co(ff,ch,1)));
            end
            title(['ERR ' FMclab{ff} chrom{ch}])
        else
            if flag_plotCCA
                boxplot([squeeze(FV_Raw{cc,rr}(ff,ch,:)), squeeze(FV_SS{cc,rr}(ff,ch,:)), squeeze(FV_CCA{cc,rr}(ff,ch,:))], 'labels', labels)
                H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,ch,1:3)));
            else
                boxplot([squeeze(FV_Raw{cc,rr}(ff,ch,:)), squeeze(FV_SS{cc,rr}(ff,ch,:))], 'labels', labels(1:2))
                H=sigstar({[1,2]},squeeze(p_co(ff,ch,1)));
            end
            
            title([FMclab{ff} chrom{ch}])
        end
        
        
        if ff<5
            ylabel('\muMol')
        end
        if ff==5
            ylabel('sec')
        end
        if ff==9
            ylabel('\muMol')
        end
    end
end


