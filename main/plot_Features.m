clear all;

malexflag = 0; % user flag
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
F_Raw_Hrf=[];
F_Raw_NoHrf=[];
F_SS_Hrf=[];
F_SS_NoHrf=[];
F_CCA_Hrf=[];
F_CCA_NoHrf=[];
for sbj = 1:numel(TTM)
    % only look at the crossvalidated test results. These are stored where cell and trial index
    % coincide (== os)
    for os = 1:numel(TTM{sbj}.tstidx)
        % for condition cc=1: sim HRF added
        for cc=1%:2
            % channel indices that have or dont have gt HRF
            idxChHrf = lstHrfAdd{sbj}(:,1);
            idxChNoHrf = arrayfun(@(x) find(lstLongAct{sbj}==x,1),squeeze(lstHrfAdd{sbj}(:,1)));
            % number of available channels
            nHrf = size(FMdc{sbj,os}(:,:,idxChHrf,:,cc));
            nNoHrf = size(FMdc{sbj,os}(:,:,idxChNoHrf,:,cc));
            % extract and append crossvalidated features (from testing trial), new dimension is F x C x I,
            % where F: # of Features, C: # Number of Chromophores, I: # of all
            % trials (epochs*channels)
            F_Raw_Hrf = cat(3, F_Raw_Hrf, FMdc{sbj,os}(:,:,idxChHrf,os,cc));
            F_Raw_NoHrf = cat(3, F_Raw_NoHrf, FMdc{sbj,os}(:,:,idxChNoHrf,os,cc));
            F_SS_Hrf = cat(3, F_SS_Hrf, FMss{sbj,os}(:,:,idxChHrf,os,cc));
            F_SS_NoHrf = cat(3, F_SS_NoHrf,FMss{sbj,os}(:,:,idxChNoHrf,os,cc));
            F_CCA_Hrf = cat(3, F_CCA_Hrf, FMcca{sbj,os}(:,:,idxChHrf,os,cc));
            F_CCA_NoHrf = cat(3, F_CCA_NoHrf, FMcca{sbj,os}(:,:,idxChNoHrf,os,cc));
        end
    end
end

%% Paired T-Tests
for ff = 1:9
    for ch=1:3
        [h_co(ff,ch,1),p_co(ff,ch,1)]= ttest(squeeze(F_Raw_Hrf(ff,ch,:)),squeeze(F_SS_Hrf(ff,ch,:)));
        [h_co(ff,ch,2),p_co(ff,ch,2)]= ttest(squeeze(F_Raw_Hrf(ff,ch,:)),squeeze(F_CCA_Hrf(ff,ch,:)));
        [h_co(ff,ch,3),p_co(ff,ch,3)]= ttest(squeeze(F_SS_Hrf(ff,ch,:)),squeeze(F_CCA_Hrf(ff,ch,:)));
    end
end

%% (Box)Plot results (normal metrics)
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
            boxplot([squeeze(F_Raw_Hrf(ff,ch,:)), squeeze(F_SS_Hrf(ff,ch,:)), squeeze(F_CCA_Hrf(ff,ch,:))], 'labels', labels)
            H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,ch,1:3)));
        else
            % without cca
            boxplot([squeeze(F_Raw_Hrf(ff,ch,:)), squeeze(F_SS_Hrf(ff,ch,:))], 'labels', labels(1:2))
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
                boxplot([abs(squeeze(F_Raw_Hrf(ff,ch,:))-FVgt.x(ff,ch)), abs(squeeze(F_SS_Hrf(ff,ch,:))-FVgt.x(ff,ch)), abs(squeeze(F_CCA_Hrf(ff,ch,:))-FVgt.x(ff,ch))], 'labels', labels)
                H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,ch,1:3)));
            else
                boxplot([abs(squeeze(F_Raw_Hrf(ff,ch,:))-FVgt.x(ff,ch)), abs(squeeze(F_SS_Hrf(ff,ch,:))-FVgt.x(ff,ch))], 'labels', labels(1:2))
                H=sigstar({[1,2]},squeeze(p_co(ff,ch,1)));
            end
            title(['ERR ' FMclab{ff} chrom{ch}])
        else
            if flag_plotCCA
                boxplot([squeeze(F_Raw_Hrf(ff,ch,:)), squeeze(F_SS_Hrf(ff,ch,:)), squeeze(F_CCA_Hrf(ff,ch,:))], 'labels', labels)
                H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,ch,1:3)));
            else
                boxplot([squeeze(F_Raw_Hrf(ff,ch,:)), squeeze(F_SS_Hrf(ff,ch,:))], 'labels', labels(1:2))
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

