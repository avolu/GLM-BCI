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


%% Sort through results and append
F_Raw_Hrf=[];
F_Raw_NoHrf=[];
F_SS_Hrf=[];
F_SS_NoHrf=[];
F_CCA_Hrf=[];
F_CCA_NoHrf=[];
for sbj = 1:numel(sbjfolder)
    for os = 1:numel(TTM{sbj}.tstidx)
        % channel indices that have or dont have gt HRF
        idxChHrf = lstHrfAdd{sbj}(:,1);
        idxChNoHrf = arrayfun(@(x) find(lstLongAct{sbj}==x,1),squeeze(lstHrfAdd{sbj}(:,1)));
        % number of available channels
        sHrf = size(FMdc{sbj}(:,:,idxChHrf,:));
        sNoHrf = size(FMdc{sbj}(:,:,idxChNoHrf,:));
        % extract and append crossvalidated features (from testing trial), new dimension is F x C x I,
        % where F: # of Features, C: # Number of Chromophores, I: # of all
        % trials (epochs*channels)
        F_Raw_Hrf = cat(3, F_Raw_Hrf, reshape(FMdc{sbj,os}(:,:,idxChHrf,:),numel(clab),3,sHrf(3)*sHrf(4)));
        F_Raw_NoHrf = cat(3, F_Raw_NoHrf, reshape(FMdc{sbj,os}(:,:,idxChNoHrf,:),numel(clab),3,sNoHrf(3)*sNoHrf(4)));
        F_SS_Hrf = cat(3, F_SS_Hrf, reshape(FMss{sbj,os}(:,:,idxChHrf,:),numel(clab),3,sHrf(3)*sHrf(4)));
        F_SS_NoHrf = cat(3, F_SS_NoHrf, reshape(FMss{sbj,os}(:,:,idxChNoHrf,:),numel(clab),3,sNoHrf(3)*sNoHrf(4)));
        F_CCA_Hrf = cat(3, F_CCA_Hrf, reshape(FMcca{sbj,os}(:,:,idxChHrf,:),numel(clab),3,sHrf(3)*sHrf(4)));
        F_CCA_NoHrf = cat(3, F_CCA_NoHrf, reshape(FMcca{sbj,os}(:,:,idxChNoHrf,:),numel(clab),3,sNoHrf(3)*sNoHrf(4)));
    end
end

%% Paired T-Tests
for ff = 1:9
    for cc=1:3
        [h_co(ff,cc,1),p_co(ff,cc,1)]= ttest(squeeze(F_Raw_Hrf(ff,cc,:)),squeeze(F_SS_Hrf(ff,cc,:)));
        [h_co(ff,cc,2),p_co(ff,cc,2)]= ttest(squeeze(F_Raw_Hrf(ff,cc,:)),squeeze(F_CCA_Hrf(ff,cc,:)));
        [h_co(ff,cc,3),p_co(ff,cc,3)]= ttest(squeeze(F_SS_Hrf(ff,cc,:)),squeeze(F_CCA_Hrf(ff,cc,:)));
    end
end

%% (Box)Plot results (normal metrics)
figure
labels = {'No GLM', 'GLM SS', 'GLM tCCA'};
chrom = {' HbO', ' HbR'};
% for all features
for ff=1:9
    % for both chromophores
    for cc=1:2
        subplot(2,9,(cc-1)*9+ff)
        xtickangle(35)
        hold on
        %% boxplots
        % with cca
        if flag_plotCCA
            boxplot([squeeze(F_Raw_Hrf(ff,cc,:)), squeeze(F_SS_Hrf(ff,cc,:)), squeeze(F_CCA_Hrf(ff,cc,:))], 'labels', labels)
            H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,cc,1:3)));
        else
            % without cca
            boxplot([squeeze(F_Raw_Hrf(ff,cc,:)), squeeze(F_SS_Hrf(ff,cc,:))], 'labels', labels(1:2))
            H=sigstar({[1,2]},squeeze(p_co(ff,cc,1)));
        end
        
        % ground truth
        if ff<8
            hAx=gca;                                   % retrieve the axes handle
            %xtk=hAx.XTick;                             % and the xtick values to plot() at...
            xtk = [0.5, 1, 2, 3, 3.5];
            hold on
            if ff==5
                FVgt.x(ff,cc) = 6; % set time to peak to 6 seconds (due to gt hrf plateau)
            end
            hL=plot(xtk, ones(numel(xtk))*FVgt.x(ff,cc),'--g');
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
        
        title([clab{ff} chrom{cc}])
    end
end


%% (Box)Plot results (metric errors)
figure
labels = {'No GLM', 'GLM SS', 'GLM tCCA'};
chrom = {' HbO', ' HbR'};
% for all features
for ff=1:9
    % for both chromophores
    for cc=1:2
        subplot(2,9,(cc-1)*9+ff)
        xtickangle(35)
        hold on
        %% boxplots
        if ff==5
            FVgt.x(ff,cc) = 6; % set time to peak to 6 seconds (due to gt hrf plateau)
        end
        
        if ff<8
            % without GLM, with GLM+SS, with GLM+CCA
            if flag_plotCCA
                boxplot([abs(squeeze(F_Raw_Hrf(ff,cc,:))-FVgt.x(ff,cc)), abs(squeeze(F_SS_Hrf(ff,cc,:))-FVgt.x(ff,cc)), abs(squeeze(F_CCA_Hrf(ff,cc,:))-FVgt.x(ff,cc))], 'labels', labels)
                H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,cc,1:3)));
            else
                boxplot([abs(squeeze(F_Raw_Hrf(ff,cc,:))-FVgt.x(ff,cc)), abs(squeeze(F_SS_Hrf(ff,cc,:))-FVgt.x(ff,cc))], 'labels', labels(1:2))
                H=sigstar({[1,2]},squeeze(p_co(ff,cc,1)));
            end
            title(['ERR ' clab{ff} chrom{cc}])
        else
            if flag_plotCCA
                boxplot([squeeze(F_Raw_Hrf(ff,cc,:)), squeeze(F_SS_Hrf(ff,cc,:)), squeeze(F_CCA_Hrf(ff,cc,:))], 'labels', labels)
                H=sigstar({[1,2],[1,3],[2,3]},squeeze(p_co(ff,cc,1:3)));
            else
                boxplot([squeeze(F_Raw_Hrf(ff,cc,:)), squeeze(F_SS_Hrf(ff,cc,:))], 'labels', labels(1:2))
                H=sigstar({[1,2]},squeeze(p_co(ff,cc,1)));
            end
            
            title([clab{ff} chrom{cc}])
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


