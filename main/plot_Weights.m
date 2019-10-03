clear all;

malexflag = 1; % user flag
if malexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\GLM_BCI_PAPER\RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\GLM_BCI_PAPER\PROCESSED_DATA\sbj_opt_tcca'; % save directory
    
    %Meryem Laptop
    %     path.code = 'C:\Users\m\Documents\GitHub\GLM-BCI'; addpath(genpath(path.code)); % code directory
    %     path.dir = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    %     path.save = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\GLM_BCI_PAPER\RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\GLM_BCI_PAPER\PROCESSED_DATA\sbj_opt_tcca\'; % save directory
end
%% load data
%load([path.save '\FV_results_std_nReg2_ldrift1_resid0_tccap1.mat'])
%load([path.save '\FV_results_std_nReg2_ldrift1_resid0_tccap1_20soffs.mat'])
%load([path.save '\FV_results_std_nReg2_ldrift1_resid1_tccap1.mat'])
%load([path.save '\FV_results_std_nReg3_ldrift0_resid1_tccap1.mat'])
%load([path.save '\FV_results_std_nReg2_ldrift0_resid0_tccap1_20soffs.mat'])
%load([path.save '\FV_results_std_nReg2_ldrift1_resid1_tccap2'])
%load([path.save '\FV_results_std_nReg2_ldrift1_resid0_tccap2'])
% load([path.save '\FV_results_std_nReg2_ldrift1_resid0_tccap2_20soffs.mat']) % good
%load([path.save '\FV_results_std_nReg2_ldrift1_resid0_tccap1.mat'])
%load([path.save '\FV_results_std_nReg2_ldrift0_resid0_tccap2_20soffs.mat'])
%load([path.save '\FV_results_std_nReg2_ldrift1_resid0_tccap2_hrf_amp50_20soffs.mat'])

% individualized parameters
load([path.save '\FV_results_std_nReg2_ldrift1_resid0_tccaIndiv_hrf_amp50_20soffs.mat'])


%% Get HRF features from all augmented channels to compare against  ground truth
%Sort through results and append
W_SS = cell(2,2);
W_CCA = cell(2,2);
for sbj = 1:numel(TTM)
    % only look at the crossvalidated test results. These are stored where cell and trial index
    % coincide (== os)
    for os = 1:numel(TTM{sbj}.tstidx)
        % channel indices that have or dont have gt HRF
        idxChHrf = lstHrfAdd{sbj}(:,1);
        idxChNoHrf = setdiff(lstLongAct{sbj},squeeze(lstHrfAdd{sbj}(:,1)));
        if size(idxChHrf,1) > size(idxChNoHrf,1)
            idxChHrf = idxChHrf(1:size(idxChNoHrf,1));
        else
            idxChNoHrf = idxChNoHrf(1:size(idxChHrf,1));
        end
        % number of available channels for HRF added
        cc=1;
        nHrf = size(FWss{sbj,os}(:,:,idxChHrf,:,cc));
        nNoHrf = size(FWss{sbj,os}(:,:,idxChNoHrf,:,cc));
        % extract and append crossvalidated weights (from testing trial), new dimension is
        % CHROMOPHORES x Concatenated weights of active CHANNELS (trials)
        for cc=1:2 % stim and resting condition
            for rr=1:2 % stim and resting hrf regressor
                W_SS{cc,rr} = cat(2, W_SS{cc,rr}, squeeze(FWss{sbj,os}(1,:,idxChHrf,os,cc,rr)));
                W_CCA{cc,rr} = cat(2, W_CCA{cc,rr}, squeeze(FWcca{sbj,os}(1,:,idxChHrf,os,cc,rr)));
            end
        end
    end
end

%% (Scatter)Plot results (weights)
% hrf regressor for condition
hrfreglab = {'STIM', 'REST'};
hblab = {'HbO', 'HbR'};
hbcol ={'or', 'ob'};
W = {W_SS, W_CCA};
wlab = {'GLM SS', 'GLM CCA'};
% Plot both methods
figure
for ww = 1:2
    % for HRF stim regressor
    for rr = 1:2
        % for both chromophores
        for hb=1:2
            % plot weights for HRF stim regressor for STIM vs REST trials (cc=1 vs 2)
            subplot(2,4,rr+(hb-1)*4+(ww-1)*2)
            scatter(abs(W{ww}{1,rr}(hb,:)),abs(W{ww}{2,rr}(hb,:)), hbcol{hb})
            hold on
            scatter(nanmean(abs(W{ww}{1,rr}(hb,:))),nanmean(abs(W{ww}{2,rr}(hb,:))), 'xk', 'Linewidth', 2)
            minmax = [min([abs(W{ww}{1,rr}(hb,:)) abs(W{ww}{2,rr}(hb,:))]) max([abs(W{ww}{1,rr}(hb,:)) abs(W{ww}{2,rr}(hb,:))])];
            plot(minmax, minmax,'k')
            xlabel('STIM Trials')
            ylabel('REST Trials')
            title([wlab{ww} ' hrf ' hrfreglab{rr} ' regressor, ' hblab{hb}])
            axis tight
            grid on
        end
    end
end


%% Histograms of weights
% hrf regressor for condition
hrfreglab = {'STIM', 'REST'};
hblab = {'HbO', 'HbR'};
hbcol ={'or', 'ob'};
W = {W_SS, W_CCA};
wlab = {'GLM SS', 'GLM CCA'};
% Plot both methods
figure
for ww = 1:2
    % for HRF stim regressor
    for rr = 1:2
        % for both chromophores
        for hb=1:2
            % plot weights for HRF stim regressor for STIM vs REST trials (cc=1 vs 2)
            subplot(2,4,rr+(hb-1)*4+(ww-1)*2)
            histogram(W{ww}{1,rr}(hb,:))
            hold on
            histogram(W{ww}{2,rr}(hb,:))
            legend('STIM Trials','REST Trials')
            title([wlab{ww} ' hrf ' hrfreglab{rr} ' regressor, ' hblab{hb}])
            axis tight
            grid on
        end
    end
end


