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

% include tcca results in plots?
flag_plotCCA = true;

%% Get HRF features from all augmented channels to compare against  ground truth
%Sort through results and append
W_SS_Hrf=[];
W_SS_NoHrf=[];
W_CCA_Hrf=[];
W_CCA_NoHrf=[];
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
            nHrf = size(FMdc{sbj,os}(:,:,idxChHrf,:,cc));
            nNoHrf = size(FMdc{sbj,os}(:,:,idxChNoHrf,:,cc));
            % extract and append crossvalidated weights (from testing trial), new dimension is 
            % CHROMOPHORES x Concatenated weights of active CHANNELS (trials) x CONDITION HRF,
            W_SS_Hrf = cat(2, W_SS_Hrf, squeeze(FWss{sbj,os}(1,:,idxChHrf,os,:)));
            W_SS_NoHrf = cat(2, W_SS_NoHrf, squeeze(FWss{sbj,os}(1,:,idxChNoHrf,os,:)));
            W_CCA_Hrf = cat(2, W_CCA_Hrf, squeeze(FWcca{sbj,os}(1,:,idxChHrf,os,:)));
            W_CCA_NoHrf = cat(2, W_CCA_NoHrf, squeeze(FWcca{sbj,os}(1,:,idxChNoHrf,os,:)));

    end
end

%% (Scatter)Plot results (weights)
% plots the weights 
figure
% hrf regressor for condition
hrfreg = 1;
hrfreglab = {'STIM', 'NO STIM'};
% hbo(1)/hbr(2)
hb = 1;
hblab = {'HbO', 'HbR'};
hbcol ={'or', 'ob'};
W = {W_SS_Hrf, W_SS_NoHrf, W_CCA_Hrf, W_CCA_NoHrf};
wlab = {'GLM SS', '', 'GLM CCA', ''};
for ww = [1,3]
    for hb=1:2
        subplot(2,2,ww+hb-1)
        scatter(abs(W{ww}(hb,:,hrfreg)),abs(W{ww+1}(hb,:,hrfreg)), hbcol{hb})
        hold on
        scatter(nanmean(abs(W{ww}(hb,:,hrfreg))),nanmean(abs(W{ww+1}(hb,:,hrfreg))), 'xk', 'Linewidth', 2)
        minmax = [min([abs(W{ww}(hb,:,hrfreg)) abs(W{ww+1}(hb,:,hrfreg))]) max([abs(W{ww}(hb,:,hrfreg)), abs(W{ww+1}(hb,:,hrfreg))])];
        plot(minmax, minmax,'k')
        xlabel('Regressor weight Condition STIM')
        ylabel('Regressor weight Condition REST')
        title(['GLM ' wlab{ww} ' hrf ' hrfreglab{hrfreg} ' regressor, ' hblab{hb}])
        axis tight
        grid on
    end
end

figure
histogram(W_SS_Hrf(1,:,1))
hold on
histogram(W_SS_NoHrf(1,:,1))

figure
% hrf regressor for condition
hrfreg = 2;
hrfreglab = {'STIM', 'NO STIM'};
% hbo(1)/hbr(2)
hb = 1;
hblab = {'HbO', 'HbR'};
hbcol ={'or', 'ob'};
W = {W_SS_Hrf, W_SS_NoHrf, W_CCA_Hrf, W_CCA_NoHrf};
wlab = {'GLM SS', '', 'GLM CCA', ''};
for ww = [1,3]
    for hb=1:2
        subplot(2,2,ww+hb-1)
        scatter(abs(W{ww}(hb,:,hrfreg)),abs(W{ww+1}(hb,:,hrfreg)), hbcol{hb})
        hold on
        scatter(nanmean(abs(W{ww}(hb,:,hrfreg))),nanmean(abs(W{ww+1}(hb,:,hrfreg))), 'xk', 'Linewidth', 2)
        minmax = [min([abs(W{ww}(hb,:,hrfreg)) abs(W{ww+1}(hb,:,hrfreg))]) max([abs(W{ww}(hb,:,hrfreg)), abs(W{ww+1}(hb,:,hrfreg))])];
        plot(minmax, minmax,'k')
        xlabel('Regressor weight Condition STIM')
        ylabel('Regressor weight Condition REST')
        title(['GLM ' wlab{ww} ' hrf ' hrfreglab{hrfreg} ' regressor, ' hblab{hb}])
        axis tight
        grid on
    end
end