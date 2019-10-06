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
load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindrift_hrf_amp50_20soffs'])


% load and init BBCI toolbox
% bbci toolbox paths
if malexflag
    %Meryem
    paths.bbciDir = 'C:\Users\mayucel\Documents\PROJECTS\CODES\bbci_public-master';
    paths.bbciDataDir = 'C:\Users\mayucel\Documents\PROJECTS\CODES\bbci_public-master\bbci_data';
    paths.bbciTmpDir = 'C:\Users\mayucel\Documents\PROJECTS\CODES\bbci_public-master\bbci_data\tmp';
    addpath(genpath(paths.bbciDir))
    cd(paths.bbciDir);
    startup_bbci_toolbox('DataDir', paths.bbciDataDir, 'TmpDir',paths.bbciTmpDir);
else
    % Alex
    paths.bbciDir = 'D:\Office\Archive Office\Toolboxes - Code Libraries\Matlab\BBCI\';
    paths.bbciDataDir = 'D:\Datasets\bbci_data';
    paths.bbciTmpDir = 'D:\Datasets\bbci_data\tmp\';
    addpath(genpath(paths.bbciDir))
    cd(paths.bbciDir);
    startup_bbci_toolbox('DataDir', paths.bbciDataDir, 'TmpDir',paths.bbciTmpDir);
end

% use hrf STIM regressor weights as features (not the REST regressor
% weights, as they are useless here)
% and transform to bbci data structure
rr = 1;
epo.className = {'STIM', 'REST'};
epo.clab = FMclab;

%% for conventional features
% select features
flab = {'min |', 'max |', 'p2p |', 'avg |', 't2p |', 'slope w1 |', 'slope w2 |'};
fsel = [4];
% for all subjects
FW = {FMdc', FMss, FMcca};
for gg = 1:3
    for sbj=1:numel(TTM)
        % for all trials
        for tt = 1:numel(TTM{sbj}.tstidx)
            if gg == 1
                cvidx = 1;
            else
                cvidx = tt;
            end
            
            xTrF{gg,sbj,tt} =[];
            xTstF{gg,sbj,tt}=[];
            yTrF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tnridx(tt,:)));
            yTstF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tstidx(tt)));
            for cc=1:2
                % train data  (from GLM with trained HRF regressor on seen training data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf = FW{gg}{sbj,cvidx}(fsel,1:2,lstLongAct{sbj},TTM{sbj}.tnridx(tt,:),cc);
                xTrF{gg,sbj,tt} = [xTrF{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2)*size(fvbuf,3),numel(TTM{sbj}.tnridx(tt,:)))];
                % generate label vector
                yTrF{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tnridx(tt,:))+1:cc*numel(TTM{sbj}.tnridx(tt,:)))=1;
                % test data (from GLM with trained HRF regressor on unseen data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf = squeeze(FW{gg}{sbj,cvidx}(fsel,1:2,lstLongAct{sbj},TTM{sbj}.tstidx(tt),cc,rr));
                xTstF{gg,sbj,tt} = [xTstF{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2)*size(fvbuf,3),numel(TTM{sbj}.tstidx(tt)))];
                % generate label vector
                yTstF{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tstidx(tt))+1:cc*numel(TTM{sbj}.tstidx(tt)))=1;
            end
        end
    end
end


%% for weight features from GLM methods
FW = {[], FWss, FWcca};
for gg = 2:3
    % for all subjects
    for sbj=1:numel(TTM)
        % for all trials
        for tt = 1:numel(TTM{sbj}.tstidx)
            xTrW{gg,sbj,tt} =[];
            xTstW{gg,sbj,tt}=[];
            yTrW{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tnridx(tt,:)));
            yTstW{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tstidx(tt)));
            for cc=1:2
                % train data  (from GLM with trained HRF regressor on seen training data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf = squeeze(FW{gg}{sbj,tt}(:,:,lstLongAct{sbj},TTM{sbj}.tnridx(tt,:),cc,rr));
                xTrW{gg,sbj,tt} = [xTrW{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2),numel(TTM{sbj}.tnridx(tt,:)))];
                % generate label vector
                yTrW{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tnridx(tt,:))+1:cc*numel(TTM{sbj}.tnridx(tt,:)))=1;
                % test data (from GLM with trained HRF regressor on unseen data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf = squeeze(FW{gg}{sbj,tt}(:,:,lstLongAct{sbj},TTM{sbj}.tstidx(tt),cc,rr));
                xTstW{gg,sbj,tt} = [xTstW{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2),numel(TTM{sbj}.tstidx(tt)))];
                % generate label vector
                yTstW{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tstidx(tt))+1:cc*numel(TTM{sbj}.tstidx(tt)))=1;
            end
        end
    end
end

%% CROSSVALIDATION using rLDA as classifier and normal features
% for all methods (1> NO GLM, 2> GLM SS, 3> GLM CCA)
for gg = 1:3
    % for all subjects
    for sbj=1:numel(TTM)
        % for all splits
        for tt=1:numel(TTM{sbj}.tstidx)
            
            %% training of rLDA
            C = train_RLDAshrink(xTrF{gg,sbj,tt}, yTrF{gg,sbj,tt});
            
            %% testing of rLDA
            fv.x = xTstF{gg,sbj,tt};
            fv.y = yTstF{gg,sbj,tt};
            out = applyClassifier(fv, C);
            % loss function
            lossF{gg,sbj}(tt,:,:)=loss_classwiseNormalized(fv.y, out, size(fv.y));
        end
        lossAvgF(gg,sbj) = mean(lossF{gg,sbj}(:));
    end
end
accuracyGLMF = 1-lossAvgF;
meanaccsF = mean(accuracyGLMF,2);


figure
bar(accuracyGLMF)
hold on
plot([.5 3.5], [0.5 0.5], '--k')
set(gca,'xtickLabel',{...
    ['no GLM (' num2str(100*meanaccsF(1), '%2.1f') '%)'], ...
    ['GLM SS (' num2str(100*meanaccsF(2),'%2.1f') '%)'], ...
    ['GLM tCCA (' num2str(100*meanaccsF(3),'%2.1f') '%)']})
ylabel('mean accuracy / subject')
title(['Normal Features ' flab{fsel}])



%% CROSSVALIDATION using rLDA as classifier and weight features
% for GLM methods ( 2> GLM SS, 3> GLM CCA)

for gg = 2:3
    % for all subjects
    for sbj=1:numel(TTM)
        % for all splits
        for tt=1:numel(TTM{sbj}.tstidx)
            
            %% training of rLDA
            C = train_RLDAshrink(xTrW{gg,sbj,tt}, yTrW{gg,sbj,tt});
            
            %% testing of rLDA
            fv.x = xTstW{gg,sbj,tt};
            fv.y = yTstW{gg,sbj,tt};
            out = applyClassifier(fv, C);
            % loss function
            lossW{gg,sbj}(tt,:,:)=loss_classwiseNormalized(fv.y, out, size(fv.y));
        end
        lossAvgW(gg,sbj) = mean(lossW{gg,sbj}(:));
    end
end
accuracyGLMW = 1-lossAvgW;
accuracyGLMW(1,:)=0;
meanaccsW = mean(accuracyGLMW,2);


figure
bar(accuracyGLMW)
hold on
plot([.5 3.5], [0.5 0.5], '--k')
set(gca,'xtickLabel',{...
    '', ...
    ['GLM SS (' num2str(100*meanaccsW(2),'%2.1f') '%)'], ...
    ['GLM tCCA (' num2str(100*meanaccsW(3),'%2.1f') '%)']})
ylabel('mean accuracy / subject')
title('Weight Features')
