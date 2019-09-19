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

load([path.save '\FV_results.mat'])

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

% for conventional NO GLM features
gg=1;
% select features
% [ min, max, p2p, avg, t2p, slope w1, slope w2]
fsel = [3,4,6,7];
% for all subjects
for sbj=1:numel(TTM)
    % for all trials
    for tt = 1:numel(TTM{sbj}.tstidx)
        xTr{gg,sbj,tt} =[];
        xTst{gg,sbj,tt}=[];
        yTr{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tnridx(tt,:)));
        yTst{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tstidx(tt)));
        for cc=1:2
            % train data  (from GLM with trained HRF regressor on seen training data)
            % append features for hbo and hbr and all channels without SS
            fvbuf = [];
            fvbuf = squeeze(FMdc{sbj}(fsel,1:2,lstLongAct{sbj},TTM{sbj}.tnridx(tt,:),cc));
            xTr{gg,sbj,tt} = [xTr{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2)*size(fvbuf,3),numel(TTM{sbj}.tnridx(tt,:)))];
            % generate label vector
            yTr{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tnridx(tt,:))+1:cc*numel(TTM{sbj}.tnridx(tt,:)))=1;
            % test data (from GLM with trained HRF regressor on unseen data)
            % append features for hbo and hbr and all channels without SS
            fvbuf = [];
            fvbuf = squeeze(FMdc{sbj}(fsel,1:2,lstLongAct{sbj},TTM{sbj}.tstidx(tt),cc,rr));
            xTst{gg,sbj,tt} = [xTst{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2)*size(fvbuf,3),numel(TTM{sbj}.tstidx(tt)))];
            % generate label vector
            yTst{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tstidx(tt))+1:cc*numel(TTM{sbj}.tstidx(tt)))=1;
        end
    end
end


% for both GLM methods
FW = {[], FWss, FWcca};
for gg = 2:3
    % for all subjects
    for sbj=1:numel(TTM)
        % for all trials
        for tt = 1:numel(TTM{sbj}.tstidx)
            xTr{gg,sbj,tt} =[];
            xTst{gg,sbj,tt}=[];
            yTr{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tnridx(tt,:)));
            yTst{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tstidx(tt)));
            for cc=1:2
                % train data  (from GLM with trained HRF regressor on seen training data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf = squeeze(FW{gg}{sbj,tt}(:,:,lstLongAct{sbj},TTM{sbj}.tnridx(tt,:),cc,rr));
                xTr{gg,sbj,tt} = [xTr{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2),numel(TTM{sbj}.tnridx(tt,:)))];
                % generate label vector
                yTr{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tnridx(tt,:))+1:cc*numel(TTM{sbj}.tnridx(tt,:)))=1;
                % test data (from GLM with trained HRF regressor on unseen data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf = squeeze(FW{gg}{sbj,tt}(:,:,lstLongAct{sbj},TTM{sbj}.tstidx(tt),cc,rr));
                xTst{gg,sbj,tt} = [xTst{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2),numel(TTM{sbj}.tstidx(tt)))];
                % generate label vector
                yTst{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tstidx(tt))+1:cc*numel(TTM{sbj}.tstidx(tt)))=1;
            end
        end
    end
end




%% CROSSVALIDATION using rLDA as classifier
% for all methods (1> NO GLM, 2> GLM SS, 3> GLM CCA)
for gg = 1:3
    % for all subjects
    for sbj=1:numel(TTM)
        % for all splits
        for tt=1:numel(TTM{sbj}.tstidx)
            
            %% training of rLDA
            C = train_RLDAshrink(xTr{gg,sbj,tt}, yTr{gg,sbj,tt});
            
            %% testing of rLDA
            fv.x = xTst{gg,sbj,tt};
            fv.y = yTst{gg,sbj,tt};
            out = applyClassifier(fv, C);
            % loss function
            loss{gg,sbj}(tt,:,:)=loss_classwiseNormalized(fv.y, out, size(fv.y));
        end
        lossAvg(gg,sbj) = mean(loss{gg,sbj}(:));
    end
end
accuracyGLM = 1-lossAvg;

meanaccs = mean(accuracyGLM,2);


figure
bar(accuracyGLM)
hold on
plot([.5 3.5], [0.5 0.5], '--k')
set(gca,'xtickLabel',{...
    ['no GLM (' num2str(100*meanaccs(1), '%2.1f') '%)'], ...
    ['GLM SS (' num2str(100*meanaccs(2),'%2.1f') '%)'], ...
    ['GLM tCCA (' num2str(100*meanaccs(3),'%2.1f') '%)']})
ylabel('mean accuracy / subject')
title('subject avg classification accuracies for each method')
