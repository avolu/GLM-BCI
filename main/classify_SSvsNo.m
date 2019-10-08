clear all;

malexflag = 1; % user flag
if malexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\GLM_BCI_PAPER\RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\GLM_BCI_PAPER\PROCESSED_DATA'; % save directory
    
    %Meryem Laptop
    %     path.code = 'C:\Users\m\Documents\GitHub\GLM-BCI'; addpath(genpath(path.code)); % code directory
    %     path.dir = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    %     path.save = 'C:\Users\m\Documents\tCCA_GLM_PAPER\FB_RESTING_DATA'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\GLM_BCI_PAPER\RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\GLM_BCI_PAPER\PROCESSED_DATA'; % save directory
end

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

%%choose HRF level
hrflab = {'HRF 100%', 'HRF 50%'};
hh = 1;


disp(['running for ' hrflab{hh} '...'])
%% load data
switch hh
    case 1 % 100%
        %load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindrift_hrf_amp100_20soffs.mat'])
        %         load([path.save '\FV_results_SSvsNo_ldrift1_resid0_tccaIndiv_hrf_amp100_20soffs.mat'])
        load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp100_20soffs.mat'])
        load([path.save '\test.mat'])
    case 2 % 50%
        %load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindrift_hrf_amp50_20soffs.mat'])
        %         load([path.save '\FV_results_SSvsNo_ldrift1_resid0_tccaIndiv_hrf_amp50_20soffs.mat'])
        load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp50_20soffs.mat'])
end


% use hrf STIM regressor weights as features (not the REST regressor
% weights, as they are useless here)
% and transform to bbci data structure
rr = 1;
epo.className = {'STIM', 'REST'};
epo.clab = FMclab;

%% sbj list
sbjl = [1:3 5:14];
%% chromophores (HbO / HbR)
chrom = [1 2];
%% channel selection
cfact=4;

%% get weight features from GLM method
%dimensionality of FWss:
%Feature type | chromophore | channels | trials | condition | regressor
FW = FWss;
% for all subjects
gg=1;
for sbj=sbjl
    %ch selection
    chsel = lstLongAct{sbj};
    chsel = chsel(1:floor(numel(chsel)/cfact));
    % for all trials
    for tt = 1:numel(TTM{sbj}.tstidx)
        xTrF{gg,sbj,tt} =[];
        xTstF{gg,sbj,tt}=[];
        yTrF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tnridx(tt,:)));
        yTstF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tstidx(tt)));
        % conditions
        for cc=1:2
            % train data  (from GLM with trained HRF regressor on seen training data)
            % append features for  chromophores(hbo and hbr) and all channels without SS
            fvbuf = [];
            fvbuf = squeeze(FW{sbj,tt}(:,chrom,chsel,TTM{sbj}.tnridx(tt,:),cc,rr));
            if numel(chrom) == 1
                xTrF{gg,sbj,tt} = [xTrF{gg,sbj,tt} fvbuf];
            else
                xTrF{gg,sbj,tt} = [xTrF{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2),numel(TTM{sbj}.tnridx(tt,:)))];
            end
            % generate label vector
            yTrF{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tnridx(tt,:))+1:cc*numel(TTM{sbj}.tnridx(tt,:)))=1;
            % test data (from GLM with trained HRF regressor on unseen data)
            % append features for hbo and hbr and all channels without SS
            fvbuf = [];
            fvbuf = squeeze(FW{sbj,tt}(:,chrom,chsel,TTM{sbj}.tstidx(tt),cc,rr));
            if numel(chrom) == 1
                xTstF{gg,sbj,tt} = [xTstF{gg,sbj,tt} fvbuf];
            else
                xTstF{gg,sbj,tt} = [xTstF{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2),numel(TTM{sbj}.tstidx(tt)))];
            end
            % generate label vector
            yTstF{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tstidx(tt))+1:cc*numel(TTM{sbj}.tstidx(tt)))=1;
        end
        %feature vector dimensionality
        nfeat(gg,sbj,tt) = size(xTrF{gg,sbj,tt},1);
    end
end
glab = {'W | GLM'};

%% for conventional features
% select features
flab = {'min ', 'max ', 'p2p ', 'avg ', 't2p ', 'slope ', 'slope w2 '};
mlab = {'no GLM', 'GLM'};
fsel = {[], [1], [1], [2], [2], [3], [3], [4], [4], [5], [5], [6], [6],...
    [3,4], [3,4], [3,5], [3,5], [3,6], [3,6], [4,5], [4,5], [4,6], [4,6], ...
    [5,6],[5,6]};
%dimensionality of FW:
%Feature type | chromophore | channels | trials | condition | regressor
FW = {FMdc', FMss};
% for all features
for gg = 2:numel(fsel)
    % for all subjects
    for sbj=sbjl
        %ch selection
    chsel = lstLongAct{sbj};
    chsel = chsel(1:floor(numel(chsel)/cfact));
        % for all trials
        for tt = 1:numel(TTM{sbj}.tstidx)
            if mod(gg,2)+1 == 1
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
                fvbuf = FW{mod(gg,2)+1}{sbj,cvidx}(fsel{gg},chrom,chsel,TTM{sbj}.tnridx(tt,:),cc);
                xTrF{gg,sbj,tt} = [xTrF{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2)*size(fvbuf,3),numel(TTM{sbj}.tnridx(tt,:)))];
                % generate label vector
                yTrF{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tnridx(tt,:))+1:cc*numel(TTM{sbj}.tnridx(tt,:)))=1;
                % test data (from GLM with trained HRF regressor on unseen data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf =FW{mod(gg,2)+1}{sbj,cvidx}(fsel{gg},chrom,chsel,TTM{sbj}.tstidx(tt),cc,rr);
                xTstF{gg,sbj,tt} = [xTstF{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2)*size(fvbuf,3),numel(TTM{sbj}.tstidx(tt)))];
                % generate label vector
                yTstF{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tstidx(tt))+1:cc*numel(TTM{sbj}.tstidx(tt)))=1;
                % feature vector dimensionality
                nfeat(gg,sbj,tt) = size(xTrF{gg,sbj,tt},1);
            end
        end
        glab{gg}= [join([join(flab(fsel{gg})) mlab{mod(gg,2)+1}])];
    end
end

%% CROSSVALIDATION using rLDA as classifier and all features
% for all methods/feature types (1> NO GLM, 2> GLM SS, 3> GLM CCA)
for gg = 1:size(xTrF,1)
    disp(['CV for all subjects and feature set ' num2str(gg) '...'])
    % for all subjects
    for sbj=sbjl
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
accuracyGLMF{hh} = 1-lossAvgF;
meanaccsF{hh} = mean(accuracyGLMF{hh}(:,sbjl),2);

%% Plot Classification results
glmacc = round(accuracyGLMF{hh}(1:2:end,sbjl)*1000)/10;
noglmacc(1,1:numel(sbjl)) = 0;
noglmacc = round([noglmacc; accuracyGLMF{hh}(2:2:end,sbjl)]*1000)/10;
% xticklabels (class accuracies)
xtckglm =  round(meanaccsF{hh}(1:2:end,:)*1000)/10;
xtcknoglm(1) = 0;
xtcknoglm = round([xtcknoglm; meanaccsF{hh}(2:2:end)]*1000)/10;
% used feature labels
fidx = fsel(1:2:end);
for i=1:numel(fidx)
    uflab{i} = string(join(flab(fidx{i})))
end
uflab{1} = '\beta GLM ';

figure
subplot(2,1,1)
hold on
bar(glmacc,'EdgeColor','none','BarWidth',1)
plot([.5 13.5], [50 50], '--k')
ylabel('sbj avg accuracy / %')
xlim([0.5 13.5])
xticks(1:13)
xticklabels(xtckglm)
% create upper x axis for feature type labels
ax1 = gca;
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlim([0.5, 13.5])
xticks(1:13)
yticks([])
xticklabels(uflab)
subplot(2,1,2)
hold on
xlim([0.5, 13.5])
xticks(1:13)
bar(noglmacc,'EdgeColor','none','BarWidth',1)
plot([.5 13.5], [50 50], '--k')
ylabel('sbj avg accuracy / %')
xticks(1:13)
xticklabels(xtcknoglm)

% Perform paired tt-tests
[h{hh},p{hh}]= ttest(accuracyGLMF{hh}(3:2:end,sbjl), accuracyGLMF{hh}(2:2:end,sbjl));



