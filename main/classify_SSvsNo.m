clear all;

malexflag = 0; % user flag
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
        %load([path.save '\FV_results_SSvsNo_ldrift1_resid0_tccaIndiv_hrf_amp100_20soffs.mat'])
        %load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp100_20soffs.mat'])
%         load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp100_mot_corr0_20soffs.mat'])
%         load([path.save '\chselInfo_FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp100_mot_corr0_20soffs_ttstchsel.mat'])
                load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp100_mot_corr0_20soffs_ttstchsel.mat'])
    case 2 % 50%
        %load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindrift_hrf_amp50_20soffs.mat'])
        %load([path.save '\FV_results_SSvsNo_ldrift1_resid0_tccaIndiv_hrf_amp50_20soffs.mat'])
        %load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp50_20soffs.mat'])
%         load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp50_mot_corr0_20soffs.mat'])
%         load([path.save '\chselInfo_FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp100_mot_corr0_20soffs_ttstchsel.mat'])
                load([path.save '\FV_results_SSvsNo_ldrift1_resid0stlindriftSWAPPED_1_hrf_amp50_mot_corr0_20soffs_ttstchsel.mat'])
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
chselType = 'r2';

switch chselType
    case {'all', 'r2'} % no selection, take all channels (r2 selection in CV down below)
        for sbj=sbjl
            for tt = 1:numel(TTM{sbj}.tstidx)
                for hh = 1:2
                    chsel{sbj,hh,tt} = lstLongAct{sbj};
                end
            end
        end 
    case 'ratio' % (fraction of available channels)
        nchaddratio = 1/4;
        for sbj=sbjl
            chselhrf{sbj} = lstHrfAdd{sbj}(1:ceil(size(lstHrfAdd{sbj},1)*nchaddratio));
            chselnohrf{sbj} = setdiff(lstLongAct{sbj}, lstHrfAdd{sbj});
            chselnohrf{sbj} = chselnohrf{sbj}(1:floor(numel(chselnohrf{sbj})*(1-nchaddratio)));
            chsel{sbj,1,1} = [ chselhrf{sbj} chselnohrf{sbj}' ];
            chsel{sbj,2,1} = chsel{sbj,1,1};
            for tt = 1:numel(TTM{sbj}.tstidx)
                chsel{sbj,1,tt} = chsel{sbj,1,1};
                chsel{sbj,2,tt} = chsel{sbj,2,1};
            end
        end
    case 'hrfChanOnly' % use only channels with known hrf modulation?
        %% channel indices that have or dont have gt HRF
        for sbj=sbjl
            chsel{sbj,1,1} = lstHrfAdd{sbj}(:,1,1);
            chsel{sbj,2,1} = chsel{sbj,1,1};
            for tt = 1:numel(TTM{sbj}.tstidx)
                chsel{sbj,1,tt} = chsel{sbj,1,1};
                chsel{sbj,2,tt} = chsel{sbj,2,1};
            end
        end
    case 'ttest' % perform channel selection based on across condition avg peak t-test
        % number of channels to restrict to
        nchan = 10;
        for sbj=sbjl
            gg=1; % raw data
            %for both chromophores
            for tt = 1:numel(TTM{sbj}.tstidx) % for all tst indices (select channels from train set)
                for hh = 1:2
                    for ch = 1:size(chselInfo{sbj}.BL_RAW,3) % for all channels
                        pkc1 = squeeze(chselInfo{sbj}.PEAK_RAW(TTM{sbj}.tnridx(tt,:),hh,ch,1)); % peak condition 1 (stim)
                        pkc2 = squeeze(chselInfo{sbj}.PEAK_RAW(TTM{sbj}.tnridx(tt,:),hh,ch,2)); % peak condition 2 (rest)
                        [h(hh,ch,tt), p(hh,ch,tt)] = ttest(pkc1, pkc2);
                        if isnan(h(hh,ch,tt))
                            h(hh,ch,tt)=0;
                        end
                    end
                end
                [b,i]=sort(p(1,:,tt),'ascend');
                chsel{sbj,gg,tt} = i(1:nchan);
                %                 chsel{sbj,gg,tt} = find(h(1,:,tt)& h(2,:,tt));
            end
            gg=2; % GLM data
            %for both chromophores
            for tt = 1:numel(TTM{sbj}.tstidx) % for all tst indices (select channels from train set)
                for hh = 1:2
                    for ch = 1:size(chselInfo{sbj}.BL_SS,3) % for all channels
                        pkc1 = squeeze(chselInfo{sbj}.PEAK_SS(tt,hh,ch,TTM{sbj}.tnridx(tt,:),1)); % peak condition 1 (stim)
                        pkc2 = squeeze(chselInfo{sbj}.PEAK_SS(tt,hh,ch,TTM{sbj}.tnridx(tt,:),2)); % peak condition 2 (rest)
                        [h(hh,ch,tt), p(hh,ch,tt)] = ttest(pkc1, pkc2);
                        if isnan(h(hh,ch,tt))
                            h(hh,ch,tt)=0;
                        end
                    end
                end
                [b,i]=sort(p(1,:,tt),'ascend');
                chsel{sbj,gg,tt} = i(1:nchan);
                %                 chsel{sbj,gg,tt} = find(h(1,:,tt)& h(2,:,tt));
            end
        end
        %% summary comparison of selected channels vs channels with added hrf
        for sbj = sbjl
            for gg = 1:2
                for tt = 1:numel(TTM{sbj}.tstidx)
                    chNsel(sbj,gg,tt) = numel(chsel{sbj,gg,tt});
                    chSelOverlap(sbj,gg,tt) = numel(intersect(chsel{sbj,gg,tt}, lstHrfAdd{sbj}))/numel(lstHrfAdd{sbj}(:,1));
                    chSelperCorr(sbj,gg,tt) = numel(intersect(chsel{sbj,gg,tt}, lstHrfAdd{sbj}))/numel(chsel{sbj,gg,tt});
                end
            end
        end
        disp(['Across trial avg overlap selected channels and channels with HRF added [%] - No GLM: ' num2str(squeeze(mean(chSelOverlap(sbjl,1,:),3))')])
        disp(['Overall avg overlap selected channels and channels with HRF added [%] - No GLM: ' num2str(mean(squeeze(mean(chSelOverlap(sbjl,1,:),3))'))])
        disp(['Overall ratio correctly selected channels to overall selected [%]- No GLM: ' num2str(nanmean(squeeze(mean(chSelperCorr(sbjl,1,:),3))'))])
        disp(['Across trial avg overlap selected channels and channels with HRF added [%] - GLM SS: ' num2str(squeeze(mean(chSelOverlap(sbjl,2,:),3))')])
        disp(['Overall avg overlap selected channels and channels with HRF added [%] - GLM SS: ' num2str(mean(squeeze(mean(chSelOverlap(sbjl,2,:),3))'))])
        disp(['Overall ratio correctly selected channels to overall selected [%] - GLM SS: ' num2str(nanmean(squeeze(mean(chSelperCorr(sbjl,2,:),3))'))])
end
clear vars h p



%% get weight features from GLM method
%dimensionality of FWss:
%Feature type | chromophore | channels | trials | condition | regressor
FW = FWss;
% for all subjects
gg=1;
for sbj=sbjl
    % for all trials
    for tt = 1:numel(TTM{sbj}.tstidx)
        xTrF{gg,sbj,tt} =[];
        xTstF{gg,sbj,tt}=[];
        yTrF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tnridx(tt,:)));
        yTstF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tstidx(tt)));
        % conditions
        for cc=1:2
            %% train data  (from GLM with trained HRF regressor on seen training data)
            % append features for  chromophores(hbo and hbr) and all channels without SS
            % select channels
            csel = chsel{sbj,2,tt};
            fvbuf = [];
            fvbuf = squeeze(FW{sbj,tt}(:,chrom,csel,TTM{sbj}.tnridx(tt,:),cc,rr));
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
            fvbuf = squeeze(FW{sbj,tt}(:,chrom,csel,TTM{sbj}.tstidx(tt),cc,rr));
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
flab = {'min ', 'max ', 'p2p ', 'avg ', 't2p ', 'slope', 'slope w2 '};
mlab = {'no GLM', 'GLM'};
fsel = {[], [1], [1], [2], [2], [3], [3], [4], [4], [6], [6],...
    [1,2], [1,2], [3,4], [3,4], [3,6], [3,6], [4,6], [4,6]};
%dimensionality of FW:
%Feature type | chromophore | channels | trials | condition | regressor
FW = {FMdc', FMss};
% for all features
for gg = 2:numel(fsel)
    % for all subjects
    for sbj=sbjl
        % for all trials
        for tt = 1:numel(TTM{sbj}.tstidx)
            if mod(gg,2)+1 == 1
                cvidx = 1;
            else
                cvidx = tt;
            end
            % channel selection
            csel = chsel{sbj,mod(gg,2)+1,tt};
            % initialize
            xTrF{gg,sbj,tt} =[];
            xTstF{gg,sbj,tt}=[];
            yTrF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tnridx(tt,:)));
            yTstF{gg,sbj,tt}=zeros(numel(epo.className),2*numel(TTM{sbj}.tstidx(tt)));
            for cc=1:2
                %% train data  (from GLM with trained HRF regressor on seen training data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf = FW{mod(gg,2)+1}{sbj,cvidx}(fsel{gg},chrom,csel,TTM{sbj}.tnridx(tt,:),cc);
                xTrF{gg,sbj,tt} = [xTrF{gg,sbj,tt} reshape(fvbuf, size(fvbuf,1)*size(fvbuf,2)*size(fvbuf,3),numel(TTM{sbj}.tnridx(tt,:)))];
                % generate label vector
                yTrF{gg,sbj,tt}(cc,(cc-1)*numel(TTM{sbj}.tnridx(tt,:))+1:cc*numel(TTM{sbj}.tnridx(tt,:)))=1;
                % test data (from GLM with trained HRF regressor on unseen data)
                % append features for hbo and hbr and all channels without SS
                fvbuf = [];
                fvbuf =FW{mod(gg,2)+1}{sbj,cvidx}(fsel{gg},chrom,csel,TTM{sbj}.tstidx(tt),cc,rr);
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
            %% R2 based channel selection?
            switch chselType 
                case 'r2'
                    nchanR2 = 25; % number of features to keep
                    fv.x = xTrF{gg,sbj,tt};
                    fv.y = yTrF{gg,sbj,tt};
                    fv.className = {'stim', 'rest'};   
                    fv = proc_rSquareSigned(fv);
                    [b,i] = sort(fv.x, 'descend');
                    if numel(b) > nchanR2
                       % keep only nchan features with highest r2
                       % training
                       xTrF{gg,sbj,tt} = xTrF{gg,sbj,tt}(i(1:nchanR2),:);
                       % testing
                       xTstF{gg,sbj,tt} = xTstF{gg,sbj,tt}(i(1:nchanR2),:);
                    end
            end
            
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

%% Perform paired tt-tests for classification results between methods
fno = floor(size(accuracyGLMF{hh},1)/2);
pvlab{1} ='';
for ff = 1:fno
    [h{hh}(ff),p{hh}(ff)]= ttest(accuracyGLMF{hh}(ff*2+1,sbjl), accuracyGLMF{hh}(ff*2,sbjl));
    if p{hh}(ff)>0.05
        pvlab{ff+1} ='';
    elseif p{hh}(ff) <=0.05
        pvlab{ff+1} ='*';
    elseif p{hh}(ff) <=1e-2
        pvlab{ff+1} ='**';
    elseif p{hh}(ff) <=1e-3
        pvlab{ff+1} ='***';
    end
end


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
    uflab{i} = string(join(flab(fidx{i})));
end
uflab{1} = '\beta GLM ';

figure
subplot(2,1,1)
hold on
bar(glmacc,'EdgeColor','none','BarWidth',1)
plot([.5 13.5], [50 50], '--k')
ylabel('sbj avg accuracy / %')
xlim([0.5 10.5])
xticks(1:13)
xticklabels(xtckglm)
% create upper x axis for feature type labels
ax1 = gca;
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlim([0.5, 10.5])
xticks(1:10)
yticks([])
xticklabels(uflab)
subplot(2,1,2)
hold on
xlim([0.5, 10.5])
xticks(1:10)
bar(noglmacc,'EdgeColor','none','BarWidth',1)
plot([.5 13.5], [50 50], '--k')
ylabel('avg accuracy / %')
xticks(1:10)
xticklabels(xtcknoglm)
% create upper x axis for pval labels
ax1 = gca;
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlim([0.5, 10.5])
xticks(1:10)
yticks([])
% significance and delta improvement
for i=2:numel(pvlab)
    pvlab{i} = [num2str(xtckglm(i)-xtcknoglm(i)) '  ' pvlab{i}]
end
xticklabels(pvlab)
set(gcf,'Position',[100 100 1000 250])

% print(gcf, '-dpdf', '-bestfit', 'classires.pdf')

disp(['Average accuracy improvement GLM vs no GLM: ' num2str(mean(xtckglm(2:end)-xtcknoglm(2:end))) '%'])