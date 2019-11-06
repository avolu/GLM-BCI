clear all
load('D:\Office\Research\Software - Scripts\Matlab\GLM-BCI\lit rev data\piechart_data.mat')

%define indices of papers to be used for charts
pidx = find(WebofScienceTop100Search == 'Yes');
% take first 100 results
pidx = pidx(1:100);
% take out reviews
pidx(find(contains(Notes(pidx),'review', 'IgnoreCase', true))) = [];

figure
%% GLM pie chart
cats = categories(GLM(pidx));
cats{end+1} = 'N/A';
undef = isundefined(GLM(pidx));
for ii = 1:numel(cats)-1
   glm(ii)= sum((GLM(pidx)==cats(ii)));
end
glm(numel(cats)) = numel(find(undef));
%flip kalman and no
b=glm(2);
glm(2) = glm(1);
glm(1)=b;
b = cats{2};
cats{2}=cats{1};
cats{1}=b;
H = pie(glm)
title('GLM-based Analysis')
T = H(strcmpi(get(H,'Type'),'text'));
P = cell2mat(get(T,'Position'));
set(T,{'Position'},num2cell(P*0.6,2))
text(P(:,1),P(:,2),cats(:))

%% Noise removal pie chart
figure
cats = categories(SystemicPhysNoiseRemoval(pidx));
cats{end+1} = 'N/A';
undef = isundefined(SystemicPhysNoiseRemoval(pidx));
for ii = 1:numel(cats)-1
   sysnoise(ii)= sum((SystemicPhysNoiseRemoval(pidx)==cats(ii)));
end
sysnoise(numel(cats)) = numel(find(undef));
H = pie(sysnoise)
title('Physiological Noise Removal beyond HP/LP Filtering')
T = H(strcmpi(get(H,'Type'),'text'));
P = cell2mat(get(T,'Position'));
set(T,{'Position'},num2cell(P*0.6,2))
text(P(:,1),P(:,2),cats(:))


figure
%% HbO/HbR pie chart
subplot(1,2,1)
cats = categories(ChromophoresusedforFeatures(pidx));
cats{end+1} = 'N/A';
undef = isundefined(ChromophoresusedforFeatures(pidx));
for ii = 1:numel(cats)-1
   chrom(ii)= sum((ChromophoresusedforFeatures(pidx)==cats(ii)));
end
chrom(numel(cats)) = numel(find(undef));
H = pie(chrom)
title('Chromophores')
T = H(strcmpi(get(H,'Type'),'text'));
P = cell2mat(get(T,'Position'));
set(T,{'Position'},num2cell(P*0.6,2))
text(P(:,1),P(:,2),cats(:))


%% Features pie chart
subplot(1,2,2)
fnames = {'Peak', 'Average', 'Slope', 'Peak2Peak', 'Time2Peak', 'Connectivity', 'Other'};
fDatVars = {FeaturePeak(pidx), FeatureAverage(pidx), FeatureSlope(pidx), FeatureRangeP2P(pidx), ...
    FeaturetimetoPeak(pidx), FeatureConnectivity(pidx), FeatureOther(pidx)};
for dd = 1:numel(fDatVars)-1
    features(dd) = numel(find(fDatVars{dd}=='Yes'))/numel(pidx)*100;
end
features(dd+1) = (numel(pidx)-numel(find(fDatVars{end}=='')))/numel(pidx)*100;
bar(features)
ylabel('feature used in % of publications')
title('Features')
xticklabels(fnames)
xtickangle(25)


%% Classifier pie chart
figure
cats = categories(Classifierused(pidx));
cats{end+1} = 'N/A';
undef = isundefined(Classifierused(pidx));
for ii = 1:numel(cats)-1
   class(ii)= sum((Classifierused(pidx)==cats(ii)));
end
class(numel(cats)) = numel(find(undef));
H = pie(class)
title('Classifiers')
T = H(strcmpi(get(H,'Type'),'text'));
P = cell2mat(get(T,'Position'));
set(T,{'Position'},num2cell(P*0.6,2))
text(P(:,1),P(:,2),cats(:))



%% SS Channel removal pie chart
figure
cats = categories(ShortSeparationChannelsused(pidx));
cats{end+1} = 'N/A';
undef = isundefined(ShortSeparationChannelsused(pidx));
for ii = 1:numel(cats)-1
   ssep(ii)= sum((ShortSeparationChannelsused(pidx)==cats(ii)));
end
ssep(numel(cats)) = numel(find(undef));
H = pie(ssep)
title('Short Separation Channels Used')
T = H(strcmpi(get(H,'Type'),'text'));
P = cell2mat(get(T,'Position'));
set(T,{'Position'},num2cell(P*0.6,2))
text(P(:,1),P(:,2),cats(:))


