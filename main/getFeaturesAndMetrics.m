function [FM, clab] = getFeaturesAndMetrics(X, param, ival, hrf)
%Calculates all features from single trial fNIRS data as well as
%correlation/mse metrics using ground truth hrf
% INPUT
% X         input data with dimensions  T x C x CH x E
%           T: # Time points, C: # chromophores, CH: # channels, E: #
%           epochs/trials
% param     parameters for feature extraction
%   .swdw:  [start, stop]   start and stop time (in ms) for window to calculate slopes
%           one row per window
% ival      Interval/window for trial in secods [start end] 
%hrf        Ground truth hrf concentration data struct, assumes same sample rate as X 
%           with dimensions T x C -- T: # Time points, C: # chromophores
% OUTPUT
% FM:       Matrix of F features+ M metrics with dimensions F+M x C x CH x E
% clab:     F+M x 1 cell array with labels for each feature/metric


%% Generate inputs for featureExtract
% dat       input data
%    .x:    data with dimensions  T x C x CH x E
%           T: # Time points, C: # chromophores, CH: # channels, E: #
%           epochs/trials
%   .t:     time, dimensions T x 1
%   .fs:    sample frequency in Hz
% param     parameters for feature extraction
%   .swdw:  [start, stop]   start and stop time (in ms) for window to calculate slopes
%           one row per window
dat.x=X;
dat.fs=25;
dat.t=ival(1):1/dat.fs:ival(2);   
%% get features 
[FV] = featureExtract(dat, param);

%% find common timebase for ground truth hrf and data, assuming gt hrf starts at t=0
t0 = find(dat.t==0);
tlen = numel(dat.t(t0:end));

%% calculate ground truth metrics
FM = FV.x;
for e=1:size(dat.x, 4)
    for ch=1:size(dat.x,3)
        for c=1:size(dat.x,2)
            % correlation
            FM(numel(FV.clab)+1,c,ch,e) = corr(squeeze(dat.x(t0:end,c,ch,e)),squeeze(hrf.hrf_conc(1:tlen,c)));
            % mse in muMol
            FM(numel(FV.clab)+2,c,ch,e) = 1e6*sqrt(nanmean((squeeze(dat.x(t0:end,c,ch,e))-squeeze(hrf.hrf_conc(1:tlen,c))).^2));
        end
    end    
end
%expand clab
clab = FV.clab ;
clab{numel(FV.clab)+1}='Corr';
clab{numel(FV.clab)+2}='MSE';

end

