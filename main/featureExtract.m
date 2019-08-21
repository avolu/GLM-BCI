function [FV] = featureExtract(dat, param)
%FEATUREEXTRACT extracts typical features from estimated HRF time series
%data
% INPUT
% dat.x:    data with dimensions  T x C x CH x E
%           T: # Time points, C: # chromophores, CH: # channels, E: #
%           epochs/trials
% dat.t     time, dimensions T x 1
% dat.fs:   sample frequency in Hz
% param:    parameters for feature extraction
%   .swdw:  [start, stop]   start and stop time (in ms) for window to calculate slopes
%           one row per window
%
% OUTPUT
% fv.x:     vector of features with dimensions F x C x CH x N
%           F: # of features, rest as input above
% fv.clab:  F x 1 cell array with labels for each feature


%% create feature labels
FV.clab = {'min','max','peak2peak','avg','time2peak'};
for i=1:size(param.swdw,1)
    FV.clab{end+1} = ['slope wdw' num2str(i)];
end

%% find window indices from times
for ss=1:size(param.swdw,1)
    [minValue,closestIndex] = min(abs(dat.t-param.swdw(ss,1)));
    wdw(ss,1)=closestIndex;
    [minValue,closestIndex] = min(abs(dat.t-param.swdw(ss,2)));
    wdw(ss,2)=closestIndex;
end

%% All concentration features in muMol for better precision
dat.x = dat.x*1e6;
nchan = size(dat.x,3);

for e=1:size(dat.x, 4)
    %% min
    FV.x(1,:,1:nchan,e) = min(dat.x(:,:,:,e));
    %% max
    FV.x(2,:,1:nchan,e) = max(dat.x(:,:,:,e));
    %% peak2peak
    FV.x(3,:,1:nchan,e) = max(dat.x(:,:,:,e))-min(dat.x(:,:,:,e));
    %% avg
    FV.x(4,:,1:nchan,e) = mean(dat.x(:,:,:,e));
    for ch=1:size(dat.x,3)
        for c=1:size(dat.x,2)
            %% time2peak
            m = NaN;
            i = NaN;
            switch c
                % HbO
                case 1
                    [m,i] = max(dat.x(:,c,ch,e));
                    % HbR
                case 2
                    [m,i] = min(dat.x(:,c,ch,e));
                    % HbTot
                case 3
                    [m,i] = max(dat.x(:,c,ch,e));
            end
            FV.x(5,c,ch,e) = dat.t(i);
            
            % verifying slope
            %             if ch == 31
            %                figure
            %                plot(dat.t,squeeze(dat.x(:,c,ch,e)))
            %             end
            %% slope
            for ss=1:size(wdw,1)
                p=polyfit(dat.t(wdw(ss,1):wdw(ss,2))', squeeze(dat.x(wdw(ss,1):wdw(ss,2),c,ch,e)),1);
                FV.x(5+ss,c,ch,e)=p(1);
                % verifying slope
                %                 if ch == 31
                %                     hold on
                %                 plot(dat.t(wdw(ss,1):wdw(ss,2)), polyval(p, dat.t(wdw(ss,1):wdw(ss,2))))
                %                 end
            end
        end
    end
end

