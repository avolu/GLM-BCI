dat.x=yavg_ss;
dat.fs=25;
dat.t=-2:1/dat.fs:17;    
param.swdw=[0,4;10,17];

% test feature extraction

[FV] = featureExtract(dat, param);

% visualize/output data
figure
plot(dat.t,squeeze(dat.x(:,1,:)))

ch = 31;
figure
plot(dat.t,squeeze(dat.x(:,1,ch)))

format shortG

FV.x(:,1,ch)
FV.clab