% user: 1 Meryem | 0 Alex
melexflag = 1;
if melexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\GLM_BCI_PAPER'; % data directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\GLM-BCI'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % data directory
end

path.fname = 'resting.nirs';
path.hrfname = 'hrf_simdat_50_shorterHRF.mat';
path.savename = 'resting_sim_50_shorterHRF.nirs';
addpath(genpath(path.code));
sbjfolder = {'Subj86','Subj91','Subj92','Subj94','Subj95','Subj96','Subj97','Subj98','Subj99','Subj100','Subj101','Subj102','Subj103','Subj104'};


for ss = 1:numel(sbjfolder) % loop across subjects
% load participant data
nirs = load([path.dir filesep 'RESTING_DATA' filesep sbjfolder{ss} filesep path.fname], '-mat');
% load simulated hrf
hrf = load([path.code filesep 'sim HRF' filesep path.hrfname]);
% hrf = load([path.dir path.hrfname]);

%% call the function
flag_prune = true;
[nirs_hrf] = addSimHRF(nirs, hrf, false, flag_prune);

save([path.dir filesep 'RESTING_DATA' filesep sbjfolder{ss} filesep path.savename], '-struct', 'nirs_hrf')
disp('data saved.')
end





