clear all;

malexflag = 0; % user flag
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

% use hrf STIM regressor weights as features (not the REST regressor
% weights, as they are useless here)
% and transform to bbci data structure
rr = 1;




for sbj=1:numel(TTM)
   
    for tt=1:numel(TTM{sbj}.tstidx
    
%% training of rLDA


%% testing of rLDA 
    


    end
    
end

