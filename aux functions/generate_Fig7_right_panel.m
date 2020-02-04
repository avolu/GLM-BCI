% First run: run_GLMcompare_noVsSSGLM.m for one subject.

  ix = 2
  epoch = 1
  % # of time points X HbO/R/T X # of channels X # of trials X conditions   % cc = 1 (stim) and cc = 2 (rest) X regressors rr=1 (hrf stim) and rr=2 (hrf rest)
  
  figure;
  % glm with ss
  plot(tHRF,squeeze(yavg_ss(:,1,lstHrfAdd{1}(ix,1),epoch,1,1)),'color',[1 0.5 0.4],'Linewidth',2);
  hold on;
  plot(tHRF,squeeze(yavg_ss(:,2,lstHrfAdd{1}(ix,1),epoch,1,1)),'color',[0.4 0.5 1],'Linewidth',2);
  
  % no glm
  plot(tHRF,(y_raw(:,1,lstHrfAdd{1}(ix,1),1,1)),'--','color',[1 0.5 0.4],'Linewidth',2);
  plot(tHRF,(y_raw(:,2,lstHrfAdd{1}(ix,1),1,1)),'--','color',[0.4 0.5 1],'Linewidth',2);
  
  % true hrf
  plot(hrf.t_hrf,hrf.hrf_conc(:,1:2),'k')
  xlabel('t / s')
  ylabel('\Delta C / Mol')
  xlim([-2 15]);
  legend('HbO - GLM with SS','HbR - GLM with SS','HbO - no GLM','HbR - no GLM')
  % get corr and mse
%   FMclab =
% 
%   1×9 cell array
% 
%     {'min'}    {'max'}    {'peak2peak'}    {'avg'}    {'time2peak'}    {'slope wdw1'}    {'slope wdw2'}    {'Corr'}    {'MSE'}

a = FMss{1,1};
%  9     3    56    15     2     2
corr_GLMss_HbO = a(8,1,lstHrfAdd{1}(ix,1),1,1,1) % Corr GLM wiht SS
corr_GLMss_HbR = a(8,2,lstHrfAdd{1}(ix,1),1,1,1) % Corr GLM wiht SS
MSE_GLMss_HbO = a(9,1,lstHrfAdd{1}(ix,1),1,1,1) % MSE GLM with SS
MSE_GLMss_HbR = a(9,2,lstHrfAdd{1}(ix,1),1,1,1) % MSE GLM with SS

b = FMdc{1,1};
%  9     3    56    15     2     2
corr_noGLM_HbO = b(8,1,lstHrfAdd{1}(ix,1),1,1,1) % Corr GLM wiht SS
corr_noGLM_HbR = b(8,2,lstHrfAdd{1}(ix,1),1,1,1) % Corr GLM wiht SS
MSE_noGLM_HbO = b(9,1,lstHrfAdd{1}(ix,1),1,1,1) % MSE GLM with SS
MSE_noGLM_HbR = b(9,2,lstHrfAdd{1}(ix,1),1,1,1) % MSE GLM with SS

%% display statistics
disp(['GLM with SS - HbO/HbR. MSE: ' num2str(MSE_GLMss_HbO) '/' num2str(MSE_GLMss_HbR) ' Corr: ' num2str(corr_GLMss_HbO) '/' num2str(corr_GLMss_HbR)])
disp(['no GLM - HbO/HbR. MSE: ' num2str(MSE_noGLM_HbO) '/' num2str(MSE_noGLM_HbR) ' Corr: ' num2str(corr_noGLM_HbO) '/' num2str(corr_noGLM_HbR)])

