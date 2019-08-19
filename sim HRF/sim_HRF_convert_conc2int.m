%% get from concentrations back to intensities
figure
plot(hrf_conc)

%convert to OD
hrf_dod=hmrConc2OD(hrf_conc, hrf_SD, [6 6]);

figure
plot(hrf_dod)

%back to intensity (natural logarithm)
hrf_d= exp(-hrf_dod);

figure
plot(hrf_d)
