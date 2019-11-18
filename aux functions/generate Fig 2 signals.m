cd('\\ad\eng\users\m\a\mayucel\Desktop\BCI_paper_figure');
% subj 37!
load block_design.nirs -mat
load clmn_bci_paper.mat
fq = 1/(t(2)-t(1));
[dc, dod] = hmrIntensity2Conc( d, SD, fq, [], [], [6  6]);
pre = 135*fq;
post = 225*fq;
ch = 2;
figure; 
subplot(7,1,1); plot(t(1:size([pre:post],2)),s(pre:post,:)); hold on;
plot(t(1:size([pre:post],2)),clmn(pre:post,:),'r')
set(gca,'xtick',[],'ytick',[]);title('Representative brain activity');

subplot(7,1,2); plot(t(1:size([pre:post],2)),aux(pre:post,2:4)./max(aux(pre:post,2:4)),'g')
set(gca,'xtick',[],'ytick',[]);title('Accelerometer'); ylim([0.95 1]);

subplot(7,1,3); plot(t(1:size([pre:post],2)),aux(pre:post,5),'g')
set(gca,'xtick',[],'ytick',[]);title('PPG');

subplot(7,1,4); plot(t(1:size([pre:post],2)),aux(pre:post,6),'g')
set(gca,'xtick',[],'ytick',[]);title('BP');

subplot(7,1,5); plot(t(1:size([pre:post],2)),aux(pre:post,7),'g')
set(gca,'xtick',[],'ytick',[]);title('RESP');

subplot(7,1,6); plot(t(1:size([pre:post],2)),dc(pre:post,1,17),'g')
set(gca,'xtick',[],'ytick',[]);title('Short separation NIRS signal');

subplot(7,1,7); plot(t(1:size([pre:post],2)),dc(pre:post,1,ch))
set(gca,'ytick',[]);title('Long separation NIRS signal'); xlabel('time (sec)');