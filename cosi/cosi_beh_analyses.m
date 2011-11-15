% COSI behavioral analyses

figFontName = 'Helvetica';

subjects = {
  'COSI001';
  'COSI002';
  'COSI003';
  'COSI004';
  'COSI005';
  'COSI006';
  'COSI007';
%   'COSI008';
%   'COSI009';
%   'COSI010';
  'COSI011';
  'COSI012';
  'COSI013';
  'COSI014';
  'COSI015';
  'COSI016';
  'COSI017';
  'COSI018';
  'COSI019';
  'COSI020';
%   'COSI021';
  'COSI022';
  'COSI023';
  'COSI024';
  'COSI025';
  'COSI026';
  'COSI027';
  'COSI028';
  'COSI029';
  'COSI030';
  'COSI031';
  'COSI032';
  'COSI033';
  'COSI034';
  'COSI035';
  };

badSub = {'COSI031','COSI002','COSI035','COSI007','COSI030'};


COSI_COLOR_RS_WIR = [0.8705 0.75 0.7296 0.7053 0.8665 0.6452 0.9738 0.952 0.7225 0.6411 0.977 0.8218 0.7244 0.6779 0.873 0.6091 0.881 0.62 0.8826 0.7529 0.8069 0.8393 0.9637 0.8875 0.6709 0.938 0.5655 0.7082 0.8312 0.7474 0.908];
COSI_COLOR_RO_WIR = [0.5705 0.8225 0.5631 0.5779 0.6018 0.505 0.6839 0.632 0.8321 0.553 0.6913 0.5306 0.4797 0.4797 0.5544 0.5309 0.5575 0.4927 0.6697 0.5349 0.5824 0.4905 0.5095 0.6191 0.4641 0.5648 0.4947 0.5567 0.6469 0.6608 0.825];
COSI_COLOR_F_WIR = [0.5476 0.5688 0.3333 0.4688 0.651 0.4736 0.616 0.5599 0.6006 0.5705 0.5361 0.5179 0.4737 0.4472 0.625 0.2512 0.5979 0.4196 0.3667 0.5226 0.5448 0.6148 0.5329 0.4868 0.5613 0.5522 0.4958 0.4583 0.5472 0.625 0.3194];

COSI_SIDE_RS_WIR = [0.9624 0.8214 0.9172 0.9109 0.9924 0.781 0.9975 0.9928 0.9782 0.9354 0.9975 0.8738 0.8755 0.8432 0.987 0.6525 0.9405 0.7568 0.909 0.9279 0.8551 0.873 0.9975 0.9558 0.7449 0.9975 0.6789 0.8728 0.9975 0.9773 0.89];
COSI_SIDE_RO_WIR = [0.5742 0.8732 0.7111 0.6234 0.8816 0.5603 0.6981 0.8333 0.5 0.5869 0.8219 0.6947 0.5278 0.515 0.6099 0.4529 0.6137 0.4457 0.5731 0.5681 0.7438 0.6281 0.7795 0.5625 0.6548 0.8818 0.5971 0.6199 0.7925 0.7787 0.8545];
COSI_SIDE_F_WIR = [0.4062 0.5987 0.5697 0.5976 0.5684 0.4643 0.6828 0.6378 0.5774 0.5033 0.4858 0.4821 0.5858 0.5856 0.5936 0.4 0.5702 0.5631 0.3618 0.6872 0.5851 0.5208 0.6377 0.4913 0.517 0.6549 0.5488 0.4118 0.6311 0.5165 0.6583];

COSI_COLOR_RS_avg = mean(COSI_COLOR_RS_WIR(~ismember(subjects,badSub)));
COSI_COLOR_RO_avg = mean(COSI_COLOR_RO_WIR(~ismember(subjects,badSub)));
COSI_COLOR_F_avg = mean(COSI_COLOR_F_WIR(~ismember(subjects,badSub)));

COSI_COLOR_RS_sem = std(COSI_COLOR_RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
COSI_COLOR_RO_sem = std(COSI_COLOR_RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
COSI_COLOR_F_sem = std(COSI_COLOR_F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));

COSI_SIDE_RS_avg = mean(COSI_SIDE_RS_WIR(~ismember(subjects,badSub)));
COSI_SIDE_RO_avg = mean(COSI_SIDE_RO_WIR(~ismember(subjects,badSub)));
COSI_SIDE_F_avg = mean(COSI_SIDE_F_WIR(~ismember(subjects,badSub)));

COSI_SIDE_RS_sem = std(COSI_SIDE_RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
COSI_SIDE_RO_sem = std(COSI_SIDE_RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
COSI_SIDE_F_sem = std(COSI_SIDE_F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));

% compare color to chance
[h,p,ci,stats] = ttest(COSI_COLOR_RS_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_COLOR_RS_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_COLOR_RS_WIR),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI_COLOR_RO_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_COLOR_RO_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_COLOR_RO_WIR),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI_COLOR_F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_COLOR_F_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_COLOR_F_WIR),stats.df,stats.tstat,p);

% compare side to chance
[h,p,ci,stats] = ttest(COSI_SIDE_RS_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_SIDE_RS_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_SIDE_RS_WIR),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI_SIDE_RO_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_SIDE_RO_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_SIDE_RO_WIR),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI_SIDE_F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_SIDE_F_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_SIDE_F_WIR),stats.df,stats.tstat,p);

% compare side and color accuracy
[h,p,ci,stats] = ttest(COSI_COLOR_F_WIR(~ismember(subjects,badSub)),COSI_SIDE_F_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('COSI: COLOR (M=%.2f) vs SIDE (M=%.2f) F_WIR: t(%d)=%.4f, p=%.10f\n',mean(COSI_COLOR_F_WIR),mean(COSI_SIDE_F_WIR),stats.df,stats.tstat,p);

%% plot
figure
bw_groupnames = {'Rem. Src.';'Rem. Other';'Familiar'};
bw_title = 'Proportion of Source Correct responses';
bw_legend = {'Color','Side'};
bw_colormap = 'gray';
bw_data = [COSI_COLOR_RS_avg, COSI_SIDE_RS_avg; COSI_COLOR_RO_avg, COSI_SIDE_RO_avg; COSI_COLOR_F_avg, COSI_SIDE_F_avg];
bw_errors = [COSI_COLOR_RS_sem, COSI_SIDE_RS_sem; COSI_COLOR_RO_sem, COSI_SIDE_RO_sem; COSI_COLOR_F_sem, COSI_SIDE_F_sem];
h = barweb(bw_data,bw_errors,[],bw_groupnames,bw_title,[],[],bw_colormap,[],bw_legend);
set(h.legend,'Location','NorthEast');
axis([0.5 3.5 0 1]);
publishfig(gcf,0,[],[],figFontName);
hold on
plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2); % horiz chance line
%print(gcf,'-dpng','~/Desktop/COSI_SIDE_COLOR_RS_RO_F_accuracy');


%% collapse across RO+F

COSI_COLOR_RO_F_WIR = mean([COSI_COLOR_RO_WIR;COSI_COLOR_F_WIR],1);
COSI_SIDE_RO_F_WIR = mean([COSI_SIDE_RO_WIR;COSI_SIDE_F_WIR],1);

% compare color to chance
[h,p,ci,stats] = ttest(COSI_COLOR_RO_F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_COLOR_RO_F_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_COLOR_RO_F_WIR),stats.df,stats.tstat,p);

% compare side to chance
[h,p,ci,stats] = ttest(COSI_SIDE_RO_F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('COSI_SIDE_RO_F_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(COSI_SIDE_RO_F_WIR),stats.df,stats.tstat,p);

% compare side and color accuracy
[h,p,ci,stats] = ttest(COSI_COLOR_RO_F_WIR(~ismember(subjects,badSub)),COSI_SIDE_RO_F_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('COSI: COLOR (M=%.2f) vs SIDE (M=%.2f) RO+F WIR: t(%d)=%.4f, p=%.10f\n',mean(COSI_COLOR_RO_F_WIR),mean(COSI_SIDE_RO_F_WIR),stats.df,stats.tstat,p);
