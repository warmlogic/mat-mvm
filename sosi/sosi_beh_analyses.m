figFontName = 'Arial';

subjects = {
  'SOSI001';
  'SOSI002';
  'SOSI003';
  'SOSI004';
  'SOSI005';
  'SOSI006';
  'SOSI007';
  'SOSI008';
  'SOSI009';
  'SOSI010';
  'SOSI011';
  'SOSI012';
  'SOSI013';
  'SOSI014';
  'SOSI015';
  'SOSI016';
  'SOSI017';
  'SOSI018';
  'SOSI020';
  'SOSI019';
  'SOSI021';
  'SOSI022';
  'SOSI023';
  'SOSI024';
  'SOSI025';
  'SOSI026';
  'SOSI027';
  'SOSI028';
  'SOSI029';
  'SOSI030';
  };
% original SOSI019 was replaced because the first didn't finish

%badSub = {};
badSub = {'SOSI011','SOSI030','SOSI007','SOSI001'};

% 005 might be bad

% all events
SOSI_RS_WIR = [0.9441 0.8333 0.9444 0.9302 0.7884 0.8621 0.5323 0.9349 0.9891 0.9878 0.9121 0.8398 0.9121 0.9789 0.8958 0.809 0.764 0.971 0.9298 0.817 0.8846 0.7308 0.9412 0.641 0.988 0.9679 0.9942 0.9767 0.8951 0.8339];
SOSI_RO_WIR = [0.7333 0.8636 0.9987 0.6667 0.7143 0.5354 0.4393 0.45 0.7391 0.8088 0.575 0.7273 0.6739 0.6792 0.625 0.6441 0.8681 0.6087 0.625 0.6053 0.6304 0.4898 0.6533 0.726 0.6825 0.5942 0.899 0.8395 0.625 0.5417];
SOSI_F_WIR = [0.7727 0.5741 0.5778 0.6094 0.6389 0.58 0.4286 0.4444 0.6842 0.5 0.0013 0.5135 0.5303 0.541 0.5 0.5067 0.6279 0.7037 0.5118 0.4667 0.5795 0.5882 0.4074 0.7368 0.5238 0.4468 0.5977 0.6 0.625 0.0013];

% plot
SOSI_RS_avg = mean(SOSI_RS_WIR(~ismember(subjects,badSub)));
SOSI_RO_avg = mean(SOSI_RO_WIR(~ismember(subjects,badSub)));
SOSI_F_avg = mean(SOSI_F_WIR(~ismember(subjects,badSub)));

SOSI_RS_sem = std(SOSI_RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
SOSI_RO_sem = std(SOSI_RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
SOSI_F_sem = std(SOSI_F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));

[h,p,ci,stats] = ttest(SOSI_RS_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('SOSI_RS_WIR: t(%d)=%.4f, p=%.10f\n',stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(SOSI_RO_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('SOSI_RO_WIR: t(%d)=%.4f, p=%.10f\n',stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(SOSI_F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('SOSI_F_WIR: t(%d)=%.4f, p=%.10f\n',stats.df,stats.tstat,p);


bw_groupnames = {'Rem. Source';'Rem. Other';'Familiar'};
bw_title = 'Proportion of Source Correct responses';
bw_legend = {'Location'};
%bw_colormap = 'gray';
bw_colormap = [.5 .5 .5];
%bw_data = [SOSI_RS_avg 0; SOSI_RO_avg 0; SOSI_F_avg 0];
%bw_errors = [SOSI_RS_sem 0; SOSI_RO_sem 0; SOSI_F_sem 0];
bw_data = [SOSI_RS_avg; SOSI_RO_avg; SOSI_F_avg];
bw_errors = [SOSI_RS_sem; SOSI_RO_sem; SOSI_F_sem];
bw_xlabel = 'RK Response';
bw_ylabel = 'Proportion Correct';

cfg_plot = [];
% set up how the lines will look
cfg_plot.linewidth = 2;
cfg_plot.errwidth = 2;
cfg_plot.errspec = 'k.';

figure
% plot the stuff
bar(bw_data,'LineWidth',cfg_plot.linewidth);
colormap(bw_colormap);
hold on
errorbar(bw_data,bw_errors,cfg_plot.errspec,'LineWidth',cfg_plot.errwidth);
% set up the information
set(gca, 'XTickLabel', bw_groupnames, 'box', 'off', 'ticklength', [0 0]);
legend(bw_legend,'Location','NorthEast');
legend boxoff
axis([0.5 3.5 0 1]);
title(bw_title);
xlabel(bw_xlabel)
ylabel(bw_ylabel)
publishfig(gcf,0,[],[],figFontName);
% horiz chance line
plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2);
hold off
% print it
print(gcf,'-dpng','~/Desktop/SOSI_RS_RO_F_accuracy');
print(gcf,'-depsc2','~/Desktop/SOSI_RS_RO_F_accuracy');

% figure
% bw_groupnames = {'Rem. Source';'Rem. Other';'Familiar'};
% bw_title = 'Proportion of Source Correct responses';
% bw_legend = {'Side'};
% %bw_colormap = 'gray';
% bw_colormap = [.5 .5 .5];
% bw_data = [SOSI_RS_avg 0; SOSI_RO_avg 0; SOSI_F_avg 0];
% bw_errors = [SOSI_RS_sem 0; SOSI_RO_sem 0; SOSI_F_sem 0];
% h = barweb(bw_data,bw_errors,[],bw_groupnames,bw_title,[],[],bw_colormap,[],bw_legend);
% set(h.legend,'Location','NorthEast');
% axis([0.5 3.5 0 1]);
% publishfig(gcf,0,[],[],figFontName);
% hold on
% plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2); % horiz chance line
% %print(gcf,'-dpng','~/Desktop/SOSI_RS_RO_F_accuracy');
