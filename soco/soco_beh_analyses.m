
subjects = {
  'SOCO001';
  'SOCO002';
  'SOCO003';
  'SOCO004';
  'SOCO005';
  'SOCO006';
  'SOCO007';
  'SOCO008';
  'SOCO009';
  'SOCO010';
  'SOCO011';
  'SOCO012';
  'SOCO013';
  'SOCO014';
  'SOCO015';
  'SOCO016';
  'SOCO017';
  'SOCO018';
  'SOCO019';
  'SOCO020';
  'SOCO021';
  'SOCO022';
  'SOCO023';
  'SOCO024';
  'SOCO025';
  'SOCO026';
  'SOCO027';
  'SOCO028';
  'SOCO029';
  'SOCO030';
  };
% SOCO002 ended early by 6(?) trials because of fire alarm

badSub = {'SOCO018','SOCO026'};

% allevents
C2_RS_WIR = [0.6304 0.7093 0.9323 0.9412 0.8333 0.4957 0.8769 0.8299 0.9322 0.7662 0.7875 0.7143 0.7674 0.5345 0.9167 0.8661 0.9975 0.6 0.9263 0.6263 0.6739 0.5455 0.9737 0.8182 0.7164 0.9975 0.7812 0.7333 0.5484 0.7857];
C2_RO_WIR = [0.3269 0.566 0.6176 0.9 0.4386 0.6923 0.5269 0.5455 0.4466 0.5455 0.5 0.5312 0.475 0.45 0.4615 0.5694 0.5849 0.5192 0.4615 0.5696 0.5882 0.5412 0.5957 0.5385 0.5122 0.6 0.6019 0.625 0.4565 0.4714];
C2_F_WIR = [0.5 0.7 0.5 0.5 0.7273 0.9975 0.5135 0.7143 0.5 0.0025 0.7 0.5806 0.561 0.6471 0.5 0.6667 0.6353 0.5 0.6 0.5 0.5217 0.6 0.5312 0.3571 0.5417 0.4091 0.25 0.5773 0.4894 0.5217];

C6_RS_WIR = [0.6119 0.8 0.9048 0.8333 0.6714 0.5091 0.9059 0.8168 0.9143 0.7415 0.9519 0.7576 0.7286 0.8904 0.8148 0.9391 0.963 0.5882 0.9024 0.6306 0.6496 0.3793 0.9556 0.88 0.8056 0.6 0.76 0.7742 0.5974 0.6364];
C6_RO_WIR = [0.4795 0.4082 0.5345 0.7561 0.5 0.4773 0.6479 0.6 0.598 0.5238 0.62 0.5283 0.5263 0.6216 0.6143 0.5745 0.5417 0.551 0.4528 0.5588 0.6207 0.5263 0.6386 0.5738 0.4643 0.5676 0.4321 0.6 0.4872 0.5];
C6_F_WIR = [0.375 0.5833 0.3571 0.5769 0.6364 0.4286 0.6111 0.75 0.4615 0.3333 0.6562 0.5147 0.4091 0.5 0.3684 0.6957 0.5068 0.625 0.2105 0.0025 0.4545 0.4545 0.4828 0.4583 0.8 0.551 0.7143 0.5242 0.5581 0.525];

% % RT < 4000 ms
% C2_RS_WIR = [0.6264 0.7059 0.9323 0.9412 0.8333 0.4957 0.8889 0.8299 0.9322 0.7662 0.7875 0.7143 0.7674 0.5345 0.9167 0.8661 0.9975 0.6 0.9263 0.6263 0.6739 0.5455 0.9737 0.8182 0.7231 0.9975 0.7812 0.7333 0.5484 0.7857];
% C2_RO_WIR = [0.3269 0.5686 0.625 0.8992 0.463 0.6923 0.5325 0.5455 0.451 0.5455 0.5 0.5397 0.4474 0.45 0.4574 0.5694 0.5849 0.52 0.4737 0.5696 0.5882 0.5412 0.5957 0.5385 0.5185 0.6 0.6019 0.625 0.4565 0.4714];
% C2_F_WIR = [0.5556 0.7 0.5 0.5 0.7273 0.9975 0.4839 0.7143 0.5 0.0025 0.7778 0.5806 0.561 0.6471 0.5217 0.6667 0.6353 0.5 0.6 0.5 0.5217 0.6 0.5312 0.3571 0.5417 0.4091 0.25 0.5729 0.5 0.5217];
% 
% C6_RS_WIR = [0.6119 0.8125 0.9109 0.825 0.6714 0.5091 0.9079 0.8168 0.913 0.7415 0.951 0.7812 0.7286 0.8889 0.8148 0.9391 0.963 0.5714 0.9012 0.6355 0.6496 0.3929 0.9556 0.88 0.8056 0.6 0.7755 0.7742 0.5974 0.6562];
% C6_RO_WIR = [0.4795 0.3958 0.5098 0.7561 0.5 0.4884 0.6136 0.6 0.598 0.5238 0.6087 0.5385 0.5263 0.6143 0.6143 0.5778 0.5417 0.5333 0.4528 0.5588 0.6207 0.5357 0.6386 0.5702 0.4643 0.5676 0.4304 0.6 0.4872 0.4921];
% C6_F_WIR = [0.375 0.5833 0.3846 0.56 0.6364 0.3333 0.6429 0.75 0.4615 0.3333 0.6333 0.5075 0.4091 0.5294 0.3684 0.6957 0.5068 0.625 0.2222 0.0026 0.4545 0.4545 0.4828 0.4583 0.8 0.551 0.6923 0.5242 0.5581 0.5278];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separate 2/6 colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C2_RS_avg = mean(C2_RS_WIR(~ismember(subjects,badSub)));
C2_RO_avg = mean(C2_RO_WIR(~ismember(subjects,badSub)));
C2_F_avg = mean(C2_F_WIR(~ismember(subjects,badSub)));

C2_RS_sem = std(C2_RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
C2_RO_sem = std(C2_RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
C2_F_sem = std(C2_F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));

C6_RS_avg = mean(C6_RS_WIR(~ismember(subjects,badSub)));
C6_RO_avg = mean(C6_RO_WIR(~ismember(subjects,badSub)));
C6_F_avg = mean(C6_F_WIR(~ismember(subjects,badSub)));

C6_RS_sem = std(C6_RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
C6_RO_sem = std(C6_RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
C6_F_sem = std(C6_F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));

[h,p,ci,stats] = ttest(C2_RS_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('C2_RS_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(C2_RS_avg),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(C2_RO_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('C2_RO_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(C2_RO_avg),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(C2_F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('C2_F_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(C2_F_avg),stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(C6_RS_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('C6_RS_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(C6_RS_avg),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(C6_RO_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('C6_RO_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(C6_RO_avg),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(C6_F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('C6_F_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(C6_F_avg),stats.df,stats.tstat,p);

figure
bw_groupnames = {'Rem. Src.';'Rem. Other';'Familiar'};
bw_title = 'Proportion of Source Correct responses';
bw_legend = {'2 Colors','6 Colors'};
bw_colormap = 'gray';
bw_data = [C2_RS_avg, C6_RS_avg; C2_RO_avg, C6_RO_avg; C2_F_avg, C6_F_avg];
bw_errors = [C2_RS_sem, C6_RS_sem; C2_RO_sem, C6_RO_sem; C2_F_sem, C6_F_sem];
h = barweb(bw_data,bw_errors,[],bw_groupnames,bw_title,[],[],bw_colormap,[],bw_legend);
set(h.legend,'Location','NorthEast');
axis([0.5 3.5 0 1]);
publishfig(gcf,0);
hold on
plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2); % horiz chance line
print(gcf,'-dpng','~/Desktop/SOCO_C6_C2_RS_RO_F_accuracy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collapse across colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RS_WIR = (C2_RS_WIR + C6_RS_WIR)./2;
RO_WIR = (C2_RO_WIR + C6_RO_WIR)./2;
F_WIR = (C2_F_WIR + C6_F_WIR)./2;

RS_avg = mean(RS_WIR(~ismember(subjects,badSub)));
RO_avg = mean(RO_WIR(~ismember(subjects,badSub)));
F_avg = mean(F_WIR(~ismember(subjects,badSub)));

RS_sem = std(RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
RO_sem = std(RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
F_sem = std(F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));

[h,p,ci,stats] = ttest(RS_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('RS_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(RS_avg),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(RO_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('RO_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(RO_avg),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(F_WIR(~ismember(subjects,badSub)),0.5*ones(1,sum(~ismember(subjects,badSub))),0.05,'both');
fprintf('F_WIR (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(F_avg),stats.df,stats.tstat,p);

bw_groupnames = {'Rem. Src.';'Rem. Other';'Familiar'};
bw_title = 'Proportion of Source Correct responses';
bw_legend = {'Color'};
%bw_colormap = 'gray';
bw_colormap = [.5 .5 .5];
%bw_data = [RS_avg 0; RO_avg 0; F_avg 0];
%bw_errors = [RS_sem 0; RO_sem 0; F_sem 0];
bw_data = [RS_avg; RO_avg; F_avg];
bw_errors = [RS_sem; RO_sem; F_sem];

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
xlabel('RK Response')
ylabel('Proportion Correct')
publishfig(gcf,0);
% horiz chance line
plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2);
hold off
% print it
print(gcf,'-dpng','~/Desktop/SOCO_RS_RO_F_accuracy');

% figure
% bw_groupnames = {'Rem. Src.';'Rem. Other';'Familiar'};
% bw_title = 'Proportion of Source Correct responses';
% bw_legend = {'Color'};
% %bw_colormap = 'gray';
% bw_colormap = [.5 .5 .5];
% bw_data = [RS_avg 0; RO_avg 0; F_avg 0];
% bw_errors = [RS_sem 0; RO_sem 0; F_sem 0];
% h = barweb(bw_data,bw_errors,[],bw_groupnames,bw_title,[],[],bw_colormap,[],bw_legend);
% set(h.legend,'Location','NorthEast');
% axis([0.5 3.5 0 1]);
% publishfig(gcf,0);
% hold on
% plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2); % horiz chance line
% print(gcf,'-dpng','~/Desktop/SOCO_RS_RO_F_accuracy');
