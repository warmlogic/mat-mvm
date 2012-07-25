figFontName = 'Helvetica';

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
  'SOSI019';
  'SOSI020';
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

% d'
SOSI_ITEM_DP = [2.4055, 1.0515, 1.9568, 1.3392, 0.875, 1.2222, 0.185, 1.9591, 1.9172, 1.1959, 2.4869, 0.9108, 1.0227, 1.705, 1.5148, 0.9685, 0.8106, 2.0084, 1.1812, 1.3761, 1.1057, 1.1521, 2.297, 0.4484, 1.4344, 1.5941, 1.398, 1.7473, 1.5984, 1.2503];
SOSI_SOURCE_DP = [2.6106, 1.2027, 1.515, 1.6033, 1.3849, 1.0842, -0.0091, 1.8163, 2.6303, 1.1873, 2.3238, 1.3773, 1.6582, 1.5506, 1.6074, 0.8615, 1.5816, 2.6045, 1.1413, 1.1, 1.4382, 0.8814, 2.1249, 1.0363, 1.6679, 1.5655, 2.2779, 2.1689, 1.6759, 1.6488];

% accuracy: all events
SOSI_RS_WIR = [0.9441 0.8333 0.9444 0.9302 0.7884 0.8621 0.5323 0.9349 0.9891 0.9878 0.9121 0.8398 0.9121 0.9789 0.8958 0.809 0.764 0.971 0.817 0.9298 0.8846 0.7308 0.9412 0.641 0.988 0.9679 0.9942 0.9767 0.8951 0.8339];
SOSI_RO_WIR = [0.7333 0.8636 0.9987 0.6667 0.7143 0.5354 0.4393 0.45 0.7391 0.8088 0.575 0.7273 0.6739 0.6792 0.625 0.6441 0.8681 0.6087 0.6053 0.625 0.6304 0.4898 0.6533 0.726 0.6825 0.5942 0.899 0.8395 0.625 0.5417];
SOSI_F_WIR = [0.7727 0.5741 0.5778 0.6094 0.6389 0.58 0.4286 0.4444 0.6842 0.5 0.0013 0.5135 0.5303 0.541 0.5 0.5067 0.6279 0.7037 0.4667 0.5118 0.5795 0.5882 0.4074 0.7368 0.5238 0.4468 0.5977 0.6 0.625 0.0013];

%% correlation between accuracy and d'

expName = 'SOSI';

% RS and Source d'
xdata = SOSI_RS_WIR(~ismember(subjects,badSub));
xstr = 'RS accuracy';
ydata = SOSI_SOURCE_DP(~ismember(subjects,badSub));
ystr = 'Source d''';

[rho,p] = corr(xdata',ydata');
fprintf('%s %s and %s:\tr=%.4f, p=%.9f',expName,xstr,ystr,rho,p);
if p < .05
  fprintf(' *\n');
else
  fprintf('\n');
end
figure;
plot(xdata,ydata,'ko');
hold on
axis([0 1 0 ceil(max(ydata))]);
title(sprintf('%s, %s and %s, r=%.3f, p=%.4f',expName,xstr,ystr,rho,p));
xlabel(sprintf('%s',xstr));
ylabel(sprintf('%s',ystr));

MB = polyfit(xdata,ydata,1);
X = [min(xdata) max(xdata)];
Ms = ones(1, 2)*MB(1);
Bs = ones(1,2)*MB(2);
Y = Ms.*X + Bs;
plot(X,Y,'r');
hold off

% RS and Item d'
xdata = SOSI_RS_WIR(~ismember(subjects,badSub));
xstr = 'RS accuracy';
ydata = SOSI_ITEM_DP(~ismember(subjects,badSub));
ystr = 'Item d''';

[rho,p] = corr(xdata',ydata');
fprintf('%s %s and %s:\tr=%.4f, p=%.9f',expName,xstr,ystr,rho,p);
if p < .05
  fprintf(' *\n');
else
  fprintf('\n');
end
figure;
plot(xdata,ydata,'ko');
hold on
axis([0 1 0 ceil(max(ydata))]);
title(sprintf('%s, %s and %s, r=%.3f, p=%.4f',expName,xstr,ystr,rho,p));
xlabel(sprintf('%s',xstr));
ylabel(sprintf('%s',ystr));

MB = polyfit(xdata,ydata,1);
X = [min(xdata) max(xdata)];
Ms = ones(1, 2)*MB(1);
Bs = ones(1,2)*MB(2);
Y = Ms.*X + Bs;
plot(X,Y,'r');
hold off

% RO and Source d'
xdata = SOSI_RO_WIR(~ismember(subjects,badSub));
xstr = 'RO accuracy';
ydata = SOSI_SOURCE_DP(~ismember(subjects,badSub));
ystr = 'Source d''';

[rho,p] = corr(xdata',ydata');
fprintf('%s %s and %s:\tr=%.4f, p=%.9f',expName,xstr,ystr,rho,p);
if p < .05
  fprintf(' *\n');
else
  fprintf('\n');
end
figure;
plot(xdata,ydata,'ko');
hold on
axis([0 1 0 ceil(max(ydata))]);
title(sprintf('%s, %s and %s, r=%.3f, p=%.4f',expName,xstr,ystr,rho,p));
xlabel(sprintf('%s',xstr));
ylabel(sprintf('%s',ystr));

MB = polyfit(xdata,ydata,1);
X = [min(xdata) max(xdata)];
Ms = ones(1, 2)*MB(1);
Bs = ones(1,2)*MB(2);
Y = Ms.*X + Bs;
plot(X,Y,'r');
hold off

% RO and Item d'
xdata = SOSI_RO_WIR(~ismember(subjects,badSub));
xstr = 'RO accuracy';
ydata = SOSI_ITEM_DP(~ismember(subjects,badSub));
ystr = 'Item d''';

[rho,p] = corr(xdata',ydata');
fprintf('%s %s and %s:\tr=%.4f, p=%.9f',expName,xstr,ystr,rho,p);
if p < .05
  fprintf(' *\n');
else
  fprintf('\n');
end
figure;
plot(xdata,ydata,'ko');
hold on
axis([0 1 0 ceil(max(ydata))]);
title(sprintf('%s, %s and %s, r=%.3f, p=%.4f',expName,xstr,ystr,rho,p));
xlabel(sprintf('%s',xstr));
ylabel(sprintf('%s',ystr));

MB = polyfit(xdata,ydata,1);
X = [min(xdata) max(xdata)];
Ms = ones(1, 2)*MB(1);
Bs = ones(1,2)*MB(2);
Y = Ms.*X + Bs;
plot(X,Y,'r');
hold off

% F and Source d'
xdata = SOSI_F_WIR(~ismember(subjects,badSub));
xstr = 'F accuracy';
ydata = SOSI_SOURCE_DP(~ismember(subjects,badSub));
ystr = 'Source d''';

[rho,p] = corr(xdata',ydata');
fprintf('%s %s and %s:\tr=%.4f, p=%.9f',expName,xstr,ystr,rho,p);
if p < .05
  fprintf(' *\n');
else
  fprintf('\n');
end
figure;
plot(xdata,ydata,'ko');
hold on
axis([0 1 0 ceil(max(ydata))]);
title(sprintf('%s, %s and %s, r=%.3f, p=%.4f',expName,xstr,ystr,rho,p));
xlabel(sprintf('%s',xstr));
ylabel(sprintf('%s',ystr));

MB = polyfit(xdata,ydata,1);
X = [min(xdata) max(xdata)];
Ms = ones(1, 2)*MB(1);
Bs = ones(1,2)*MB(2);
Y = Ms.*X + Bs;
plot(X,Y,'r');
hold off

% F and Item d'
xdata = SOSI_F_WIR(~ismember(subjects,badSub));
xstr = 'F accuracy';
ydata = SOSI_ITEM_DP(~ismember(subjects,badSub));
ystr = 'Item d''';

[rho,p] = corr(xdata',ydata');
fprintf('%s %s and %s:\tr=%.4f, p=%.9f',expName,xstr,ystr,rho,p);
if p < .05
  fprintf(' *\n');
else
  fprintf('\n');
end
figure;
plot(xdata,ydata,'ko');
hold on
axis([0 1 0 ceil(max(ydata))]);
title(sprintf('%s, %s and %s, r=%.3f, p=%.4f',expName,xstr,ystr,rho,p));
xlabel(sprintf('%s',xstr));
ylabel(sprintf('%s',ystr));

MB = polyfit(xdata,ydata,1);
X = [min(xdata) max(xdata)];
Ms = ones(1, 2)*MB(1);
Bs = ones(1,2)*MB(2);
Y = Ms.*X + Bs;
plot(X,Y,'r');
hold off


%% plot
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
%bw_colormap = [.5 .5 .5];
bw_colormap = [.7 .7 .7];
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
%print(gcf,'-dpng','~/Desktop/SOSI_RS_RO_F_accuracy');
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
