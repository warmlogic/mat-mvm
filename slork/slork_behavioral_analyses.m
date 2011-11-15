% Can I also calculate F from the IRK equation for experiments 1 and 4 using
% confidence? What would Wixted say?

figFontName = 'Helvetica';

subjects = {
  'SLORK002'...
  'SLORK003'...
  'SLORK004'...
  'SLORK005'...
  'SLORK006'...
  'SLORK007'...
  'SLORK008'...
  'SLORK009'...
  'SLORK010'...
  'SLORK011'...
  'SLORK012'...
  'SLORK013'...
  'SLORK014'...
  'SLORK015'...
  'SLORK016'...
  'SLORK017'...
  'SLORK018'...
  'SLORK019'...
  'SLORK020'...
  'SLORK022'...
  'SLORK023'...
  'SLORK024'...
  'SLORK025'...
  'SLORK026'...
  'SLORK027'...
  'SLORK029'...
  'SLORK030'...
  'SLORK031'...
  'SLORK032'...
  'SLORK033'...
  };

badSub = {};

% source accuracy artRej
%badSub = {'SLORK009', 'SLORK015', 'SLORK016', 'SLORK024', 'SLORK026', 'SLORK032', 'SLORK033'};
% RKN artRej
%badSub = {'SLORK003','SLORK009', 'SLORK013', 'SLORK015', 'SLORK018', 'SLORK029'};

% ROKN artRej
%badSub = {'SLORK003','SLORK009', 'SLORK013', 'SLORK014', 'SLORK015', 'SLORK017', 'SLORK018', 'SLORK023', 'SLORK029'};

% % Familiarity values calculated using IRK: F = "know"/(1-R)
% S_RH_SrcCor_fam = [0.22 0.06 0.45 0.18 0.23 0.14 0.08 0.10 0.17 0.28 0.12 0.00 0.28 0.05 0.21 0.42 0.02 0.21 0.14 0.13 0.33 0.13 0.09 0.76 0.14 0.00 0.08 0.14 0.10 0.25];
% S_RH_SrcInc_fam = [0.08 0.02 0.10 0.06 0.08 0.15 0.03 0.01 0.09 0.18 0.08 0.00 0.19 0.00 0.07 0.28 0.00 0.11 0.04 0.07 0.08 0.03 0.10 0.05 0.05 0.01 0.11 0.09 0.11 0.07];
% Q_RH_SrcCor_fam = [0.06 0.01 0.29 0.12 0.16 0.17 0.09 0.03 0.10 0.28 0.10 0.00 0.27 0.00 0.27 0.50 0.05 0.25 0.08 0.06 0.16 0.14 0.07 0.49 0.10 0.00 0.14 0.08 0.11 0.16];
% Q_RH_SrcInc_fam = [0.13 0.01 0.13 0.10 0.12 0.26 0.06 0.02 0.11 0.21 0.04 0.00 0.17 0.00 0.06 0.34 0.04 0.13 0.08 0.08 0.03 0.02 0.11 0.11 0.06 0.00 0.16 0.02 0.04 0.10];
% 
% nSub = length(Q_RH_SrcCor_fam);
% 
% rates = {S_RH_SrcCor_fam S_RH_SrcInc_fam;...
%   Q_RH_SrcCor_fam Q_RH_SrcInc_fam};
% 
% cfg.task = {'Side','Question'};
% cfg.accur = {'Correct','Incorrect'};
% 
% cfg.alpha = .05;
% cfg.showtable = 1;
% cfg.printTable_tex = 1;
% if length(cfg.task) > 2 || length(cfg.accur) > 2
%   cfg.calcGGHF = 1;
% else
%   cfg.calcGGHF = 0;
% end
% 
% fprintf('====================== RMAOV2_mod =========================\n');
% fprintf('IV1: Task (%d;%s), IV2: Source Accuracy (%d;%s)\n',length(cfg.task),sprintf(repmat(' %s',1,length(cfg.task)),cfg.task{:}),length(cfg.accur),sprintf(repmat(' %s',1,length(cfg.accur)),cfg.accur{:}));
% 
% % col 1 = dependent variable (rate)
% % col 2 = independent variable 1 (task: side or question)
% % col 3 = independent variable 2 (accuracy: correct or incorrect)
% % col 4 = subject number
% 
% anovamat = [];
% for t = 1:length(cfg.task)
%   for a = 1:length(cfg.accur)
%     subCount = 0;
%     for sub = 1:nSub
%       if ~ismember(subjects(sub),badSub)
%         subCount = subCount + 1;
%         anovamat = [anovamat; rates{t,a}(sub) t a subCount];
%       end
%     end
%   end
% end
% p_vals = RMAOV2_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

% % allSub
% 
% % There was a main effect of source accuracy [F(1,29)=15.464, MSE=0.012,
% % p<.001]. Familiarity contributed more to responses with correct source
% % retrieval than to those with incorrect source retrieval.
% 
% % There was a Task X Source Accuracy interaction [F(1,29)=10.963, MSE=0.002,
% % p<.01]. Familiarity contributed significantly more to correct source
% % retrieval in the side task than in the question task.
% 
% 
% % Excluding badSub
% 
% % There was a main effect of source accuracy [F(1,22)=15.369, MSE=0.005,
% % p<.001]. Familiarity contributed more to responses with correct source
% % retrieval than to those with incorrect source retrieval.
% 
% % There was a Task X Source Accuracy interaction [F(1,22)=13.256, MSE=0.001,
% % p<.005]. Familiarity contributed significantly more to correct source
% % retrieval in the side task than in the question task.

% RH_SrcCor_fam = mean([S_RH_SrcCor_fam;Q_RH_SrcCor_fam],1);
% RH_SrcInc_fam = mean([S_RH_SrcInc_fam;Q_RH_SrcInc_fam],1);
% 
% [h,p,ci,stats] = ttest(RH_SrcCor_fam(~ismember(subjects,badSub)),RH_SrcInc_fam(~ismember(subjects,badSub)),0.05,'both');
% fprintf('SrcCor (M=%.3f) vs. SrcInc (M=%.3f): t(%d)=%.3f, p=%.5f\n',mean(RH_SrcCor_fam(~ismember(subjects,badSub))),mean(RH_SrcInc_fam(~ismember(subjects,badSub))),stats.df,stats.tstat,p);
% 
% [h,p,ci,stats] = ttest(S_RH_SrcCor_fam(~ismember(subjects,badSub)),Q_RH_SrcCor_fam(~ismember(subjects,badSub)),0.05,'both');
% fprintf('Side SC (M=%.3f) vs. Question SC (M=%.3f): t(%d)=%.3f, p=%.5f\n',mean(S_RH_SrcCor_fam(~ismember(subjects,badSub))),mean(Q_RH_SrcCor_fam(~ismember(subjects,badSub))),stats.df,stats.tstat,p);
% %[h,p,ci,stats] = ttest(S_RH_SrcInc_fam(~ismember(subjects,badSub)),Q_RH_SrcInc_fam(~ismember(subjects,badSub)),0.05,'both');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% within-response two-way ANOVA for all three responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_RS_WIR = [0.8765, 0.8621, 0.9474, 0.9067, 0.8545, 0.75, 0.8889, 0.8519, 0.6324, 0.7273, 0.8525, 0.8571, 0.85, 0.964, 0.85, 0.9091, 0.8191, 0.8, 0.7708, 0.9, 0.837, 0.9123, 0.45, 0.9965, 0.8095, 0.7042, 0.425, 0.7, 0.6061, 0.4815];
S_RO_WIR = [0.5714, 0.6957, 0.7692, 0.7045, 0.7455, 0.5205, 0.64, 0.7895, 0.4815, 0.6, 0.5938, 0.5714, 0.0035, 0.7826, 0.7679, 0.8333, 0.6667, 0.6842, 0.6667, 0.5789, 0.9965, 0.8824, 0.6316, 0.9, 0.5714, 0.3333, 0.4762, 0.5, 0.4898, 0.5833];
S_F_WIR = [0.5833, 0.5, 0.72, 0.5333, 0.5652, 0.4667, 0.6, 0.8889, 0.625, 0.5893, 0.4737, 0.0035, 0.5098, 0.9965, 0.6667, 0.5556, 0.9965, 0.6216, 0.6875, 0.4375, 0.6774, 0.7143, 0.4545, 0.9, 0.6471, 0.0035, 0.4286, 0.5238, 0.4444, 0.7778];

Q_RS_WIR = [0.7711, 0.7558, 0.9394, 0.7727, 0.7636, 0.5909, 0.8125, 0.8167, 0.6415, 0.875, 0.8333, 0.8125, 0.7368, 0.8348, 0.8837, 0.8182, 0.7941, 0.8462, 0.6604, 0.7869, 0.6923, 0.8571, 0.5323, 0.7432, 0.8235, 0.5781, 0.6818, 0.7952, 0.6184, 0.7222];
Q_RO_WIR = [0.2857, 0.3846, 0.5, 0.6111, 0.6087, 0.4062, 0.7073, 0.4444, 0.4, 0.6538, 0.6333, 0.3846, 0.5, 0.4545, 0.9394, 0.0035, 0.6667, 0.75, 0.5, 0.5455, 0.6667, 0.6667, 0.2857, 0.7143, 0.5333, 0.0035, 0.375, 0.5455, 0.5532, 0.72];
Q_F_WIR = [0.25, 0.5, 0.575, 0.45, 0.4615, 0.4186, 0.5333, 0.5, 0.4737, 0.5484, 0.6154, 0.0035, 0.551, 0.0035, 0.7143, 0.5826, 0.375, 0.5909, 0.4667, 0.3571, 0.7857, 0.7692, 0.3889, 0.7455, 0.5, 0.0035, 0.4615, 0.7143, 0.6667, 0.5517];

nSub = length(S_RS_WIR);

respaccur = {S_RS_WIR, S_RO_WIR, S_F_WIR;...
  Q_RS_WIR, Q_RO_WIR, Q_F_WIR};

%respaccur = {S_RS_WIR, S_F_WIR;...
%  Q_RS_WIR, Q_F_WIR};

cfg.task = {'Side','Question'};
cfg.resp = {'RSrc','RO','F'};
%cfg.resp = {'RSrc','F'};

cfg.alpha = .05;
cfg.showtable = 1;
cfg.printTable_tex = 1;
if length(cfg.task) > 2 || length(cfg.resp) > 2
  cfg.calcGGHF = 1;
else
  cfg.calcGGHF = 0;
end

fprintf('====================== RMAOV2_mod =========================\n');
fprintf('IV1: Task (%d;%s), IV2: Response (%d;%s)\n',length(cfg.task),sprintf(repmat(' %s',1,length(cfg.task)),cfg.task{:}),length(cfg.resp),sprintf(repmat(' %s',1,length(cfg.resp)),cfg.resp{:}));

% col 1 = dependent variable (rate)
% col 2 = independent variable 1 (task: side or question)
% col 3 = independent variable 2 (response: RSrc, RO, F)
% col 4 = subject number

anovamat = [];
for t = 1:length(cfg.task)
  for r = 1:length(cfg.resp)
    subCount = 0;
    for sub = 1:nSub
      if ~ismember(subjects(sub),badSub)
        subCount = subCount + 1;
        anovamat = [anovamat; respaccur{t,r}(sub) t r subCount];
      end
    end
  end
end
p_vals = RMAOV2_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

% IV1
S_WIR = mean([S_RS_WIR;S_RO_WIR;S_F_WIR],1);
Q_WIR = mean([Q_RS_WIR;Q_RO_WIR;Q_F_WIR],1);
[h,p,ci,stats] = ttest(S_WIR(~ismember(subjects,badSub)),Q_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('S (M=%.3f) vs. Q (M=%.3f): t(%d)=%.3f, p=%.10f\n',mean(S_WIR(~ismember(subjects,badSub))),mean(Q_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);


% IV2
RS_WIR = mean([S_RS_WIR;Q_RS_WIR],1);
RO_WIR = mean([S_RO_WIR;Q_RO_WIR],1);
F_WIR = mean([S_F_WIR;Q_F_WIR],1);

[h,p,ci,stats] = ttest(RS_WIR(~ismember(subjects,badSub)),RO_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('RS (M=%.3f) vs. RO (M=%.3f): t(%d)=%.3f, p=%.10f\n',mean(RS_WIR(~ismember(subjects,badSub))),mean(RO_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(RS_WIR(~ismember(subjects,badSub)),F_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('RS (M=%.3f) vs. F (M=%.3f): t(%d)=%.3f, p=%.10f\n',mean(RS_WIR(~ismember(subjects,badSub))),mean(F_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(RO_WIR(~ismember(subjects,badSub)),F_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('RO (M=%.3f) vs. F (M=%.3f): t(%d)=%.3f, p=%.10f\n',mean(RO_WIR(~ismember(subjects,badSub))),mean(F_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);


% IV1 x IV2
[h,p,ci,stats] = ttest(S_RS_WIR(~ismember(subjects,badSub)),Q_RS_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('S_RS (M=%.3f) vs. Q_RS (M=%.3f): t(%d)=%.3f, p=%.5f\n',mean(S_RS_WIR(~ismember(subjects,badSub))),mean(Q_RS_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(S_RO_WIR(~ismember(subjects,badSub)),Q_RO_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('S_RO (M=%.3f) vs. Q_RO (M=%.3f): t(%d)=%.3f, p=%.5f\n',mean(S_RO_WIR(~ismember(subjects,badSub))),mean(Q_RO_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(S_F_WIR(~ismember(subjects,badSub)),Q_F_WIR(~ismember(subjects,badSub)),0.05,'both');
fprintf('S_F (M=%.3f) vs. Q_F (M=%.3f): t(%d)=%.3f, p=%.5f\n',mean(S_F_WIR(~ismember(subjects,badSub))),mean(Q_F_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);

% vs chance
[h,p,ci,stats] = ttest(S_F_WIR(~ismember(subjects,badSub)),repmat(.5,size(S_F_WIR(~ismember(subjects,badSub)))),0.05,'both');
fprintf('S_F (M=%.3f) vs. chance: t(%d)=%.3f, p=%.5f\n',mean(S_F_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(Q_F_WIR(~ismember(subjects,badSub)),repmat(.5,size(S_F_WIR(~ismember(subjects,badSub)))),0.05,'both');
fprintf('Q_F (M=%.3f) vs. chance: t(%d)=%.3f, p=%.5f\n',mean(Q_F_WIR(~ismember(subjects,badSub))),stats.df,stats.tstat,p);


% plot
S_RS_avg = mean(S_RS_WIR(~ismember(subjects,badSub)));
Q_RS_avg = mean(Q_RS_WIR(~ismember(subjects,badSub)));
S_RO_avg = mean(S_RO_WIR(~ismember(subjects,badSub)));
Q_RO_avg = mean(Q_RO_WIR(~ismember(subjects,badSub)));
S_F_avg = mean(S_F_WIR(~ismember(subjects,badSub)));
Q_F_avg = mean(Q_F_WIR(~ismember(subjects,badSub)));

S_RS_sem = std(S_RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
Q_RS_sem = std(Q_RS_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
S_RO_sem = std(S_RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
Q_RO_sem = std(Q_RO_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
S_F_sem = std(S_F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));
Q_F_sem = std(Q_F_WIR(~ismember(subjects,badSub)))/sqrt(sum(~ismember(subjects,badSub)));

bw_groupnames = {'RS','RO','F'};
bw_title = 'Proportion of Source Correct responses';
bw_legend = {'Side','Question'};
bw_colormap = 'gray';
h = barweb([S_RS_avg, Q_RS_avg; S_RO_avg, Q_RO_avg; S_F_avg, Q_F_avg],[S_RS_sem, Q_RS_sem; S_RO_sem, Q_RO_sem; S_F_sem, Q_F_sem],[],bw_groupnames,bw_title,[],[],bw_colormap,[],bw_legend);
set(h.legend,'Location','NorthEast');
axis([0.5 3.5 0 1]);
publishfig(gcf,0,[],[],figFontName);
hold on
plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2); % horiz chance line
print(gcf,'-depsc2','~/Desktop/RS_RO_F_accuracy');
