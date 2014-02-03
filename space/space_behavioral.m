% space behavioral analyses

%% load the analysis details

procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';

subjects = {
  'SPACE001';
  'SPACE002';
  'SPACE003';
  'SPACE004';
  'SPACE005';
  'SPACE006';
  'SPACE007';
  %'SPACE008'; % didn't perform task correctly, didn't perform well
  'SPACE009';
  'SPACE010';
  'SPACE011';
  'SPACE012';
  'SPACE013';
  'SPACE014';
  'SPACE015';
  'SPACE016';
  'SPACE017'; % previous assessment: really noisy EEG, half of ICA components rejected
  'SPACE018';
  'SPACE019';
  'SPACE020';
  'SPACE021';
  'SPACE022';
  'SPACE027';
  'SPACE029';
  };

% only one cell, with all session names
sesNames = {'session_1'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = false;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {{}};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

%% load the behavioral data

behfile = fullfile(getenv('HOME'),'data','SPACE','Behavioral','Sessions','SPACE_behav_results.mat');

load(behfile);

%% ttest stuff

alpha = 0.05;
tails = 'both';

%% recog

ses = 'oneDay';
phase = 'cued_recall';
test = 'recog';
measure = 'hr';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data1_str = sprintf('%s %s %s',manip,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data2_str = sprintf('%s %s %s',manip,test,measure);
d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
  data1_str, ...
  mean(data1), ...
  std(data1) / sqrt(length(data1)), ...
  data2_str, ...
  mean(data2), ...
  std(data2) / sqrt(length(data2)), ...
  stats.df, ...
  stats.tstat, ...
  d, ...
  std(data1 - data2),...
  std(data1 - data2) / sqrt(length(data1)),...
  p);

%% recall

ses = 'oneDay';
phase = 'cued_recall';
test = 'recall';
measure = 'hr';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data1_str = sprintf('%s %s %s',manip,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data2_str = sprintf('%s %s %s',manip,test,measure);

d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
  data1_str, ...
  mean(data1), ...
  std(data1) / sqrt(length(data1)), ...
  data2_str, ...
  mean(data2), ...
  std(data2) / sqrt(length(data2)), ...
  stats.df, ...
  stats.tstat, ...
  d, ...
  std(data1 - data2),...
  std(data1 - data2) / sqrt(length(data1)),...
  p);

%% recall Faces

categ = 'Faces';

ses = 'oneDay';
phase = 'cued_recall';
test = 'recall';
measure = 'hr';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data2_str = sprintf('%s %s %s %s',manip,categ,test,measure);

d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
  data1_str, ...
  mean(data1), ...
  std(data1) / sqrt(length(data1)), ...
  data2_str, ...
  mean(data2), ...
  std(data2) / sqrt(length(data2)), ...
  stats.df, ...
  stats.tstat, ...
  d, ...
  std(data1 - data2),...
  std(data1 - data2) / sqrt(length(data1)),...
  p);

%% recall HouseInside

categ = 'HouseInside';

ses = 'oneDay';
phase = 'cued_recall';
test = 'recall';
measure = 'hr';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data2_str = sprintf('%s %s %s %s',manip,categ,test,measure);

d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
  data1_str, ...
  mean(data1), ...
  std(data1) / sqrt(length(data1)), ...
  data2_str, ...
  mean(data2), ...
  std(data2) / sqrt(length(data2)), ...
  stats.df, ...
  stats.tstat, ...
  d, ...
  std(data1 - data2),...
  std(data1 - data2) / sqrt(length(data1)),...
  p);

%% recall Faces vs HouseInside spaced

ses = 'oneDay';
phase = 'cued_recall';
test = 'recall';
measure = 'hr';

manip = 'spaced';

categ = 'Faces';

data1 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

categ = 'HouseInside';

data2 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data2_str = sprintf('%s %s %s %s',manip,categ,test,measure);

d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
  data1_str, ...
  mean(data1), ...
  std(data1) / sqrt(length(data1)), ...
  data2_str, ...
  mean(data2), ...
  std(data2) / sqrt(length(data2)), ...
  stats.df, ...
  stats.tstat, ...
  d, ...
  std(data1 - data2),...
  std(data1 - data2) / sqrt(length(data1)),...
  p);

%% recall Faces vs HouseInside massed

ses = 'oneDay';
phase = 'cued_recall';
test = 'recall';
measure = 'hr';

manip = 'massed';

categ = 'Faces';

data1 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

categ = 'HouseInside';

data2 = results.(ses).(phase).(manip).(categ).(test).(sprintf('%s_%s',test,measure))(~exper.badSub);
data2_str = sprintf('%s %s %s %s',manip,categ,test,measure);

d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
  data1_str, ...
  mean(data1), ...
  std(data1) / sqrt(length(data1)), ...
  data2_str, ...
  mean(data2), ...
  std(data2) / sqrt(length(data2)), ...
  stats.df, ...
  stats.tstat, ...
  d, ...
  std(data1 - data2),...
  std(data1 - data2) / sqrt(length(data1)),...
  p);

