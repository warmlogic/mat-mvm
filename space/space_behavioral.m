% space behavioral analyses

%% load the analysis details

subDir = '';
dataDir = fullfile('SPACE','EEG','Sessions','ftpp',subDir);
% Possible locations of the data files (dataroot)
serverDir = fullfile(filesep,'Volumes','curranlab','Data');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dreamDir = fullfile(filesep,'data','projects','curranlab');
localDir = fullfile(getenv('HOME'),'data');

% pick the right dataroot
if exist('serverDir','var') && exist(serverDir,'dir')
  dataroot = serverDir;
  %runLocally = 1;
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
  %runLocally = 1;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
  %runLocally = 0;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

% procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');

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
  'SPACE037';
  };

% only one cell, with all session names
sesNames = {'session_1'};

allowRecallSynonyms = true;

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%% load the behavioral data

behfile = fullfile(getenv('HOME'),'data','SPACE','Behavioral','Sessions','SPACE_behav_results.mat');

load(behfile);

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {{'SPACE001','SPACE017','SPACE019'}};

% % exclude subjects with low event counts
% [exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

% exper.badSub = zeros(size(subjects));
exper.badSub = ismember(subjects,exper.badBehSub{1});

%% ttest stuff

alpha = 0.05;
tails = 'both';

%% recog

ses = 'oneDay';
phase = 'cued_recall';
test = 'recog';
measure = 'recog_hr';
% measure = 'recog_dp';
% measure = 'recog_rt';
% measure = 'recog_rt_hit';
% measure = 'recog_rt_miss';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);
data1_str = sprintf('%s %s %s',manip,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);
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
measure = 'recall_hr';
% measure = 'recall_rt';
% measure = 'recall_rt_hit';
% measure = 'recall_rt_miss';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);
data1_str = sprintf('%s %s %s',manip,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);
data2_str = sprintf('%s %s %s',manip,test,measure);

d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.10f\n', ...
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
measure = 'recall_hr';
% measure = 'recall_rt';
% measure = 'recall_rt_hit';
% measure = 'recall_rt_miss';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
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
measure = 'recall_hr';
% measure = 'recall_rt';
% measure = 'recall_rt_hit';
% measure = 'recall_rt_miss';

manip = 'spaced';

data1 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

manip = 'massed';

data2 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
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
measure = 'recall_hr';
% measure = 'recall_rt';
% measure = 'recall_rt_hit';
% measure = 'recall_rt_miss';

manip = 'spaced';

categ = 'Faces';

data1 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

categ = 'HouseInside';

data2 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
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
measure = 'recall_hr';
% measure = 'recall_rt';
% measure = 'recall_rt_hit';
% measure = 'recall_rt_miss';

manip = 'massed';

categ = 'Faces';

data1 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
data1_str = sprintf('%s %s %s %s',manip,categ,test,measure);

categ = 'HouseInside';

data2 = results.(ses).(phase).(manip).(categ).(test).(measure)(~exper.badSub);
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

