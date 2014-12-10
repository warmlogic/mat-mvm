% space behavioral analyses

%% load the analysis details

expName = 'SPACE2';

% beh_dir = 'behavioral_pilot';
beh_dir = 'Behavioral';

subDir = '';
behDir = fullfile(expName,beh_dir,'Sessions',subDir);
eegDir = fullfile(expName,'EEG','Sessions','ftpp',subDir);
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

behDir = fullfile(dataroot,behDir);

% procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
% procDir = fullfile(dataroot,eegDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');
% procDir = fullfile(dataroot,eegDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');

subjects = {
  'SPACE2001'
  'SPACE2002' % really noisy EEG, finished both sessions
  'SPACE2003' % DNF session 2
  'SPACE2004'
  'SPACE2005'
  'SPACE2006'
  'SPACE2007' % terrible performance
  'SPACE2008'
  'SPACE2009' % DNF session 2
  'SPACE2010'
  'SPACE2011';
  'SPACE2012';
  %'SPACE2013'; % didn't record EEG, stopped session 1 in middle
  'SPACE2014';
  'SPACE2015';
  'SPACE2016';
  'SPACE2017'; % terrible performance
  'SPACE2018';
  'SPACE2019';
  'SPACE2020'; % DNF session 2
  'SPACE2021';
  'SPACE2022';
  'SPACE2023'; % no ses2
  'SPACE2024'; % no ses2
  'SPACE2025'; % terrible performance
  'SPACE2026';
  'SPACE2027'; % really noisy EEG, DNF session 2
  'SPACE2028'; % no ses2
  'SPACE2029';
  'SPACE2030'; % no ses2
  'SPACE2031'; % no ses2
  'SPACE2032'; % no ses2
  };

% subjects = {
%   'SPACE2001'
%   'SPACE2002'
%   'SPACE2003'
%   'SPACE2004'
%   'SPACE2005'
%   'SPACE2006'
%   'SPACE2007'
%   'SPACE2008'
%   'SPACE2009'
%   'SPACE2010'
%   'SPACE2011'
%   'SPACE2012'
%   'SPACE2013'
%   'SPACE2014'
%   'SPACE2015' % really poor performance
%   'SPACE2016'
%   'SPACE2017' % really poor performance
%   'SPACE2018'
%   'SPACE2019'
%   'SPACE2020'
%   'SPACE2021';
%   'SPACE2022';
%   'SPACE2023';
%   'SPACE2024';
%     'SPACE2025';
%     'SPACE2026';
%     'SPACE2027';
%     'SPACE2028';
%     'SPACE2029';
%     'SPACE2029-2';
%     'SPACE2030';
%     'SPACE2031';
%     'SPACE2032';
%     'SPACE2033';
%     'SPACE2034';
%     'SPACE2035';
%     'SPACE2036';
%   };

% only one cell, with all session names
% sesNames = {'session_1'};
% 
% allowRecallSynonyms = true;
% 
% % replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
% replaceDataroot = true;
% 
% [exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

saveFigs = true;
if saveFigs
  %figsDir = fullfile(dataroot,dirs.behDir,'figs');
  figsDir = fullfile(behDir,'figs');
  if ~exist(figsDir,'dir')
    mkdir(figsDir);
  end
  %figFormat = '-djpg90';
  figFormat = '-dpng';
  %figFormat = '-deps2';
  %figFormat = '-depsc2';
  
  figRes = '-r150';
end

%% collapse phases?

collapsePhases = true;

if collapsePhases
  collapseStr = '_collapsed';
else
  collapseStr = '';
end

%% split into quantile divisions?

nDivisions = 1;
% nDivisions = 2;
% nDivisions = 3;
% nDivisions = 4;

if nDivisions > 1
  quantStr = sprintf('_%dquantileDiv',nDivisions);
else
  quantStr = '';
end

%% load the behavioral data

%resultsFile = fullfile(dataroot,dirs.behDir,sprintf('%s_behav_results%s%s.mat',expName,quantStr,collapseStr));
resultsFile = fullfile(behDir,sprintf('%s_behav_results%s%s.mat',expName,quantStr,collapseStr));

fprintf('Loading %s...',resultsFile);
load(resultsFile);
fprintf('Done.\n');

%% decide who to kick out based on trial counts

exper.subjects = subjects;

% Subjects with bad behavior
exper.badBehSub = {{}};

% exper.badBehSub = {{
%   'SPACE2001' % 5??? blocks
%   'SPACE2002'
%   'SPACE2003'
%   'SPACE2004'
%   'SPACE2005'
%   'SPACE2006'
%   'SPACE2007'
%   'SPACE2008'
%   'SPACE2009'
%   'SPACE2010'
%   'SPACE2011' % 7??? blocks
%   'SPACE2012'
%   'SPACE2013'
%   'SPACE2014'
%   'SPACE2015' % really poor performance
%   'SPACE2016'
%   'SPACE2017' % really poor performance
%   'SPACE2018'
%   'SPACE2019'
%   'SPACE2020'
%   'SPACE2021';
%   'SPACE2022';
%   'SPACE2023';
%   'SPACE2024';
%     'SPACE2025';
%     'SPACE2026';
% %     'SPACE2027'; % started with 9 blocks here
% %     'SPACE2028';
% %     'SPACE2029';
% %     'SPACE2029-2';
% %     'SPACE2030';
% %     'SPACE2031';
% %     'SPACE2032';
% %     'SPACE2033';
% %     'SPACE2034';
% %     'SPACE2035';
% %     'SPACE2036';
% }};

% % exclude subjects with low event counts
% [exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

% exper.badSub = zeros(size(subjects));
exper.badSub = ismember(subjects,exper.badBehSub{1});

%% ttest stuff

alpha = 0.05;
tails = 'both';

% %% expo - rating
% 
% ses = 'day1';
% phase = 'expo';
% test = 'rating';
% measure = 'rating_resp';
% 
% manip = 'Faces';
% data1 = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);
% data1_str = sprintf('%s %s %s',manip,test,measure);
% manip = 'HouseInside';
% data2 = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);
% data2_str = sprintf('%s %s %s',manip,test,measure);
% 
% d = mm_effect_size('within',data1,data2);
% [h, p, ci, stats] = ttest(data1,data2,alpha,tails);
% 
% fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.11f\n', ...
%   data1_str, ...
%   mean(data1), ...
%   std(data1) / sqrt(length(data1)), ...
%   data2_str, ...
%   mean(data2), ...
%   std(data2) / sqrt(length(data2)), ...
%   stats.df, ...
%   stats.tstat, ...
%   d, ...
%   std(data1 - data2),...
%   std(data1 - data2) / sqrt(length(data1)),...
%   p);

%% ANOVA - test - recall - spacing (spaced/massed) X category (Faces/HouseInside)

ses = 'day1';
% ses = 'day2';
phase = 'cued_recall_only';
test = 'recall';
measure = 'recall_hr';
% measure = 'recall_rt';
% measure = 'recall_rt_hit';
% measure = 'recall_rt_miss';
% measure = 'recall_nHit';

spacings = {'once','massed','lag2','lag12','lag32'};
stimCats = {'Faces', 'HouseInside'};

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~exper.badSub(sub,:)
    %for ses = 1:length(exper.sesStr)
    theseData = [];
    
    for sp = 1:length(spacings)
      for st = 1:length(stimCats)
        theseData = cat(2,theseData,results.(ses).(phase).(spacings{sp}).(stimCats{st}).(test).(measure)(sub));
      end
    end
    %end
    anovaData = cat(1,anovaData,theseData);
  end
end

% levelnames = {{'img','word'}, {'rc', 'fo'}, {'spac','mass'}, latStr};
% varnames = {'stimType','subseqMem','spacing','time'};
% O = teg_repeated_measures_ANOVA(anovaData, [2 2 2 length(latInd(1):latInd(2))], varnames,[],[],[],[],[],[],levelnames);

levelnames = {spacings, stimCats};
varnames = {'spacing','img_cat'};
O = teg_repeated_measures_ANOVA(anovaData, [length(spacings),length(stimCats)], varnames,[],[],[],[],[],[],levelnames);



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


%% Plot: spaced vs massed

sesName = 'oneDay';
% test = 'recog';
test = 'recall';

if collapsePhases
  phases = {'cued_recall'};
else
  phases = {'cued_recall_1','cued_recall_2','cued_recall_3','cued_recall_4','cued_recall_5','cued_recall_6'};
end

data = struct;
data.spaced = nan(sum(~exper.badSub),length(phases));
data.massed = nan(sum(~exper.badSub),length(phases));

% % d' (recog only)
% dataMeasure = sprintf('%s_dp',test);ylimits = [0 4];
% dataLabel = 'd''';

% Hit rate / Accuracy (recog and recall)
dataMeasure = sprintf('%s_hr',test);ylimits = [0 1];
if strcmp(test,'recog')
  dataLabel = 'Hit Rate';
elseif strcmp(test,'recall')
  dataLabel = 'Accuracy';
end

% % Response times (recog and recall)
% % dataMeasure = sprintf('%s_rt',test);
% % dataLabel = 'Response Time (ms)';
% dataMeasure = sprintf('%s_rt_hit',test);
% dataLabel = 'Response Time: Hits (ms)';
% % dataMeasure = sprintf('%s_rt_miss',test);
% % dataLabel = 'Response Time: Misses (ms)';
% if strcmp(test,'recog')
%   ylimits = [500 2500];
% elseif strcmp(test,'recall')
%   ylimits = [500 4000];
% end

tpCounter = 0;
for p = 1:length(phases)
  if isfield(results.(sesName),phases{p})
    tpCounter = tpCounter + 1;
    
    manip = 'spaced';
    data.spaced(:,tpCounter) = results.(sesName).(phases{p}).(manip).(test).(dataMeasure)(~exper.badSub);
    
    manip = 'massed';
    data.massed(:,tpCounter) = results.(sesName).(phases{p}).(manip).(test).(dataMeasure)(~exper.badSub);
  end
end

spaced_mean = nanmean(data.spaced,1);
spaced_sem = nanstd(data.spaced,1) ./ sqrt(sum(~isnan(data.spaced)));

massed_mean = nanmean(data.massed,1);
massed_sem = nanstd(data.massed,1) ./ sqrt(sum(~isnan(data.massed)));

if collapsePhases
  %bw_title = sprintf('Test phases (collapsed): %s',strrep(dataMeasure,'_','\_'));
  if strcmp(test,'recog')
    bw_title = 'Test: Recognition (collapsed)';
  elseif strcmp(test,'recall')
    bw_title = 'Test: Cued Recall (collapsed)';
  end
else
  %bw_title = sprintf('Test phases: %s',strrep(dataMeasure,'_','\_'));
  if strcmp(test,'recog')
    bw_title = 'Test: Recognition';
  elseif strcmp(test,'recall')
    bw_title = 'Test: Cued Recall';
  end
end

bw_legend = {'Spaced','Massed'};
if length(phases) > 1
  bw_width = 1;
  bw_groupnames = 1:length(phases);
  bw_xlabel = 'Block number';
else
  bw_width = 0.75;
  bw_groupnames = [];
  bw_xlabel = [];
end
bw_ylabel = dataLabel;
bw_colormap = 'gray';
% bw_colormap = 'linspecer';
bw_data = [spaced_mean;massed_mean]';
bw_errors = [spaced_sem;massed_sem]';
bw_legend_type = 'plot';
% bw_legend_type = 'axis';
h = barweb(bw_data,bw_errors,bw_width,bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend,[],bw_legend_type);
if strcmp(bw_legend_type,'plot')
  set(h.legend,'Location','NorthEast');
end
if length(phases) > 1
  axis([0.5 (length(phases)+0.5) ylimits(1) ylimits(2)]);
else
  axis([0 (length(phases) + 1) ylimits(1) ylimits(2)]);
end
axis square
% ylim([ylimits(1) ylimits(2)]);

publishfig(gcf,0);

if saveFigs
  if collapsePhases
    collapseStr = '_collapsed';
  else
    collapseStr = '';
  end
  print(gcf,figFormat,figRes,fullfile(figsDir,sprintf('test_%s%s',dataMeasure,collapseStr)));
end
