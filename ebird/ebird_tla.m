%% load the analysis details

expName = 'EBIRD';

subDir = '';
dataDir = fullfile(expName,'EEG','Sessions','ftpp',subDir);
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

procDir = fullfile(dataroot,dataDir,'ft_data/data_art_nsClassic_ftAuto/tla');

subjects = {
  %'EBIRD049'; % Pilot. (due to short ses1 match, missing ses2 name)
  %'EBIRD002'; % Pilot. (due to short ses1 match, missing ses2 name)
  %'EBIRD003'; % Pilot. (due to missing ses7 name) - NB: LAST PILOT TO BE REPLACED
  %'EBIRD004'; % DNF. Dropout. Last session: 8.
  'EBIRD005';
  %'EBIRD006'; % DNF. Dropout. Last session: 2.
  'EBIRD007';
  'EBIRD008';
  'EBIRD009';
  'EBIRD010';
  'EBIRD011';
  'EBIRD012';
  %'EBIRD013'; % DNF. Dropout. Last session: 5. Lost session 6 in HD crash.
  %'EBIRD014'; % DNF. Rejected. Last session: 1.
  %'EBIRD015'; % DNF. Lost in HD crash.
  %'EBIRD016'; % DNF. Lost in HD crash.
  %'EBIRD017'; % DNF. Lost in HD crash.
  'EBIRD018';
  'EBIRD019';
  'EBIRD020';
  'EBIRD021';
  %'EBIRD022'; % DNF. Dropout. Last session: 8.
  %'EBIRD023'; % DNF. Dropout. Last session: 1.
  'EBIRD024';
  'EBIRD025';
  'EBIRD027';
  'EBIRD029';
  'EBIRD032';
  'EBIRD034';
  'EBIRD042';
  };

% only one cell, with all session names
sesNames = {'session_1','session_8','session_9'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = false;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);


% adFile = '/Users/matt/data/EBIRD/EEG/Sessions/ftpp/ft_data/match_stim_eq0_art_ftManual_ftICA/tla/analysisDetails.mat';
% % [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails(adFile,adFile5,true,true,true);
% 
% % % server_adFile = '';
% % if exist(server_adFile,'file')
% %   mm_mergeAnalysisDetails(adFile,server_adFile,true,false,false);
% % end
% 
% % [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
% [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,false);

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% list the event values to analyze; specific to each experiment

% this is useful for when there are multiple types of event values, for
% example, hits and CRs in two conditions. You don't have to enter anything
% if you just want all events from exper.eventValues together in a single
% cell because it will get set to {exper.eventValues}, but it needs to be a
% cell containing a cell of eventValue strings

% this is only used by mm_ft_checkCondComps to create pairwise combinations
% either within event types {'all_within_types'} or across all event types
% {'all_across_types'}; mm_ft_checkCondComps is called within subsequent
% analysis functions

% match

% trained
% sameSpecies



ana.eventValues = repmat({{'match_stim'}},1,length(exper.sessions));
ana.eventValuesSplit = repmat({...
  {{'stim1_basic_color', 'stim1_basic_g', 'stim1_basic_g_hi8', 'stim1_basic_g_lo8', 'stim1_basic_norm', ...
  'stim1_subord_color', 'stim1_subord_g', 'stim1_subord_g_hi8', 'stim1_subord_g_lo8', 'stim1_subord_norm', ...
  'stim2_basic_color', 'stim2_basic_g', 'stim2_basic_g_hi8', 'stim2_basic_g_lo8', 'stim2_basic_norm', ...
  'stim2_subord_color', 'stim2_subord_g', 'stim2_subord_g_hi8', 'stim2_subord_g_lo8', 'stim2_subord_norm'} ...
  }},1,length(exper.sessions));

ana.trl_expr = cell(1,length(exper.sessions));

for ses = 1:length(exper.sessions)
  if ismember(exper.sessions{ses},{'session_1','session_8','session_9'})
    ana.trl_expr{ses} = {...
      {...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & imgCond == 1 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & imgCond == 2 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & imgCond == 3 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & imgCond == 4 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & imgCond == 5 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & imgCond == 1 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & imgCond == 2 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & imgCond == 3 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & imgCond == 4 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & imgCond == 5 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & imgCond == 1 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & imgCond == 2 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & imgCond == 3 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & imgCond == 4 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & imgCond == 5 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & imgCond == 1 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & imgCond == 2 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & imgCond == 3 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & imgCond == 4 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & imgCond == 5 & response ~= 0 & rt < 2000 & acc == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
      }...
      };
  end
end

% % make sure ana.eventValues is set properly
% if ~iscell(ana.eventValues{1})
%   ana.eventValues = {ana.eventValues};
% end
% if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
%   ana.eventValues = {exper.eventValues};
% end

%% load in the subject data

% [data_tla,exper] = mm_ft_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');

keeptrials = false;
[data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');

% %% get rid of the bad channels
% 
% cfg = [];
% cfg.printRoi = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% [data_tla] = mm_rmBadChan(cfg,exper,ana,data_tla);

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {{},{},{}};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,10,[],'vert');
% [exper,ana] = mm_threshSubs(exper,ana,1);

%% Test plots to make sure data look ok

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.ylim = [-15 15];
cfg_ft.layout = ft_prepare_layout([],ana);
sub = 1;
ses = 1;
for i = 1:length(ana.eventValues{1})
  figure
  ft_multiplotER(cfg_ft,data_tla.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(strrep(ana.eventValues{1}{i},'_','-'));
end

% cfg_ft = [];
% cfg_ft.channel = {'E20'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% figure
% ft_singleplotER(cfg_ft,data_tla.(ana.eventValues{1}{1}).sub(1).ses(1).data,data_tla.(ana.eventValues{1}{2}).sub(1).ses(1).data);

% cfg_ft = [];
% cfg_ft.channel = {'E20'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% %figure
% %ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data,data_tla.(exper.eventValues{2}).sub(1).ses(1).data,data_tla.(exper.eventValues{3}).sub(1).ses(1).data);
% %legend(strrep(exper.eventValues,'_','-'),'Location','SouthEast');
% figure
% cfg_ft.graphcolor = 'b';
% ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% legend(strrep(exper.eventValues{1},'_','-'),'Location','SouthEast');
% hold on
% plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
% plot([0 0],[-5 5],'k--'); % vertical
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cfg_ft = [];
% cfg_ft.channel = {'all'};
% cfg_ft.latency = [.5 .8];
% data_topo = ft_timelockanalysis(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% cfg_ft = [];
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.interactive = 'yes';
% cfg_ft.colormap = 'jet';
% %cfg_ft.colormap = 'hot';
% cfg_ft.colorbar = 'yes';
% cfg_ft.xlim = 'maxmin';
% cfg_ft.layout = ft_prepare_layout([],ana);
% figure
% ft_topoplotER(cfg_ft,data_topo);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Remove the dof field. I'm not sure what degrees of freedom refers to, but the tutorial says to.
% for evVal = 1:length(exper.eventValues)
%   for sub = 1:length(exper.subjects)
%     for ses = 1:length(exper.sesStr)
%       if isfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof')
%         data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data = rmfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof');
%       end
%     end
%   end
% end

%% get the grand average

ga_tla = struct;

for ses = 1:length(exper.sesStr)
  % set up strings to put in grand average function
  cfg_ana = [];
  cfg_ana.is_ga = 0;
  %cfg_ana.conditions = ana.eventValues{ses};
  cfg_ana.data_str = 'data_tla';
  % cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);
  %cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
  
  cfg_ft = [];
  cfg_ft.keepindividual = 'no';
  %cfg_ft.parameter = 'trial';
  for typ = 1:length(ana.eventValues{ses})
    cfg_ana.conditions = ana.eventValues{ses}{typ};
    cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
    
    for evVal = 1:length(ana.eventValues{ses}{typ})
      %tic
      fprintf('Running ft_timelockgrandaverage on %s...',ana.eventValues{ses}{typ}{evVal});
      %ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{evVal}){ses}));
      ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
      fprintf('Done.\n');
      %toc
    end
  end
end

% % turn keeptrial data into average for statistical functions because proper
% % processing of dimord is currently broken
% 
% data_tla_avg = struct;
% 
% cfg = [];
% cfg.keeptrials = 'no';
% 
% for ses = 1:length(exper.sesStr)
%   for evVal = 1:length(ana.eventValues{ses})
%     for sub = 1:length(exper.subjects)
%       fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{evVal});
%       if isfield(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}).sub(sub).data,'avg')
%         data_tla_avg.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}).sub(sub).data = ft_timelockanalysis(cfg,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}).sub(sub).data);
%       end
%     end
%   end
% end

%% plot the conditions - simple

cfg_ft = [];
cfg_ft.xlim = [-0.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% %cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

% % cfg_plot.rois = {{'LAS'},{'LPS'}};
% cfg_plot.rois = {{'posterior'}};
% % cfg_plot.rois = {{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-8 8; -8 8];
% cfg_plot.legendlocs = {'SouthEast','NorthWest'};

% same as Scott et al. (2008)
cfg_plot.rois = {{'LPI2'},{'RPI2'}};
cfg_plot.ylims = [-2 6; -2 6];
cfg_plot.legendlocs = {'NorthEast','NorthEast'};

% cfg_ft.xlim = [-0.2 1.0];
% cfg_plot.rois = {{'E70'},{'E83'}};
% cfg_plot.ylims = [-10 10; -10 10];
% cfg_plot.legendlocs = {'NorthEast','NorthEast'};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

%cfg_plot.condByTypeByROI = repmat({{{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'cond1' 'cond2'}},size(cfg_plot.rois));

% sesNum = [1 2];
sesNum = [2];
% cfg_plot.condByROI = repmat(ana.eventValues{sesNum},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm'}},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm' 'stim2_basic_g' 'stim2_subord_g' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim1_basic_g' 'stim1_subord_g' 'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim1_basic_g_lo8' 'stim1_subord_g_lo8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_g' 'stim2_subord_g'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_g_hi8' 'stim2_subord_g_hi8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim2_basic_norm' 'stim2_subord_norm'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g' 'stim1_subord_g' 'stim2_basic_g' 'stim2_subord_g'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g_lo8' 'stim1_subord_g_lo8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_color' 'stim1_subord_color' 'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));


for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  %mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
  mm_ft_simpleplotER_multiSes(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
  %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
end

%% find the peak of N170

cfg = [];
cfg.conditions = cellflat(ana.eventValues{1});
cfg.data_str = 'ga_tla';
cfg.is_ga = true;
cfg.excludeBadSub = true;

sesNum = 1;

ana_str = mm_catSubStr_multiSes2(cfg,exper,sesNum);

cfg_ft = [];
% cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'posterior'})});
% cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS', 'RPS', 'LPI', 'RPI', 'LPS2', 'RPS2', 'LPI2', 'RPI2', 'PS', 'PI'})});
cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPI', 'RPI', 'LPI2', 'RPI2', 'PI'})});
cfg_ft.channel{end+1} = 'E81';

cfg_ft.latency = [0 0.5];


cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'RPI2'})});

ga_allCond = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',ana_str));

figure
ft_singleplotER(cfg_ft,ga_allCond);


%% find peaks

cfg = [];
cfg.latency = [0.155 0.211];

ga_n170 = ft_selectdata_new(cfg,ga_allCond);

[y,i] = sort(mean(ga_n170.avg,2));

ga_n170.label(i)'
y'

cfg = [];
cfg.latency = [0.23 0.33];

ga_n250 = ft_selectdata_new(cfg,ga_allCond);

[y,i] = sort(mean(ga_n250.avg,2));

ga_n250.label(i)'
y'

%% plot peaks

cfg_ft = [];
cfg_ft.showlabels = 'yes';
% cfg_ft.interactive = 'yes';
cfg_ft.interactive = 'no';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.ylim = [-2 5];
cfg_ft.layout = ft_prepare_layout([],ana);

figure
ft_multiplotER(cfg_ft,ga_allCond);

%% plot 90

cfg_ft = [];
cfg_ft.channel = {'E90'};

figure
ft_singleplotER(cfg_ft,ga_allCond);

%% peak to peak

pospeak = [0.1 0.15];

negpeak = [0.15 0.2];

cfg = [];
cfg.latency = pospeak;
ga_pos = ft_selectdata_new(cfg,ga_allCond);


cfg = [];
cfg.latency = negpeak;
ga_neg = ft_selectdata_new(cfg,ga_allCond);

pos_max = max(ga_pos.avg,[],2);
neg_min = min(ga_neg.avg,[],2);

posnegdiff = pos_max - abs(neg_min);

[y,i] = sort(posnegdiff);

ga_allCond.label(i)'
y'

%% BIG ANOVA

% test day (3) x image condition (5) x stim1/stim2 (2) x training basic/suborde (2) x hemi (2)

sessions = {'session_1', 'session_8', 'session_9'};
% sessions = {'session_1', 'session_8'};

imgConds = {'norm','g','g_hi8','g_lo8','color'};
training = {'basic', 'subord'};

stimNum = [1 2];

hemis = {'LPI2','RPI2'};

allBadSub = logical(sum(exper.badSub,2));

cfg = [];
% cfg.latency = [0.155 0.211]; % N170
% cfg.latency = [0.152 0.212]; % N170
cfg.latency = [0.23 0.33]; % N250

tbeg = nearest(data_tla.(exper.sesStr{1}).(sprintf('stim%d_%s_%s',stimNum(1),training{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(exper.sesStr{1}).(sprintf('stim%d_%s_%s',stimNum(1),training{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  %if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
      for im = 1:length(imgConds)
        
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%d_%s_%s',stimNum(sn),training{tr},imgConds{im});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis(h))});
              %dat = ft_selectdata_new(cfg,data_tla.(exper.sesStr{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(exper.sesStr{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(exper.sesStr{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  %end
end

varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord','Hemisphere'};
O = teg_repeated_measures_ANOVA(anovaData, [length(exper.sesStr) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord'};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training)], varnames);

%% subplots of each subject's ERPs

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.rois = {{'LAS'},{'LPS'}};
% cfg_plot.rois = {{'E70'},{'E83'}};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-0.2 1.0];
cfg_plot.ylim = [-10 10];
cfg_plot.parameter = 'avg';

% cfg_plot.rois = {{'E83'}};
% cfg_plot.xlim = [-0.2 1.0];
% cfg_plot.ylim = [-10 10];

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% % abbreviations for the condition types
% cfg_plot.typesByROI = {...
%   {'C2','C6'},...
%   {'C2','C6'}};

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotER(cfg_plot,ana,exper,data_tla);
end

% %% plot the conditions
% 
% cfg_ft = [];
% cfg_ft.xlim = [-0.2 1.5];
% cfg_ft.parameter = 'avg';
% 
% cfg_plot = [];
% 
% %cfg_plot.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -1 6; -1 6; -1 6];
% % vertical solid lines to plot
% cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
% cfg_plot.plotLegend = 0;
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
% cfg_plot.plotTitle = 0;
% 
% cfg_plot.is_ga = 1;
% cfg_plot.excludeBadSub = 1;
% 
% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% 
% % cfg_plot.condByTypeByROI = {...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   };
% % cfg_plot.typesByROI = repmat({{'C2','C6'}},size(cfg_plot.condByTypeByROI));
% 
% %cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));
% 
% for r = 1:length(cfg_plot.rois)
%   cfg_plot.roi = cfg_plot.rois{r};
%   %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
%   cfg_plot.conditions = cfg_plot.condByROI{r};
%   cfg_ft.ylim = cfg_plot.ylims(r,:);
%   cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
%   if cfg_plot.plotLegend
%     cfg_plot.legendloc = cfg_plot.legendlocs{r};
%   end
%   
%   mm_ft_plotERP(cfg_ft,cfg_plot,ana,exper,files,dirs,ga_tla);
% end

%% RSA - very basic
thisROI = 'PI';
elecInd = ismember(data_tla.(ana.eventValues{1}{1}).sub(1).ses(1).data.label,ana.elecGroups{ismember(ana.elecGroupsStr,thisROI)});

dataType = 'word_RgH_rc_spac';
% dataType = 'img_RgH_rc_spac';
% dataType = 'word_RgH_rc_mass';
% dataType = 'img_RgH_rc_mass';

% dataType = 'word_RgH_fo_spac';
% dataType = 'img_RgH_fo_spac';
% dataType = 'word_RgH_fo_mass';
% dataType = 'img_RgH_fo_mass';

phaseCountCol = 4;
stimNumCol = 6;
categNumCol = 7;

distanceMetric = 'euclidean';
% distanceMetric = 'seuclidean';
% distanceMetric = 'spearman';
% distanceMetric = 'cosine';
% distanceMetric = 'correlation';

if strcmp(distanceMetric,'euclidean')
  distanceScale = [0 100];
elseif strcmp(distanceMetric,'spearman') || strcmp(distanceMetric,'correlation') || strcmp(distanceMetric,'cosine')
  distanceScale = [0 2];
elseif strcmp(distanceMetric,'seuclidean')
  distanceScale = [0 20];
else
  distanceScale = [];
end

% timeS = [0.1 0.8];
timeS = [0.0 1.0];
timeInd = data_tla.(ana.eventValues{1}{1}).sub(1).ses(1).data.time >= timeS(1) & data_tla.(ana.eventValues{1}{1}).sub(1).ses(1).data.time <= timeS(2);

for i = 1:size(data_tla.(sprintf('%s_p1',dataType)).sub(1).ses(1).data.trial,1)
  %p1_trlInd = 1;
  p1_trlInd = i;
  p1_phaseCount = data_tla.(sprintf('%s_p1',dataType)).sub(1).ses(1).data.trialinfo(p1_trlInd,phaseCountCol);
  p1_stimNum = data_tla.(sprintf('%s_p1',dataType)).sub(1).ses(1).data.trialinfo(p1_trlInd,stimNumCol);
  p1_categNum = data_tla.(sprintf('%s_p1',dataType)).sub(1).ses(1).data.trialinfo(p1_trlInd,categNumCol);
  
  p2_trlInd = find(...
    data_tla.(sprintf('%s_p2',dataType)).sub(1).ses(1).data.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
    data_tla.(sprintf('%s_p2',dataType)).sub(1).ses(1).data.trialinfo(:,stimNumCol) == p1_stimNum & ...
    data_tla.(sprintf('%s_p2',dataType)).sub(1).ses(1).data.trialinfo(:,categNumCol) == p1_categNum);
  %p2_trlInd = 1;
  
  if ~isempty(p2_trlInd)
    % pdist2: rows (dim 1) are observations; columns (dim 2) are variables;
    % distances are measured between observations
    
    % rows = samples; cols = channels
    p1_data = squeeze(data_tla.(sprintf('%s_p1',dataType)).sub(1).ses(1).data.trial(p1_trlInd,elecInd,timeInd))';
    p2_data = squeeze(data_tla.(sprintf('%s_p2',dataType)).sub(1).ses(1).data.trial(p2_trlInd,elecInd,timeInd))';
    
    % rows = channels; cols = samples
    %p1_data = squeeze(data_tla.(sprintf('%s_p1',dataType)).sub(1).ses(1).data.trial(p1_trlInd,elecInd,timeInd));
    %p2_data = squeeze(data_tla.(sprintf('%s_p2',dataType)).sub(1).ses(1).data.trial(p2_trlInd,elecInd,timeInd));
    
    D = pdist2(p1_data,p2_data,distanceMetric);
    
    figure;
    if exist('distanceScale','var') && ~isempty(distanceScale)
      imagesc(linspace(timeS(1),timeS(2),size(p1_data,1)),linspace(timeS(1),timeS(2),size(p2_data,1)),D,distanceScale);
    else
      imagesc(linspace(timeS(1),timeS(2),size(p1_data,1)),linspace(timeS(1),timeS(2),size(p2_data,1)),D);
    end
    h = colorbar;
    set(get(h,'YLabel'),'string','Dissimilarity');
    title(sprintf('%s (%d vs %d): phaseCount=%d stimNum=%d categNum=%d\n',strrep(dataType,'_','-'),p1_trlInd,p2_trlInd,p1_phaseCount,p1_stimNum,p1_categNum));
    xlabel('P1');
    ylabel('P2');
  else
    fprintf('%s: No p2 found for p1 phaseCount=%d stimNum=%d categNum=%d\n',dataType,p1_phaseCount,p1_stimNum,p1_categNum);
  end
end

%% make some GA plots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [-0.2 1.0]; % time
cfg_ft.xlim = [-0.2 1.5]; % time
%cfg_ft.xlim = [-0.2 2.0]; % time

cfg_ft.parameter = 'avg';

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

cfg_plot.ftFxn = 'ft_singleplotER';
cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5; -2 5];
cfg_plot.rois = {{'FC'}};
cfg_plot.ylims = [-7.5 2];
cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.plotLegend = 0;
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};

% cfg_plot.xlabel = 'Time (s)';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

% cfg_plot.ftFxn = 'ft_topoplotER';
% cfg_plot.ylims = [-6 6];
% %cfg_plot.ylims = 'maxmin';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% %cfg_ft.comment = 'no';
% % cfg_plot.rois = {'all'};
% % cfg_ft.xlim = [0 1.2]; % time
% % cfg_plot.subplot = 1;
% cfg_plot.rois = {{'LAS'}};
% cfg_ft.xlim = [0.3 0.5]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [0.5 0.8]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [1.0 1.5]; % time

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% %cfg_plot.rois = {{'FS'},{'LAS','RAS'},{'LPS','RPS'}};
% %cfg_plot.rois = {{'FS'},{'PS'}};
% %cfg_plot.rois = {'E71'};
% cfg_plot.rois = {'all'};
% cfg_plot.ylims = [-5 5]; % voltage in multiplot
% %cfg_plot.ylims = repmat('maxmin',size(cfg_plot.rois,2),1); % voltage in multiplot

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
% cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{{'FSC','FSI','N'}}},size(cfg_plot.rois));
cfg_plot.rename_condByROI = repmat({{{'FSC','FSI','CR'}}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));
%cfg_plot.rename_condByROI = repmat({{{'SC','SI','CR'}}},size(cfg_plot.rois));

% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% cfg_plot.condByTypeByROI = repmat({{{'HSC2','HSI2','CR2'},{'HSC6','HSI6','CR6'}}},size(cfg_plot.rois));
% % cfg_plot.condByTypeByROI = {...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   };
% cfg_plot.typesByROI = repmat({{'C2','C6'}},size(cfg_plot.condByTypeByROI));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.types = cfg_plot.typesByROI{r};
  cfg_plot.rename_conditions = cfg_plot.rename_condByROI{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  
  if strcmp(cfg_plot.ftFxn,'ft_singleplotER')
    cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
    if cfg_plot.plotLegend
      cfg_plot.legendloc = cfg_plot.legendlocs{r};
    end
  end
  
  mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_ft = [];
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
cfg_ft.colormap = 'jet';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'no';

%cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';
%cfg_plot.conditions = {{'HSC2','CR2'},{'HSI2','CR2'},{'HSC2','HSI2'},{'HSC6','CR6'},{'HSI6','CR6'},{'HSC6','HSI6'}}; % {'H2','CR2'}, {'H6','CR6'},
%cfg_plot.conditions = {{'RH','RCR'},{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}};
% cfg_plot.conditions = {{'SC','CR'},{'SC','SI'},{'SI','CR'}};
%cfg_plot.conditions = {{'RHSC','RHSI'}};
%cfg_plot.conditions = {{'FSC','RSSI'}};
cfg_plot.conditions = {{'Face','House'}};


cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-5.5 5.5]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
% cfg_ft.comment = 'no';
cfg_ft.comment = 'xlim';
cfg_ft.commentpos = 'middletop';

cfg_plot.roi = {'E73'};
%cfg_plot.roi = {'LAS'};
% cfg_ft.xlim = [0.01 0.8]; % time

% cfg_plot.roi = {'LPS','RPS'};
% cfg_ft.xlim = [0.5 0.8]; % time

%cfg_plot.roi = {'RAS'};
%cfg_ft.xlim = [1.1 1.9]; % time
% cfg_ft.xlim = [0 1.5]; % time
% cfg_plot.roi = {'all'};

% cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time
cfg_ft.xlim = (0:0.05:1.0); % time

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.xlim = [-0.2 1.5]; % time
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.ylim = [-1 1]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_ft.xlim = [-0.2 1.5]; % time
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.ylim = [-1 1]; % volt

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
% and the times that correspond to each set of ROIs

%cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'},{'LPS','RPS'},{'LPS','RPS'},{'LAS','RAS'},{'LPS','RPS'},{'LPS','RPS'},{'LPS','RPS'}};
%cfg_ana.latencies = [0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];

% cfg_ana.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_ana.latencies = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];

cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_ana.rois = {{'FC'}};
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% % LF O/N
% cfg_ana.rois = {{'RAS'},{'RAS'},{'RAI'},{'RAI'}};
% cfg_ana.latencies = [1.1 1.4; 1.4 1.9; 1.1 1.4; 1.4 1.9];

% % LPN
% cfg_ana.rois = {{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_ana.latencies = [1.2 1.8; 1.2 1.8; 1.2 1.8];


% cfg_ana.conditions = {'all'};
%cfg_ana.conditions = {{'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},{'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'}};
%cfg_ana.conditions = {{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}}; % {'RH','RCR'},
%cfg_ana.conditions = {{'SC','CR'},{'SI','CR'},{'SC','SI'}}; % {'RH','CR'},
cfg_ana.conditions = {{'Face','House'}};
cfg_ana.conditions = {{'RgH','CR'}};

%cfg_ana.conditions = {{'all'}};
%cfg_ana.conditions = {{'all_within_types'}};
%cfg_ana.conditions = {{'all_across_types'}};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.parameter = 'avg';
cfg_ft.correctm = 'fdr';

% line plot parameters
cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
%cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 1 4; 1 4; 1 4; -4 -1; 1 4; 1 4; 1 4];
%cfg_plot.ylims = [-4 -1; 1 4; 1 4];
% cfg_plot.ylims = [-5 -2; -5 -2; -5 -2; 1.5 4.5; 1.5 4.5;];
cfg_plot.ylims = [-5 -2; 1.5 4.5];
%cfg_plot.plot_order = {'CR2','H2','HSC2','HSI2','CR6','H6','HSC6','HSI6'};
%cfg_plot.plot_order = {'RHSC','RHSI','RCR'};
% cfg_plot.plot_order = {'SC','SI','CR'};
%cfg_plot.rename_conditions = {'SC','SI','CR'};
% cfg_plot.xlabel = 'Condition';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  % use data_tla_avg because FieldTrip doesn't deal with dimord properly
  % when single-trials exist
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla_avg);
end

%% output some values

cfg = [];

cfg.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg.latencies = [0.3 0.5; 0.5 0.8];
% cfg.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg.latencies = [0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];
cfg.condByROI = repmat({{'SC','SI','CR'}},size(cfg.rois));

% cfg.rois = {{'FC'}};
% cfg.latencies = [0.3 0.5];
% cfg.condByROI = repmat({{'FSC','FSI','N'}},size(cfg.rois));

cfg.parameter = 'avg';

cfg.excludeBadSub = false;

%cfg.direction = 'columns';
cfg.direction = 'rows';

mm_printDataToText(cfg,exper,ana,dirs,data_tla);

%% 3-way ANOVA: Hemisphere x Block Type x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% IV2: abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% IV3: outermost cell holds one cell for each ROI; each ROI cell holds one
% cell for each event type; each event type cell holds strings for its
% conditions
cfg_ana.condByTypeByROI = {...
  %{{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
  {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}},...
  {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

% For each ROI, what's common among the conditions in each type
cfg_ana.condCommonByROI = {...
  %{'CR','H','HSC','HSI'},...
  {'CR','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Block Type','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByTypeByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov33ER(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Hemisphere x Condition - separate colors

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% % IV2: define the conditions tested for each set of ROIs
% cfg_ana.condByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% 
% % For each ROI, what's common among the conditions in each type
% cfg_ana.condCommonByROI = {...
%   {'CR','H','HSC','HSI'},...
%   {'CR','HSC','HSI'}};
% 
% % abbreviations for the condition types
% cfg_ana.typesByROI = {...
%   {'C2','C6'},...
%   {'C2','C6'}};

% IV2: define the conditions tested for each set of ROIs
%cfg_ana.condByROI = {{'RCR','RH'},{'RCR','RHSC','RHSI'}};
cfg_ana.condByROI = {{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
% for the one-way ANOVAs
cfg_ana.condCommonByROI = {...
  %{'CR','H'},...
  {'CR','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  %cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% IV2: define the conditions tested for each set of ROIs
%cfg_ana.condByROI = {{'RH','RCR'},{'RCR','RHSC','RHSI'}};
%cfg_ana.condByROI = {{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
cfg_ana.condByROI = {exper.eventValues,exper.eventValues};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% cfg_ana.condCommonByROI = {...
%   {'CR','HSC','HSI'},...
%   {'CR','HSC','HSI'}};
cfg_ana.condCommonByROI = {...
  exper.eventValues,...
  exper.eventValues};

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% 1-way ANOVA: Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 0;

% define which regions to average across for the test
cfg_ana.rois = {{'FC'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5];

cfg_ana.condByROI = repmat({{'FSC', 'FSI', 'N'}},size(cfg_ana.rois));
% Define the IVs (type: event, roi, latency)
cfg_ana.IV1.name = 'FSC/FSI/CR';
cfg_ana.IV1.cond = {'FSC', 'FSI', 'N'};
cfg_ana.IV1.type = 'event';

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  
  mm_ft_rmaov1ER_spec(cfg_ana,exper,ana,data_tla);
end

%% cluster statistics

cfg_ft = [];
%cfg_ft.avgovertime = 'no';
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'no';

cfg_ft.parameter = 'avg';

cfg_ana = [];
cfg_ana.roi = 'all';
%cfg_ana.latencies = [0 1.0; 1.0 2.0];
%cfg_ana.latencies = [0.2 0.6];
cfg_ana.latencies = [0.3 0.5];

cfg_ft.numrandomization = 500;
% cfg_ft.clusteralpha = 0.05;
cfg_ft.clusteralpha = 0.1;
cfg_ft.alpha = 0.05;

% extra directory info
cfg_ana.dirStr = '';
if strcmp(cfg_ft.avgovertime,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end
if strcmp(cfg_ft.avgoverchan,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgC'];
end

% cfg_ana.conditions = {...
%   {'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},...
%   {'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'},...
%   {'CR2','CR6'},{'H2','H6'},{'HSC2','HSC6'},{'HSI2','HSI6'}};
cfg_ana.conditions = {'all_within_types'};
%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  
  stat_clus = mm_ft_clusterstatER(cfg_ft,cfg_ana,exper,ana,dirs,data_tla);
end

%% plot the cluster statistics

%files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = 0.1;

cfg_plot = [];
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.dirStr = cfg_ana.dirStr;

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  
  mm_ft_clusterplotER(cfg_ft,cfg_plot,ana,files,dirs);
end

%% let me know that it's done
emailme = 1;
if emailme
  subject = sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    };
  send_gmail(subject,mail_message);
end

%% correlations

cfg_ana = [];

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.dpTypesByROI = {...
  {'Item','Source'},...
  {'Item','Source'}};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds two cells, one for each
% d' type; each d' cell contains strings for its conditions
cfg_ana.condByROI = {...
  {{{'CR2','H2'},{'HSC2','HSI2'}}, {{'CR6','H6'},{'HSC6','HSI6'}}}...
  {{{'CR2','H2'},{'HSC2','HSI2'}}, {{'CR6','H6'},{'HSC6','HSI6'}}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% C2 d' values
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_item =  abs([]);
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_source =  abs([]);
% C6 d' values
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_item =  abs([]);
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_source =  abs([]);

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data_tla);
end
