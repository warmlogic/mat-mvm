%% load the analysis details

expName = 'EBUG';

subDir = '';
% subDir = 'new_art500';
if ~isempty(subDir)
  warning('Loading data from subDir ''%s''...',subDir);
end
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
    'EBUG010';
    'EBUG012';
%     'EBUG016';
%     'EBUG017';
%     'EBUG018';
%     'EBUG019';
    'EBUG020';
%     'EBUG022';
%     'EBUG025';
%     'EBUG027';
%     'EBUG029';
%     'EBUG030';
%     'EBUG032';
%     'EBUG045';
%     'EBUG047';
%     'EBUG051';
%     'EBUG052';
%     'EBUG054';
%     'EBUG055';
%     'EBUG061';
  };

% only one cell, with all session names
sesNames = {'session_1'};
% sesNames = {'session_1','session_8','session_9'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot,[],subDir);


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
  {{'stim1_basic_untrained', 'stim1_subord_untrained','stim2_basic_untrained', 'stim2_subord_untrained' ...
  'stim1_basic_trained', 'stim1_subord_trained','stim2_basic_trained', 'stim2_subord_trained'} ...
  }},1,length(exper.sessions));

% ana.eventValuesSplit = repmat({...
%   {{'stim1_basic_color', 'stim1_basic_g', 'stim1_basic_g_hi8', 'stim1_basic_g_lo8', 'stim1_basic_norm', ...
%   'stim1_subord_color', 'stim1_subord_g', 'stim1_subord_g_hi8', 'stim1_subord_g_lo8', 'stim1_subord_norm'} ...
%   }},1,length(exper.sessions));

ana.trl_expr = cell(1,length(exper.sessions));

% for session 1, could collapse trained/untrained because this is first
% exposure
ses = 1;
ana.trl_expr{ses} = {...
  {...
  sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
  }...
  };

% ses = 2;
% ana.trl_expr{ses} = {...
%   {...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   }...
%   };
% 
% ses = 3;
% ana.trl_expr{ses} = {...
%   {...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 0',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 1 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 0 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   sprintf('eventNumber == %d & stimNum == 2 & isSubord == 1 & response ~= 0 & rt < 2900 & trained == 1',find(ismember(exper.eventValues{ses},'match_stim'))) ...
%   }...
%   };

%% load in the subject data

% [data_tla,exper] = mm_ft_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');

keeptrials = false;
[data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo',true);

% %% get rid of the bad channels
% 
% cfg = [];
% cfg.printRoi = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% [data_tla] = mm_rmBadChan(cfg,exper,ana,data_tla);

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

%% save

saveDir = dirs.saveDirProc;
save(fullfile(saveDir,'ebug_match_data_tla.mat'),'data_tla','exper','ana','dirs','files','-v7.3');
% clear data_tla

%% load

subDir = '';
dataDir = fullfile('EBUG','EEG','Sessions','ftpp',subDir);
% Possible locations of the data files (dataroot)
serverDir = fullfile(filesep,'Volumes','curranlab','Data');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dreamDir = fullfile(filesep,'data','projects','curranlab');
localDir = fullfile(getenv('HOME'),'data');

% pick the right dataroot
if exist('serverDir','var') && exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
else
  error('Data directory not found.');
end

loadDir = fullfile(dataroot,dataDir,'ft_data/data_art_nsClassic_ftAuto/tla');

load(fullfile(loadDir,'ebug_match_data_tla.mat'));

replaceDataType = {};

[dirs] = mm_checkDirs(dirs,replaceDataType,subDir);

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {{},{},{}};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,14,[],'vert');
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
for typ = 1:length(ana.eventValues{ses})
  for ev = 1:length(ana.eventValues{ses}{typ})
    figure
    ft_multiplotER(cfg_ft,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{ev}).sub(sub).data);
    title(strrep(ana.eventValues{ses}{typ}{ev},'_','-'));
  end
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

%% plot the conditions - simple

cfg_ft = [];
% cfg_ft.xlim = [0 1.0];
cfg_ft.xlim = [-0.05 1.0];
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
% cfg_plot.rois = {{'LPI2'},{'RPI2'}};
cfg_plot.rois = {{'LPI2','RPI2'}};
cfg_plot.ylims = [-0.5 6; -0.5 6];
% cfg_plot.legendlocs = {'NorthEast','NorthEast'};
cfg_plot.legendlocs = {'SouthEast','SouthEast'};

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
% sesNum = [1 3];
% sesNum = [2 3];
sesNum = [1];
% sesNum = [2];
% sesNum = [3];
% cfg_plot.condByROI = repmat(ana.eventValues{sesNum},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm' 'stim2_basic_g' 'stim2_subord_g' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8' 'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm' 'stim2_basic_color' 'stim2_subord_color' 'stim2_basic_g' 'stim2_subord_g'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_g' 'stim2_subord_g' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim1_basic_g' 'stim1_subord_g' 'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim1_basic_g_lo8' 'stim1_subord_g_lo8' 'stim1_basic_color' 'stim1_subord_color'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim1_basic_color' 'stim1_subord_color' 'stim1_basic_g' 'stim1_subord_g'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g' 'stim1_subord_g' 'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim1_basic_g_lo8' 'stim1_subord_g_lo8'}},size(cfg_plot.rois));

cfg_plot.condByROI = repmat({{'stim1_basic_trained' 'stim1_subord_trained'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_color' 'stim1_subord_color'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g' 'stim1_subord_g'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g_hi8' 'stim1_subord_g_hi8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g_lo8' 'stim1_subord_g_lo8'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_g' 'stim2_subord_g'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_g_hi8' 'stim2_subord_g_hi8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_basic_g' 'stim2_basic_g_hi8' 'stim2_basic_g_lo8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_subord_norm' 'stim2_subord_g' 'stim2_subord_g_hi8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim2_basic_norm' 'stim2_subord_norm'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g' 'stim1_subord_g' 'stim2_basic_g' 'stim2_subord_g'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_g_lo8' 'stim1_subord_g_lo8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim1_basic_color' 'stim1_subord_color' 'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));

% cfg_ft.graphcolor = 'rmbckgrykbrmbckgrykb';
% cfg_ft.linestyle = {'-','-','-','-','-','-','--','--','--','--', '-.','-.','-.','-.','-.','-.','-.','-.','-.','-.'};

% cfg_ft.graphcolor = 'rmbckgrmbckg';
cfg_ft.graphcolor = 'mrcbgkmrcbgkmrcbgk';
cfg_ft.linestyle = {'-','-','-','-','-','-','--','--','--','--','--','--','-.','-.','-.','-.','-.','-.'};

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  

  %mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
  mm_ft_simpleplotER_multiSes(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
  %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
  
  stimStr = 'Stim1';
%   imgCond = 'Incongruent';
%   imgCond = 'Gray-AllSF';
%   imgCond = 'HiSF';
%   imgCond = 'LoSF';
  
  legendLoc = 'SouthEast';
  legend({sprintf('%s Basic Pretest %s',stimStr) sprintf('%s Subord Pretest %s',stimStr) sprintf('%s Basic 1-Day Posttest %s',stimStr) sprintf('%s Subord 1-Day Posttest %s',stimStr)},'Location',legendLoc);
  %legend({sprintf('%s Basic %s Pretest',stimStr,imgCond) sprintf('%s Subord %s Pretest',stimStr,imgCond) sprintf('%s Basic %s 1-Day Posttest',stimStr,imgCond) sprintf('%s Subord %s 1-Day Posttest',stimStr,imgCond)},'Location',legendLoc);
  %legend off
  
  saveFigs = true;
  
  figFileName = sprintf('ERP_%s_%s_basic_subord_prepost',imgCond,stimStr);
  
  publishfig(gcf,1,14);
  if saveFigs
    print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
  end
end


%% make average voltage plots with this function

% mm_ft_avgplotER_multiSes


%% lineplot of the conditions

% don't use this, move ahead to next two sections

% cfg_ft = [];
% % cfg_ft.xlim = [-0.2 1.0];
% % cfg_ft.parameter = 'avg';
% 
% cfg_plot = [];
% 
% cfg_ft.latency = [0.156 0.208]; % N170 (use this)
% % cfg_ft.latency = [0.234 0.334]; % N250
% 
% cfg_plot.is_ga = 0;
% cfg_plot.excludeBadSub = 1;
% 
% % %cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% % cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% % cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5];
% % cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};
% 
% % % cfg_plot.rois = {{'LAS'},{'LPS'}};
% % cfg_plot.rois = {{'posterior'}};
% % % cfg_plot.rois = {{'LPS'},{'RPS'}};
% % cfg_plot.ylims = [-8 8; -8 8];
% % cfg_plot.legendlocs = {'SouthEast','NorthWest'};
% 
% % same as Scott et al. (2008)
% cfg_plot.rois = {{'LPI2'},{'RPI2'}};
% % cfg_plot.rois = {{'LPI2','RPI2'}};
% cfg_plot.ylims = [-2 6; -2 6];
% cfg_plot.legendlocs = {'NorthEast','NorthEast'};
% 
% % cfg_ft.xlim = [-0.2 1.0];
% % cfg_plot.rois = {{'E70'},{'E83'}};
% % cfg_plot.ylims = [-10 10; -10 10];
% % cfg_plot.legendlocs = {'NorthEast','NorthEast'};
% 
% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% 
% % cfg_plot.condByTypeByROI = {...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% 
% %cfg_plot.condByTypeByROI = repmat({{{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}},size(cfg_plot.rois));
% 
% % cfg_plot.condByROI = repmat({{'cond1' 'cond2'}},size(cfg_plot.rois));
% 
% % sesNum = [1 2];
% % sesNum = [1 3];
% % sesNum = [2 3];
% sesNum = [1];
% % sesNum = [2];
% % sesNum = [3];
% % cfg_plot.condByROI = repmat(ana.eventValues{sesNum},size(cfg_plot.rois));
% 
% % cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm'}},size(cfg_plot.rois));
% 
% % cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm' 'stim2_basic_g' 'stim2_subord_g' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8' 'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_subord_norm' 'stim2_basic_g' 'stim2_subord_g' 'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim2_basic_g' 'stim2_subord_g' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% 
% % cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim1_basic_g' 'stim1_subord_g' 'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim1_basic_g_lo8' 'stim1_subord_g_lo8' 'stim1_basic_color' 'stim1_subord_color'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim1_basic_g' 'stim1_subord_g' 'stim1_basic_color' 'stim1_subord_color'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim1_basic_g' 'stim1_subord_g' 'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim1_basic_g_lo8' 'stim1_subord_g_lo8'}},size(cfg_plot.rois));
% 
% % cfg_plot.condByROI = repmat({{'stim2_basic_g' 'stim2_subord_g'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim2_basic_g_hi8' 'stim2_subord_g_hi8'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));
% 
% 
% % cfg_plot.condByROI = repmat({{'stim2_basic_norm' 'stim2_basic_g' 'stim2_basic_g_hi8' 'stim2_basic_g_lo8'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim2_subord_norm' 'stim2_subord_g' 'stim2_subord_g_hi8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% 
% 
% % cfg_plot.condByROI = repmat({{'stim1_basic_norm' 'stim1_subord_norm' 'stim2_basic_norm' 'stim2_subord_norm'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim1_basic_g' 'stim1_subord_g' 'stim2_basic_g' 'stim2_subord_g'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim1_basic_g_hi8' 'stim1_subord_g_hi8' 'stim2_basic_g_hi8' 'stim2_subord_g_hi8'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim1_basic_g_lo8' 'stim1_subord_g_lo8' 'stim2_basic_g_lo8' 'stim2_subord_g_lo8'}},size(cfg_plot.rois));
% % cfg_plot.condByROI = repmat({{'stim1_basic_color' 'stim1_subord_color' 'stim2_basic_color' 'stim2_subord_color'}},size(cfg_plot.rois));
% 
% % cfg_ft.graphcolor = 'rmbckgrykbrmbckgrykb';
% % cfg_ft.linestyle = {'-','-','-','-','-','-','--','--','--','--', '-.','-.','-.','-.','-.','-.','-.','-.','-.','-.'};
% 
% % cfg_ft.graphcolor = 'rmbckgrmbckg';
% cfg_ft.graphcolor = 'mrcbgkmrcbgkmrcbgk';
% cfg_ft.linestyle = {'-','-','-','-','-','-','--','--','--','--','--','--','-.','-.','-.','-.','-.','-.'};
% 
% collapseSessions = false;
% if collapseSessions
%   exper.badSub = logical(sum(exper.badSub,2));
%   exper.badSub = repmat(exper.badSub,1,length(exper.sessions));
% end
% 
% for r = 1:length(cfg_plot.rois)
%   cfg_plot.roi = cfg_plot.rois{r};
%   cfg_plot.legendloc = cfg_plot.legendlocs{r};
%   cfg_ft.ylim = cfg_plot.ylims(r,:);
%   %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
%   cfg_plot.conditions = cfg_plot.condByROI{r};
%   
% 
%   %mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
% %   mm_ft_lineplotER_multiSes(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
%   mm_ft_lineplotER_multiSes(cfg_ft,cfg_plot,ana,exper,sesNum,data_tla);
%   %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
% end

%% LINES: break pre/post tests into different graphs - training x img cond

% pretest only will collapse across training

% dataMeasure = 'dp';
% dataLabel = 'd''';
% ylimits = [0 4];

% dataMeasure = 'rt_hit';
% dataLabel = 'Response Time: Hits';
% ylimits = [0 5000];

% dataMeasure = 'hr';
% dataLabel = 'Hit Rate';
% ylimits = [0 1];

% dataMeasure = 'c';
% dataLabel = 'Response bias (criterion; c)';
% ylimits = [-0.6 0.6];
% % positive/conservative bias indicates a tendency to say 'new', whereas
% % negative/liberal bias indicates a tendency to say 'old'

% dataMeasure = 'Br';
% dataLabel = 'Response bias index (Br)';
% ylimits = [0 1];

% dataMeasure = 'Pr';
% dataLabel = 'Discrimination index (Pr)';
% ylimits = [0 1];

sessions = {'session_1', 'session_8', 'session_9'};
% sessions_str = {'Pretest', '1-Day Posttest','Delay'};
sessions_str = {'Pretest', '1-Day Posttest','1-Week Posttest'};

% if collapsePhases
phases = {'match'};
% else
%   phases = {'match_1'};
% end
% training = {'trained','untrained'};

% naming = {'basic','subord'};
% naming_str = {'Basic','Subord'};
naming = {'subord'};
naming_str = {'Subord'};

% training = {'TT','UU','TU','UT'};
% trainStr = {'Trained','Untrained','TU','UT'};

% training = {'TT','UU'};
% trainingStr = {'Trained', 'Untrained'};

stimNum = {'1','2'};
% stimNum = {'1'};
% stimNum = {'2'};
% stimNum_str = {'Stim1','Stim2'};
stimNum_str = {'1','2'};

% % separate hemi
% hemis = {'LPI2','RPI2'};
% hemis_str = {'Left','Right'};

% collapse hemi
hemis = {{'LPI2','RPI2'}};
hemis_str = {'LR'};

imgConds = {'norm','color','g'};
imgConds_str = {'Congruent','Incongruent','Gray'};
groupname = 'Color';

% imgConds = {'g','g_hi8','g_lo8'};
% % imgConds_str = {'Gray','Hi8','Lo8'};
% imgConds_str = {'AllSF','HiSF','LoSF'};
% groupname = 'SpatialFreq';

allBadSub = logical(sum(exper.badSub,2));

cfg = [];
% cfg.latency = [0.155 0.211]; % N170 (Scott)
% cfg.latency = [0.152 0.212]; % N170 (initial guess)
cfg.latency = [0.156 0.208]; % N170 (use this)
component_str = 'N170';

% cfg.latency = [0.234 0.334]; % N250
% component_str = 'N250';

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},naming{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},naming{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

%anovaData = [];
%theseData = [];

data_sub = nan(length(sessions),length(imgConds),length(naming),length(stimNum),length(hemis),length(subjects));

for ses = 1:length(sessions)
  for im = 1:length(imgConds)
    
    for n = 1:length(naming)
      for sn = 1:length(stimNum)
        
        condition = sprintf('stim%s_%s_%s',stimNum{sn},naming{n},imgConds{im});
        
        for h = 1:length(hemis)
          for sub = 1:length(exper.subjects)
            if ~allBadSub(sub)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis{h})});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              data_sub(ses,im,n,sn,h,sub) = mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2);
              
            end
          end
        end
      end
    end
    %anovaData = cat(1,anovaData,theseData);
  end
end

% for s = 1:length(sessions)
%   for p = 1:length(phases)
%     for i = 1:length(imgConds)
%       
%       %collapseData = [];
%       
%       for t = 1:length(training)
%         for n = 1:length(naming)
%           data(t,i,s,p,:) = results.(sessions{s}).(phases{p}).(training{t}).(imgConds{i}).(naming{n}).(dataMeasure);
%         end
%       end
%       
%     end
%   end
% end

% % stats
% [h, p, ci, stats] = ttest(squeeze(data.(dataMeasure).(naming{1})(:,1,1,1,1)),squeeze(data.(dataMeasure).(naming{1})(:,2,1,1,1)));

%% Make line plots - separate session, hemi, stimNum

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;
if ~verLessThan('matlab', '8.4')
  cfg_plot.removeErrBarEnds = false;
end

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 7];
end

cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};

if strcmp(groupname,'Color')
  cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
elseif strcmp(groupname,'SpatialFreq')
  cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    for sn = 1:length(stimNum)
      
      data_mean = [];
      data_sem = [];
      %data_mean = nan(length(naming),length(imgConds));
      %data_sem = nan(length(naming),length(imgConds));
      
      thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      figFileName = sprintf('%s_%s_%s_%s_s%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      
      cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
      cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));
      
      poCount = 0;
      for im = 1:length(imgConds)
        for n = 1:length(naming)
          poCount = poCount + 1;
          
          condition = sprintf('%s_%s',imgConds{im},naming{n});
          cfg_plot.plot_order{poCount} = condition;
          %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
          cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));
          
          data_mean.(condition) = nanmean(data_sub(s,im,n,sn,h,:),6);
          data_sem.(condition) = nanstd(data_sub(s,im,n,sn,h,:),0,6) ./ sqrt(sum(~allBadSub));
        end
      end
      
      hpm = nan(1,length(cfg_plot.plot_order));
      
      figure
      % plot the lines
      eval(sprintf('plot([%s],cfg_plot.linespec{1},''LineWidth'',cfg_plot.linewidth);',sprintf(repmat('data_mean.%s ',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:})));
      hold on
      for c = 1:length(cfg_plot.plot_order)
        % errorbars
        heb = errorbar(c,mean(data_mean.(cfg_plot.plot_order{c}),1),data_sem.(cfg_plot.plot_order{c}),cfg_plot.linespec{c},'LineWidth',cfg_plot.errwidth);
        % remove errorbar ends
        if cfg_plot.removeErrBarEnds
          chil = get(heb,'Children');
          xdata = get(chil(2),'XData');
          ydata = get(chil(2),'YData');
          xdata(cfg_plot.errBarEndMarkerInd) = NaN;
          ydata(cfg_plot.errBarEndMarkerInd) = NaN;
          set(chil(2),'XData',xdata);
          set(chil(2),'YData',ydata);
          set(heb,'Children',chil);
        end
        % plot the markers
        hpm(c) = plot(c,mean(data_mean.(cfg_plot.plot_order{c}),1),cfg_plot.linespec{c},'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor{c});
        
        text(c,  cfg_plot.ylim(1)-.2, cfg_plot.rename_conditions{c}, 'Rotation', 15, 'FontSize', 18, 'HorizontalAlignment', 'center');
      end
      if cfg_plot.plotLegend
        legend(hpm(1:length(cfg_plot.legendtext)),cfg_plot.legendtext);
      end
      
      cfg_plot.xlim = [0.5 (length(cfg_plot.rename_conditions) + .5)];
      
      plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal
      
      hold off
      
      set(gcf,'Name',sprintf('%s',thisTitle))
      
      title(thisTitle);
      
      % make it look good
      axis([cfg_plot.xlim(1) cfg_plot.xlim(2) cfg_plot.ylim(1) cfg_plot.ylim(2)])
      xlabel(cfg_plot.xlabel);
      ylabel(cfg_plot.ylabel);
      set(gca,'XTick',(1:length(cfg_plot.rename_conditions)))
      set(gca,'XTickLabel',[])
      %set(gca,'XTickLabel',strrep(cfg_plot.rename_conditions,'_',''))
      set(gca,'YTick',(cfg_plot.ylim(1):.5:cfg_plot.ylim(2)))
      axis square
      publishfig(gcf,0);
      
      if saveFigs
        print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
      end
      
    end
  end
end

%% bar plots of subordinate - collapse stim order

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;
if ~verLessThan('matlab', '8.4')
  cfg_plot.removeErrBarEnds = false;
end

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
  cfg_plot.ylim = [0 3];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 6];
end

% cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};
% 
% if strcmp(groupname,'Color')
%   cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
% elseif strcmp(groupname,'SpatialFreq')
%   cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
% end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

ylimits = cfg_plot.ylim;

cfg_plot.axisxy = false;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    %for sn = 1:length(stimNum)
      
      %thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      %figFileName = sprintf('%s_%s_%s_%s_s%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      thisTitle = sprintf('%s: %s: %s, %s hemi',groupname,component_str,sessions_str{s},hemis_str{h});
      figFileName = sprintf('bar_%s_%s_%s_%s',groupname,component_str,sessions_str{s},hemis_str{h});
      
      cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
      cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));
      
      data_mean = [];
      data_sem = [];
      %data_mean = nan(length(naming),length(imgConds));
      %data_sem = nan(length(naming),length(imgConds));
      
      poCount = 0;
      for im = 1:length(imgConds)
        for n = 1:length(naming)
          poCount = poCount + 1;
          
          condition = sprintf('%s_%s',imgConds{im},naming{n});
          cfg_plot.plot_order{poCount} = condition;
          %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
          cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));
          
          %data_mean.(condition) = nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6);
          %data_sem.(condition) = nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub));
          
          data_mean = cat(2,data_mean,nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6));
          data_sem = cat(2,data_sem,nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub)));
        end
      end
      
      bw_legend = imgConds_str;
      
      bw_title = thisTitle;
      %bw_title = sprintf('%s: %s, %s',groupname,sesStr{s},naming{n});
      %bw_title = sprintf('%s%s: %s: %s',upper(naming{n}(1)),naming{n}(2:end),strrep(phases{p},'_','\_'),strrep(imgConds{i},'_','\_'));
      %bw_groupnames = {'Pretest', 'Posttest', 'Delay'};
      if strcmp(sessions{s},'session_1')
        bw_groupnames = [];
        axis_x = 1.5;
      else
        %bw_groupnames = trainingStr;
        bw_groupnames = [];
        
        %axis_x = length(training) + 1.5;
        %axis_x = length(training) + 0.5;
        axis_x = 1.5;
      end
      %bw_xlabel = 'Test day';
      bw_xlabel = cfg_plot.xlabel;
      bw_ylabel = cfg_plot.ylabel;
      %if exist('linspecer','file')
      %  bw_colormap = 'linspecer';
      %else
      %  bw_colormap = 'gray';
      %end
      
      bw_colormap = flipud(jet);
      bw_data = data_mean;
      bw_errors = data_sem;
      bw_width = 0.75;
      hbp = barweb(bw_data,bw_errors,bw_width,bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend,[],'plot');
      set(hbp.legend,'Location','NorthEast');
      axis([0.5 axis_x ylimits(1) ylimits(2)]);
      
      %if length(sessions) == 1 && strcmp(sessions{s},'pretest')
      %  xlabel('Collapsed');
      %end
      publishfig(gcf,0);
      
      if cfg_plot.axisxy
        %axis xy;
        set(gca,'YDir','reverse');
        set(hbp.legend,'Location','SouthWest');
        figFileName = sprintf('%s_flip',figFileName);
      end
      
      figFileName = strrep(strrep(figFileName,' ','_'),'-','_');
      
      if saveFigs
        print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
        %print(gcf,figFormat,figRes,fullfile(figsDir,sprintf('prepost_trainUn_%s_%s_%s_%s_%s',phases{p},dataMeasure,sessions{s},naming{n},groupname)));
      end
      
  end
end

%% bar plots of subordinate - separate stim order

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;
if ~verLessThan('matlab', '8.4')
  cfg_plot.removeErrBarEnds = false;
end

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
  cfg_plot.ylim = [0 3];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 6];
end

% cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};
% 
% if strcmp(groupname,'Color')
%   cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
% elseif strcmp(groupname,'SpatialFreq')
%   cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
% end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

ylimits = cfg_plot.ylim;

cfg_plot.axisxy = true;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    for sn = 1:length(stimNum)
      
      thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      figFileName = sprintf('bar_%s_%s_%s_%s_s%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      %thisTitle = sprintf('%s: %s: %s, %s hemi',groupname,component_str,sessions_str{s},hemis_str{h});
      %figFileName = sprintf('bar_%s_%s_%s_%s',groupname,component_str,sessions_str{s},hemis_str{h});
      
      cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
      cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));
      
      data_mean = [];
      data_sem = [];
      %data_mean = nan(length(naming),length(imgConds));
      %data_sem = nan(length(naming),length(imgConds));
      
      poCount = 0;
      for im = 1:length(imgConds)
        for n = 1:length(naming)
          poCount = poCount + 1;
          
          condition = sprintf('%s_%s',imgConds{im},naming{n});
          cfg_plot.plot_order{poCount} = condition;
          %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
          cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));
          
          %data_mean.(condition) = nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6);
          %data_sem.(condition) = nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub));
          
          data_mean = cat(2,data_mean,nanmean(data_sub(s,im,n,sn,h,:),6));
          data_sem = cat(2,data_sem,nanstd(data_sub(s,im,n,sn,h,:),0,6) ./ sqrt(sum(~allBadSub)));
        end
      end
      
      bw_legend = imgConds_str;
      
      bw_title = thisTitle;
      %bw_title = sprintf('%s: %s, %s',groupname,sesStr{s},naming{n});
      %bw_title = sprintf('%s%s: %s: %s',upper(naming{n}(1)),naming{n}(2:end),strrep(phases{p},'_','\_'),strrep(imgConds{i},'_','\_'));
      %bw_groupnames = {'Pretest', 'Posttest', 'Delay'};
      if strcmp(sessions{s},'session_1')
        bw_groupnames = [];
        axis_x = 1.5;
      else
        %bw_groupnames = trainingStr;
        bw_groupnames = [];
        
        %axis_x = length(training) + 1.5;
        %axis_x = length(training) + 0.5;
        axis_x = 1.5;
      end
      %bw_xlabel = 'Test day';
      bw_xlabel = cfg_plot.xlabel;
      bw_ylabel = cfg_plot.ylabel;
      %if exist('linspecer','file')
      %  bw_colormap = 'linspecer';
      %else
      %  bw_colormap = 'gray';
      %end
      
      bw_colormap = flipud(jet);
      bw_data = data_mean;
      bw_errors = data_sem;
      bw_width = 0.75;
      hbp = barweb(bw_data,bw_errors,bw_width,bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend,[],'plot');
      set(hbp.legend,'Location','NorthEast');
      axis([0.5 axis_x ylimits(1) ylimits(2)]);
      
      %if length(sessions) == 1 && strcmp(sessions{s},'pretest')
      %  xlabel('Collapsed');
      %end
      publishfig(gcf,0);
      
      if cfg_plot.axisxy
        %axis xy;
        set(gca,'YDir','reverse');
        set(hbp.legend,'Location','SouthWest');
        figFileName = sprintf('%s_flip',figFileName);
      end
      
      figFileName = strrep(strrep(figFileName,' ','_'),'-','_');
      
      if saveFigs
        print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
        %print(gcf,figFormat,figRes,fullfile(figsDir,sprintf('prepost_trainUn_%s_%s_%s_%s_%s',phases{p},dataMeasure,sessions{s},naming{n},groupname)));
      end
      
    end
  end
end

%% Make line plots - separate session (collapse hemi, stimNum)

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;
if ~verLessThan('matlab', '8.4')
  cfg_plot.removeErrBarEnds = false;
end

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 7];
end

cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};

if strcmp(groupname,'Color')
  cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
elseif strcmp(groupname,'SpatialFreq')
  cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    %for sn = 1:length(stimNum)
      
      data_mean = [];
      data_sem = [];
      %data_mean = nan(length(naming),length(imgConds));
      %data_sem = nan(length(naming),length(imgConds));
      
      %thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      %figFileName = sprintf('%s_%s_%s_%s_s%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      thisTitle = sprintf('%s: %s: %s, %s hemi',groupname,component_str,sessions_str{s},hemis_str{h});
      figFileName = sprintf('%s_%s_%s_%s',groupname,component_str,sessions_str{s},hemis_str{h});
      
      cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
      cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));
      
      poCount = 0;
      for im = 1:length(imgConds)
        for n = 1:length(naming)
          poCount = poCount + 1;
          
          condition = sprintf('%s_%s',imgConds{im},naming{n});
          cfg_plot.plot_order{poCount} = condition;
          %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
          cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));
          
          data_mean.(condition) = nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6);
          data_sem.(condition) = nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub));
        end
      end
      
      hpm = nan(1,length(cfg_plot.plot_order));
      
      figure
      % plot the lines
      eval(sprintf('plot([%s],cfg_plot.linespec{1},''LineWidth'',cfg_plot.linewidth);',sprintf(repmat('data_mean.%s ',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:})));
      hold on
      for c = 1:length(cfg_plot.plot_order)
        % errorbars
        heb = errorbar(c,mean(data_mean.(cfg_plot.plot_order{c}),1),data_sem.(cfg_plot.plot_order{c}),cfg_plot.linespec{c},'LineWidth',cfg_plot.errwidth);
        % remove errorbar ends
        if cfg_plot.removeErrBarEnds
          chil = get(heb,'Children');
          xdata = get(chil(2),'XData');
          ydata = get(chil(2),'YData');
          xdata(cfg_plot.errBarEndMarkerInd) = NaN;
          ydata(cfg_plot.errBarEndMarkerInd) = NaN;
          set(chil(2),'XData',xdata);
          set(chil(2),'YData',ydata);
          set(heb,'Children',chil);
        end
        % plot the markers
        hpm(c) = plot(c,mean(data_mean.(cfg_plot.plot_order{c}),1),cfg_plot.linespec{c},'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor{c});
        
        text(c,  cfg_plot.ylim(1)-.2, cfg_plot.rename_conditions{c}, 'Rotation', 15, 'FontSize', 18, 'HorizontalAlignment', 'center');
      end
      if cfg_plot.plotLegend
        legend(hpm(1:length(cfg_plot.legendtext)),cfg_plot.legendtext);
      end
      
      cfg_plot.xlim = [0.5 (length(cfg_plot.rename_conditions) + .5)];
      
      plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal
      
      hold off
      
      set(gcf,'Name',sprintf('%s',thisTitle))
      
      title(thisTitle);
      
      % make it look good
      axis([cfg_plot.xlim(1) cfg_plot.xlim(2) cfg_plot.ylim(1) cfg_plot.ylim(2)])
      xlabel(cfg_plot.xlabel);
      ylabel(cfg_plot.ylabel);
      set(gca,'XTick',(1:length(cfg_plot.rename_conditions)))
      set(gca,'XTickLabel',[])
      %set(gca,'XTickLabel',strrep(cfg_plot.rename_conditions,'_',''))
      set(gca,'YTick',(cfg_plot.ylim(1):.5:cfg_plot.ylim(2)))
      axis square
      publishfig(gcf,0);
      
      if saveFigs
        print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
      end
      
    %end
  end
end

%% BIG ANOVA

% test day (3) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

cfg = [];

cfg.latency = [0.156 0.208]; % N170 (use this)
component_str = 'N170';

cfg.latency = [0.234 0.334]; % N250
component_str = 'N250';

sessions = {'session_1', 'session_8', 'session_9'};
lev_sessions = {'pre','post','delay'};

% sessions = {'session_1', 'session_8'};
% lev_sessions = {'pre','post'};

% sessions = {'session_1', 'session_9'};
% lev_sessions = {'pre','delay'};

% sessions = {'session_8', 'session_9'};
% lev_sessions = {'post','delay'};

% groupname = 'SpatialFreq';
% imgConds = {'g','g_hi8','g_lo8'};
% lev_imgConds = {'gray','g_hi8','g_lo8'};

training = {'basic', 'subord'};

stimNum = {'1' '2'};
lev_stimNum = {'s1' 's2'};
% stimNum = {'2'};
% lev_stimNum = {'s2'};

hemis = {'LPI2','RPI2'};
lev_hemis = {'left','right'};

allBadSub = logical(sum(exper.badSub,2));

fprintf('ANOVA: %s %s\n',component_str);

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
        
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%s_%s_%s',stimNum{sn},training{tr});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis(h))});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'TestDay','StimNum','Basic/Subord','Hemisphere'};

levelnames = {lev_sessions lev_stimNum training lev_hemis};

O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord'};
% levelnames = {sessions imgConds stimNum training};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training)], varnames);

%% BIG ANOVA - stim1 only, collapse hemi

% test day (3) x image condition (5) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

cfg = [];

cfg.latency = [0.156 0.208]; % N170 (use this)
component_str = 'N170';

cfg.latency = [0.234 0.334]; % N250
component_str = 'N250';

sessions = {'session_1', 'session_8', 'session_9'};
lev_sessions = {'pre','post','delay'};

% sessions = {'session_1', 'session_8'};
% lev_sessions = {'pre','post'};

% sessions = {'session_1', 'session_9'};
% lev_sessions = {'pre','delay'};

% sessions = {'session_8', 'session_9'};
% lev_sessions = {'post','delay'};

groupname = 'Color';
imgConds = {'norm','color','g'};
lev_imgConds = {'cong','incon','gray'};

groupname = 'SpatialFreq';
imgConds = {'g','g_hi8','g_lo8'};
lev_imgConds = {'gray','g_hi8','g_lo8'};

training = {'basic', 'subord'};

% stimNum = {'1' '2'};
% lev_stimNum = {'s1' 's2'};
stimNum = {'1'};
lev_stimNum = {'s1'};

% hemis = {'LPI2','RPI2'};
% lev_hemis = {'left','right'};

hemis = {{'LPI2','RPI2'}};
lev_hemis = {'LR'};

allBadSub = logical(sum(exper.badSub,2));

fprintf('ANOVA: %s %s\n',groupname,component_str);

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
      for im = 1:length(imgConds)
        
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%s_%s_%s',stimNum{sn},training{tr},imgConds{im});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis{h})});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

% varnames = {'TestDay', 'ImgCond','StimNum','Basic/Subord','Hemisphere'};
% 
% levelnames = {lev_sessions lev_imgConds lev_stimNum training lev_hemis};
% 
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% % O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);

varnames = {'testDay', 'ImgCond','Basic/Subord'};
levelnames = {lev_sessions lev_imgConds training};
O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training)], varnames,[],[],[],[],[],[],levelnames);

%% normal condition only ANOVA

% test day (3) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

sessions = {'session_1', 'session_8', 'session_9'};
% sessions = {'session_1', 'session_8'};
% sessions = {'session_1', 'session_9'};
% sessions = {'session_8', 'session_9'};


training = {'basic', 'subord'};

stimNum = {'1' '2'};
% stimNum = {'2'};

hemis = {'LPI2','RPI2'};

allBadSub = logical(sum(exper.badSub,2));

cfg = [];
% cfg.latency = [0.155 0.211]; % N170 (Scott)
% cfg.latency = [0.152 0.212]; % N170 (initial guess)
cfg.latency = [0.156 0.208]; % N170 (use this)

% cfg.latency = [0.234 0.334]; % N250

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
      
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%s_%s_%s',stimNum{sn},training{tr});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis(h))});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
        end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'testDay','StimNum','Basic/Subord','Hemisphere'};

lev_sessions = {'pre','post','delay'};
% lev_sessions = {'pre','post'};
% lev_sessions = {'pre','delay'};
% lev_sessions = {'post','delay'};

lev_stimNum = {'s1' 's2'};
% lev_stimNum = {'s2'};

lev_hemis = {'left','right'};

levelnames = {lev_sessions lev_stimNum training lev_hemis};

O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord'};
% levelnames = {sessions imgConds stimNum training};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training)], varnames);

%% Pretest ANOVA

% test day (3) x image condition (5) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

sessions = {'session_1'};

imgConds = {'norm','g','color'};
% imgConds = {'g','g_hi8','g_lo8'};

training = {'basic', 'subord'};

stimNum = {'1' '2'};
% stimNum = {'2'};

hemis = {'LPI2','RPI2'};

allBadSub = logical(sum(exper.badSub,2));

cfg = [];
% cfg.latency = [0.155 0.211]; % N170 (Scott)
% cfg.latency = [0.152 0.212]; % N170 (initial guess)
cfg.latency = [0.156 0.208]; % N170 (use this)

cfg.latency = [0.234 0.334]; % N250

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
      for im = 1:length(imgConds)
        
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%s_%s_%s',stimNum{sn},training{tr},imgConds{im});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis(h))});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'ImgCond','StimNum','Basic/Subord','Hemisphere'};

%lev_sessions = {'pre'};

lev_imgConds = {'cong','gray','incon'};
% lev_imgConds = {'gray','g_hi8','g_lo8'};

lev_stimNum = {'s1' 's2'};
% lev_stimNum = {'s2'};

lev_hemis = {'left','right'};

levelnames = {lev_imgConds lev_stimNum training lev_hemis};

O = teg_repeated_measures_ANOVA(anovaData, [length(imgConds) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord'};
% levelnames = {sessions imgConds stimNum training};
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

%% make some GA plots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [-0.2 1.0]; % time
% cfg_ft.xlim = [-0.2 1.5]; % time
%cfg_ft.xlim = [-0.2 2.0]; % time
cfg_ft.xlim = [-0.5 0.6]; % time

cfg_ft.parameter = 'avg';

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_plot.rois = {{'LPI2'},{'RPI2'}};
cfg_plot.rois = {{'LPI2','RPI2'}};
cfg_plot.ylims = [-0.5 6; -0.5 6];
% cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.x_bounds = [];
cfg_plot.plotLegend = 0;
cfg_plot.legendlocs = {'SouthEast','SouthEast'};

cfg_plot.xlabel = 'Time (s)';
cfg_plot.ylabel = 'Voltage (\muV)';
% cfg_plot.xlabel = '';
% cfg_plot.ylabel = '';

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

% cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_ana.rois = {{'FC'}};

cfg_ana.rois = {{'LPI2','RPI2'},{'LPI2','RPI2'},{'LPI2','RPI2'}};

cfg_ana.latencies = [0.156 0.208; 0.156 0.208; 0.156 0.208];

% cfg.latency = [0.156 0.208]; % N170 (use this)
% 
% cfg.latency = [0.234 0.334]; % N250


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
% cfg_ana.conditions = {{'Face','House'}};
% cfg_ana.conditions = {{'RgH','CR'}};

cfg_ana.conditions = {{'stim2_basic_norm' 'stim2_subord_norm'} {'stim2_basic_g' 'stim2_subord_g'} {'stim2_basic_color' 'stim2_subord_color'}};


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
cfg_plot.line_plots = 1;
%cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 1 4; 1 4; 1 4; -4 -1; 1 4; 1 4; 1 4];
%cfg_plot.ylims = [-4 -1; 1 4; 1 4];
% cfg_plot.ylims = [-5 -2; -5 -2; -5 -2; 1.5 4.5; 1.5 4.5;];
cfg_plot.ylims = [-2 6; -2 6; -2 6];
%cfg_plot.plot_order = {'CR2','H2','HSC2','HSI2','CR6','H6','HSC6','HSI6'};
%cfg_plot.plot_order = {'RHSC','RHSI','RCR'};
% cfg_plot.plot_order = {'SC','SI','CR'};
%cfg_plot.rename_conditions = {'SC','SI','CR'};
% cfg_plot.xlabel = 'Condition';
cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
% cfg_plot.ylabel = '';

sesNum = 1;

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  % use data_tla_avg because FieldTrip doesn't deal with dimord properly
  % when single-trials exist
  %mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla_avg,sesNum);
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla,sesNum);
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
