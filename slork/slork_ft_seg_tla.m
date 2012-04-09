% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SLORK';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
%exper.eventValues = sort({'SCR','SHSC','SHSI','TCR','THSC','THSI'});
exper.eventValues = sort({'SF','SN','SRO','SRS','TF','TN','TRO','TRS'});

% combine some events into higher-level categories
%exper.eventValuesExtra.toCombine = {{'SHSC','SHSI'},{'THSC','THSI'}};
%exper.eventValuesExtra.newValue = {{'SH'},{'TH'}};
exper.eventValuesExtra.toCombine = {{'SF','TF'},{'SN','TN'},{'SRO','TRO'},{'SRS','TRS'}};
exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 0;

% make sure eventValuesExtra exists because create_ft_struct uses it
if ~isfield(exper,'eventValuesExtra')
  exper.eventValuesExtra = {};
  exper.eventValuesExtra.onlyKeepExtras = 0;
end

% equate the number of trials across event values?
exper.equateTrials = 0;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 2.0];

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
exper.eegFileExt = 'egis';
%exper.eegFileExt = 'raw';

% name of the folder to save the FT data in
if isempty(exper.eventValuesExtra)
  evStr = sprintf(repmat('%s_',1,length(exper.eventValues)),exper.eventValues{:});
elseif isempty(exper.eventValuesExtra) && ~isfield(exper.eventValuesExtra,'onlyKeepExtras')
  evStr = sprintf(repmat('%s_',1,length(exper.eventValues)),exper.eventValues{:});
elseif ~isempty(exper.eventValuesExtra) && isfield(exper.eventValuesExtra,'onlyKeepExtras') && exper.eventValuesExtra.onlyKeepExtras == 0
  evStr = cat(2,exper.eventValues,cat(2,exper.eventValuesExtra.newValue{:}));
  evStr = sprintf(repmat('%s_',1,length(evStr)),evStr{:});
elseif ~isempty(exper.eventValuesExtra) && isfield(exper.eventValuesExtra,'onlyKeepExtras') && exper.eventValuesExtra.onlyKeepExtras == 1
  evStr = cat(2,cat(2,exper.eventValuesExtra.newValue{:}));
  evStr = sprintf(repmat('%s_',1,length(evStr)),evStr{:});
end
dirs.saveDirName = fullfile('ft_data',sprintf('tla_%seq%d',evStr,exper.equateTrials));

exper.subjects = {
  'SLORK002';
  'SLORK003';
  'SLORK004';
  'SLORK005';
  'SLORK006';
  'SLORK007';
  'SLORK008';
  'SLORK009';
  'SLORK010';
  'SLORK011';
  'SLORK012';
  'SLORK013';
  'SLORK014';
  'SLORK015';
  'SLORK016';
  'SLORK017';
  'SLORK018';
  'SLORK019';
  'SLORK020';
  'SLORK022';
  'SLORK023';
  'SLORK024';
  'SLORK025';
  'SLORK026';
  'SLORK027';
  'SLORK029';
  'SLORK030';
  'SLORK031';
  'SLORK032';
  'SLORK033'; % master's analyses included up to 33
  'SLORK034';
  'SLORK035';
  'SLORK036';
  'SLORK038';
  'SLORK039';
  'SLORK040';
  'SLORK041';
  'SLORK042';
  'SLORK043';
  };
  % 'SLORK001'; % exper was not set up correctly
  % 'SLORK021'; % didn't record
  % 'SLORK028'; % didn't record
  % 37 does not exist

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_0'};

%% set up file and directory handling parameters

dirs.homeDir = getenv('HOME');

% directory where the data to read is located
%dirs.subDir = 'RK';
dirs.subDir = '';
dirs.dataDir = fullfile('eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
%dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));
% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile('/Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile('/Volumes','RAID','curranlab','Data');
dirs.localDir = fullfile(getenv('HOME'),'data');

% pick the right dirs.dataroot
if exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
elseif exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
elseif exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
else
  error('Data directory not found.');
end

% directory to save the data
dirs.saveDirProc = fullfile(dirs.dataroot,dirs.saveDirName);
if ~exist(dirs.saveDirProc,'dir')
  mkdir(dirs.saveDirProc)
end

% Assumes we have chan locs file in ~/Documents/MATLAB/mat_mvm/eeg/; use
% "short" because the Fiduciary points in the FT locs screws up some plots
files.elecfile = fullfile(dirs.homeDir,'Documents/MATLAB/mat_mvm/eeg/GSN_HydroCel_129_short.sfp');
%files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

% figure printing options - see mm_ft_setSaveDirs for other options
files.saveFigs = 0;
files.figPrintFormat = 'png';
%files.figPrintFormat = 'epsc2';
if strcmp(files.figFileExt,'eps')
  files.figPrintFormat = '-depsc2';
elseif strcmp(files.figFileExt,'pdf')
  files.figPrintFormat = '-dpdf';
elseif strcmp(files.figFileExt,'png')
  files.figPrintFormat = '-dpng';
elseif strcmp(files.figFileExt,'jpg')
  files.figPrintFormat = '-djpg90';
end

% directory to save figures
dirs.saveDirFigs = fullfile(dirs.saveDirProc,'figs');
if ~exist(dirs.saveDirFigs,'dir')
  mkdir(dirs.saveDirFigs)
end

% %% add NS's artifact information to the event structure
% nsEvFilters.eventValues = exper.eventValues;
% % SCR
% nsEvFilters.SCR.type = 'LURE_PRES';
% nsEvFilters.SCR.testType = 'side';
% nsEvFilters.SCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % SHSC
% nsEvFilters.SHSC.type = 'TARG_PRES';
% nsEvFilters.SHSC.testType = 'side';
% nsEvFilters.SHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % SHSI
% nsEvFilters.SHSI.type = 'TARG_PRES';
% nsEvFilters.SHSI.testType = 'side';
% nsEvFilters.SHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% % TCR
% nsEvFilters.TCR.type = 'LURE_PRES';
% nsEvFilters.TCR.testType = 'task';
% nsEvFilters.TCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % THSC
% nsEvFilters.THSC.type = 'TARG_PRES';
% nsEvFilters.THSC.testType = 'task';
% nsEvFilters.THSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % THSI
% nsEvFilters.THSI.type = 'TARG_PRES';
% nsEvFilters.THSI.testType = 'task';
% nsEvFilters.THSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% 
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sessions)
%     overwriteArtFields = 1;
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,overwriteArtFields);
%   end
% end

%% Convert the data to FieldTrip structs - excludes NS artifact trials
ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_timelockanalysis';

% any preprocessing?
cfg_pp = [];
% single precision to save space
cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.keeptrials = 'no';
[data_tla,exper] = create_ft_struct(ana,cfg_pp,cfg_proc,exper,dirs,files);

% save the structs for loading in later
if strcmp(cfg_proc.keeptrials,'no')
  saveFile = fullfile(dirs.saveDirProc,sprintf('data_tla_avg_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000));
elseif strcmp(cfg_proc.keeptrials,'yes')
  saveFile = fullfile(dirs.saveDirProc,sprintf('data_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000));
end
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'data_tla');
  fprintf('Done.\n');
else
  error('Not saving! %s already exists.\n',saveFile);
end

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

% list the values separated by types: Side, Question (Task)
%ana.eventValues = {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}};
%ana.eventValues = {{'SF','SN','SRO','SRS'},{'TF','TN','TRO','TRS'}};
ana.eventValues = {{'F','N','RO','RS'}};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end

%% save the analysis details

saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000));
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  fprintf('Done.\n');
else
  error('Not saving! %s already exists.\n',saveFile);
end

%% if already saved and not yet loaded, load the ft_timelockanalysis files

if ~exist('cfg_proc','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,'analysisDetails*.mat'));
  if length(savedFiles) == 1
    load(fullfile(dirs.saveDirProc,savedFiles.name));
  elseif length(savedFiles) > 1
    error('Multiple analysisDetails*.mat files found in %s!',dirs.saveDirProc)
  elseif isempty(savedFiles)
    error('analysisDetails*.mat not found in %s!',dirs.saveDirProc)
  end
end

if ~exist('data_tla','var')
  if strcmp(cfg_proc.keeptrials,'no')
    savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('data_tla_avg_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000)));
  elseif strcmp(cfg_proc.keeptrials,'yes')
    savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('data_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000)));
  end
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDirProc,savedFiles(sf).name));
    fprintf('Done.\n');
  end
  % get all the exper.eventValues and exper.eventValuesExtra together; make sure the extra event values aren't in the list
  if ~isempty(exper.eventValuesExtra)
    if exper.eventValuesExtra.onlyKeepExtras
      exper.eventValues = cat(2,exper.eventValuesExtra.newValue{:});
    else
      for nVal = 1:length(exper.eventValuesExtra.newValue)
        if ~ismember(exper.eventValuesExtra.newValue{nVal},exper.eventValues)
          exper.eventValues = cat(2,exper.eventValues,exper.eventValuesExtra.newValue{nVal});
        else
          fprintf('%s is already in the event value list!\n',exper.eventValuesExtra.newValue{nVal}{1});
        end
      end
    end
    exper.eventValues = sort(exper.eventValues);
  end
end

%% Test plots to make sure data look ok

% cfg_ft = [];
% cfg_ft.showlabels = 'yes';
% cfg_ft.interactive = 'yes';
% cfg_ft.showoutline = 'yes';
% cfg_ft.fontsize = 9;
% cfg_ft.layout = ft_prepare_layout([],ana);
% figure
% ft_multiplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% 
% cfg_ft = [];
% cfg_ft.channel = {'E20'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% %figure
% %ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data,data_tla.(exper.eventValues{2}).sub(1).ses(1).data,data_tla.(exper.eventValues{3}).sub(1).ses(1).data);
% %legend(strrep(exper.eventValues,'_',''),'Location','SouthEast');
% figure
% cfg_ft.graphcolor = 'b';
% ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% legend(strrep(exper.eventValues{1},'_',''),'Location','SouthEast');
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
% cfg_ft.layout = ft_prepare_layout([],data_tla);
% figure
% ft_topoplotER(cfg_ft,data_topo);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Remove the dof field. I'm not sure what degrees of freedom refers to, but the tutorial says to.
% for evVal = 1:length(exper.eventValues)
%   for sub = 1:length(exper.subjects)
%     for ses = 1:length(exper.sessions)
%       if isfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof')
%         data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data = rmfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof');
%       end
%     end
%   end
% end

%% decide who to kick out based on trial counts

numEv = struct;
numEv.thresh = 15;

% Subjects with bad behavior
exper.badBehSub = {};

% exclude subjects with low event counts - use events in ana.eventValues
[exper,numEv] = mm_threshSubs(exper,ana,numEv,data_tla);

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = exper.eventValues;
cfg_ana.data_str = 'data_tla';
cfg_ana.excludeBadSub = 1;
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sessions)
  for evVal = 1:length(exper.eventValues)
    %tic
    fprintf('Running ft_timelockgrandaverage on %s...',exper.eventValues{evVal});
    ga_tla.(exper.eventValues{evVal})(ses) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(exper.eventValues{evVal}){ses}));
    fprintf('Done.\n');
    %toc
  end
end

% cfg_ft = [];
% cfg_ft.keepindividual = 'no';
% for ses = 1:length(exper.sessions)
%   typ = 1;
%   for evVal = 1:length(ana.eventValues{typ})
%     %tic
%     fprintf('Running ft_timelockgrandaverage on %s...',ana.eventValues{typ}{evVal});
%     ga_tla.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
%     fprintf('Done.\n');
%     %toc
%   end
% end

%% save grand average file

saveFile = fullfile(dirs.saveDirProc,sprintf('ga_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000));
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'ga_tla');
  fprintf('Done.\n');
else
  error('Not saving! %s already exists.\n',saveFile);
end

%% (re)save the analysis details

saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000));
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','numEv');
else
  fprintf('Appending to %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','numEv','-append');
end
fprintf('Done.\n');

%% let me know that it's done
emailme = 1;
if emailme
  subject = sprintf('Done with%s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with%s %s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    sprintf('%s',saveFile),...
    };
  send_gmail(subject,mail_message);
end

%% if already saved and not yet loaded, load the grand average file

if ~exist('ga_tla','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('ga_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000)));
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDirProc,savedFiles(sf).name));
    fprintf('Done.\n');
  end
end

%% plot the conditions - simple

cfg_ft = [];
cfg_ft.xlim = [-.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];
cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.ylims = [-5 2; -2 5];
cfg_plot.legendlocs = {'SouthEast','NorthWest'};

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
% cfg_plot.condByTypeByROI = {...
%   {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}},...
%   {{'SCR','SHSC','SHSI'},{'TCR','THSC','THSI'}}};
cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,ga_tla);
end

%% subplots of each subject's ERPs

cfg_plot = [];
cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
%cfg_plot.roi = {'E124'};
%cfg_plot.roi = {'RAS'};
%cfg_plot.roi = {'LPS','RPS'};
%cfg_plot.roi = {'LPS'};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-.2 1.0];
cfg_plot.ylim = [-10 10];
cfg_plot.parameter = 'avg';

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_plot.condByTypeByROI = {...
  {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}},...
  {{'SCR','SHSC','SHSI'},{'TCR','THSC','THSI'}}};
% abbreviations for the condition types
cfg_plot.typesByROI = {...
  {'Side','Question'},...
  {'Side','Question'}};

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.types = cfg_plot.typesByROI{r};
  
  mm_ft_subjplotER(cfg_plot,ana,exper,numEv,data_tla);
end

%% plot the conditions

cfg_ft = [];
cfg_ft.xlim = [-0.2 2.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

%cfg_plot.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
cfg_plot.rois = {{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-5 2; -5 2; -5 2; -1 6; -1 6; -1 6];
% vertical solid lines to plot
cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.plotLegend = 1;
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
cfg_plot.plotTitle = 1;

cfg_plot.is_ga = 1;
% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_plot.condByTypeByROI = {...
  {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}},...
  {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}},...
  {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}},...
  {{'SCR','SHSC','SHSI'},{'TCR','THSC','THSI'}}...
  {{'SCR','SHSC','SHSI'},{'TCR','THSC','THSI'}}...
  {{'SCR','SHSC','SHSI'},{'TCR','THSC','THSI'}}...
  };

cfg_plot.typesByROI = repmat({{'Side','Question'}},size(cfg_plot.condByTypeByROI));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
  cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.types = cfg_plot.typesByROI{r};
  
  mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 1;

cfg_ft = [];
cfg_ft.xlim = [-0.2 1]; % time
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'yes';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'yes';

cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.condMethod = 'pairwise';

cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-2 2]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.comment = 'no';
cfg_ft.xlim = [0.5 0.8]; % time
%cfg_ft.xlim = [0 1.0]; % time
%cfg_ft.xlim = (0:0.05:1.0); % time
%cfg_plot.subplot = 1;
cfg_plot.roi = {'PS'};

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.ylim = [-1 1]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.ylim = [-1 1]; % volt

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);

%% Plot average time contrast topoplots - old function
% 
% cfg_plot = [];
% cfg_plot.plotTitle = 1;
% cfg_plot.plotColorbar = 1;
% 
% % define which regions to highlight in the plot
% cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% % define the times that correspond to each set of rois
% cfg_plot.xlims = [0.3 0.5; 0.5 0.8];
% cfg_plot.conditions = {{'all_within_types'}};
% 
% cfg_ft = [];
% cfg_ft.interactive = 'no';
% %cfg_ft.colorbar = 'yes';
% cfg_ft.highlight = 'on';
% cfg_ft.highlightsize = 10;
% %cfg_ft.comment = 'xlim';
% %cfg_ft.commentpos = 'title';
% cfg_ft.comment = 'no';
% %cfg_ft.colormap = 'hot'; % dark to light; better for b&w printers
% cfg_ft.zlim = [-1 1];
% 
% for r = 1:length(cfg_plot.rois)
%   cfg_plot.roi = cfg_plot.rois{r};
%   cfg_ft.xlim = cfg_plot.xlims(r,:);
%   
%   mm_ft_contrasttopoER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);
% end

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.conditions = {{'SCR','SH'},{'SCR','SHSC'},{'SCR','SHSI'},{'SHSC','SHSI'},{'TCR','TH'},{'TCR','THSC'},{'TCR','THSI'},{'THSC','THSI'}};

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
cfg_plot.ylims = [-4 -1; 1 4];
cfg_plot.plot_order = {'SCR','SH','SHSC','SHSI','TCR','TH','THSC','THSI'};

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla);
end

%% 3-way ANOVA: Hemisphere x Block Type x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_ana.condByTypeByROI = {...
  {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}},...
  {{'SCR','SHSC','SHSI'},{'TCR','THSC','THSI'}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'Side','Question'},...
  {'Side','Question'}};

% For each ROI, what's common among the conditions in each type
cfg_ana.condCommonByROI = {...
  {'CR','H','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Block Type','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByTypeByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  
  mm_ft_rmaov33ER(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% IV2: define the conditions tested for each set of ROIs
cfg_ana.condByROI = {...
  {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}},...
  {{'SCR','SH','SHSC','SHSI'},{'TCR','TH','THSC','THSI'}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'Side','Question'},...
  {'Side','Question'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% cluster statistics

cfg_ft = [];
cfg_ft.avgovertime = 'no';

cfg_ft.parameter = 'avg';

cfg_ana = [];
cfg_ana.roi = 'all';
cfg_ana.latencies = [0 1.0; 1.0 2.0];

%cfg_ana.conditions = {{'SHSC','THSC'}};
%cfg_ana.conditions = {'all_within_types'};
cfg_ana.conditions = {{'SCR','SH'},{'SCR','SHSC'},{'SCR','SHSI'},{'SHSC','SHSI'},{'TCR','TH'},{'TCR','THSC'},{'TCR','THSI'},{'THSC','THSI'}};

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  
  stat_clus = mm_ft_clusterstatER(cfg_ft,cfg_ana,exper,ana,dirs,data_tla);
end

%% plot the cluster statistics

files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = .1;

cfg_plot = [];
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.conditions = cfg_ana.conditions;

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
  {{{'SCR','SH'},{'SHSC','SHSI'}}, {{'TCR','TH'},{'THSC','THSI'}}}...
  {{{'SCR','SH'},{'SHSC','SHSI'}}, {{'TCR','TH'},{'THSC','THSI'}}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'Side','Question'},...
  {'Side','Question'}};

% Side d' values
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_item = abs([1.6145, 2.0823, 2.6991, 1.7986, 1.9783, 1.0971, 2.0902, 2.0895, 0.6578, 0.6442, 1.678, 1.7853, 1.6032, 2.4474, 2.3005, 1.7255, 1.566, 1.2529, 1.0774, 3.0318, 1.9626, 1.5029, 0.8762, 2.7775, 1.553, 1.3498, 0.8494, 0.844, 0.7109, 1.2076, 0.3676, 3.1351, 0.6876, 2.1477, 1.5546, 1.729, 1.729, 2.5606, 1.6701]);
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_source = abs([1.5528, 1.7537, 1.9537, 1.7535, 1.687, 0.0788, 1.2857, 2.0013, 0.4771, 0.6571, 1.1276, 1.8787, 1.1546, 3.0483, 1.4379, 0.7955, 1.7102, 0.8845, 1.1945, 1.7385, 1.7149, 2.3272, 0.0968, 3.2844, 1.2291, 0.9366, -0.6108, 0.6999, 0.2037, 0.622, -0.038, 2.5101, 0.0176, 0.3024, 1.6043, 0.6393, 0.3983, 2.6648, 0.4885]);
% Question d' values
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_item = abs([1.8692, 2.0083, 3.0126, 1.1594, 1.6155, 1.1885, 1.6652, 1.7839, 0.9898, 1.0621, 1.5352, 2.1099, 1.4867, 2.6831, 2.3211, 1.7398, 1.7557, 1.5432, 1.4436, 2.4062, 2.2162, 1.7094, 1.39, 2.8908, 1.301, 1.0106, 1.0331, 1.2649, 1.6038, 1.5393, 0.9931, 2.5151, 0.6482, 2.2035, 1.5761, 1.5546, 1.8078, 2.2108, 1.7255]);
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_source = abs([0.7439, 1.0587, 1.5439, 0.7418, 0.8039, -0.3197, 1.0668, 1.317, 0.0736, 0.5171, 1.3582, 1.2274, 0.738, 1.6176, 2.1666, 0.5236, 1.4651, 1.1363, 0.3559, 0.8078, 1.0654, 1.8929, -0.2507, 1.3646, 1.2062, 0.3971, -0.1053, 1.2216, 0.511, 0.8607, -0.1276, 0.7271, -0.2012, 1.187, 1.0989, 1.6598, 1.3599, 1.401, 0.9302]);

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data_tla);
end
