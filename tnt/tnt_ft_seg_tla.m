% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for a info on nonparametric statistics


% contrasts:

% ANOVA:

% T/NT x 1/2 x R/F
% T/NT/B x 1/2

% ttest
% T vs NT (whole session)
% T vs Nt (2)
% NTF vs NTR (2)

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'TNT';

exper.sampleRate = 250;

% pre- and post-stimulus times used to segment NS files, in seconds (pre is
% negative)
exper.prepost = [-1.0 1.7];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
exper.eegFileExt = 'egis';
%exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
exper.eventValues = sort({'B1of2','B2of2','NT1of2For','NT1of2Rec','NT2of2For','NT2of2Rec','T1of2For','T1of2Rec','T2of2For','T2of2Rec'});
%exper.eventValues = sort({'NT1of2For','NT1of2Rec','NT2of2For','NT2of2Rec','T1of2For','T1of2Rec','T2of2For','T2of2Rec'});
%exper.eventValues = sort({'NT1of2Rec','T1of2Rec'});

% combine some events into higher-level categories

% % B and T/NT x 1/2 (collapse across R/F)
% exper.eventValuesExtra.toCombine = {...
%   {'B1of2'},{'B2of2'}...
%   {'NT1of2For','NT1of2Rec'},{'NT2of2For','NT2of2Rec'}...
%   {'T1of2For','T1of2Rec'},{'T2of2For','T2of2Rec'}...
%   };
% exper.eventValuesExtra.newValue = {...
%   {'B1'},{'B2'}...
%   {'NT1'},{'NT2'}...
%   {'T1'},{'T2'}...
%   };

% B and T/NT x R/F (collapse across 1/2)
exper.eventValuesExtra.toCombine = {...
  {'B1of2','B2of2'}...
  {'NT1of2For','NT2of2For'},{'NT1of2Rec','NT2of2Rec'}...
  {'T1of2For','T2of2For'},{'T1of2Rec','T2of2Rec'}...
  };
exper.eventValuesExtra.newValue = {...
  {'B'}...
  {'NTF'},{'NTR'}...
  {'TF'},{'TR'}...
  };

% % B, NT, T (collapse across 1/2 and R/F)
% exper.eventValuesExtra.toCombine = {...
%   {'B1of2','B2of2'}...
%   {'NT1of2For','NT2of2For','NT1of2Rec','NT2of2Rec'},...
%   {'T1of2For','T2of2For','T1of2Rec','T2of2Rec'}...
%   };
% exper.eventValuesExtra.newValue = {...
%   {'B'}...
%   {'NT'}...
%   {'T'}...
%   };

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 1;

% From Nick: Data: /Volumes/RAID/curranlab/Data/TNT/TNT_matt. Here are the
% subjects I've used in the latest analysis (note that unused subjects are
% commented out):
exper.subjects = {
  % 'TNT 06';
   'TNT 07';
  % 'TNT 08';
   'TNT 09';
  % 'TNT 11';
  % 'TNT 13';
   'TNT 14';
   'TNT 15';
   'TNT 17';
   'TNT 19';
   'TNT 20';
   'TNT 21';
   'TNT 22';
   'TNT 23';
   'TNT 25';
   'TNT 26';
   'TNT 27';
   'TNT 28';
   'TNT 30';
   'TNT 32';
   'TNT 33';
   'TNT 35';
  % 'TNT 39';
   'TNT 41';
   'TNT 42';
   'TNT 44';
   'TNT 45';
  % 'TNT 46';
   'TNT 47';
  % 'TNT 48';
   'TNT 49';
  % 'TNT 50';
   'TNT 51';
   'TNT 52';
   'TNT 53';
   'TNT 54';
   };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'ses1'};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.dataDir = fullfile(exper.name,'TNT_matt','eeg',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));

% if save directory is different from read directory, can set
% dirs.saveDirStem; note that this currently needs to exist on
% dirs.dataroot, which is chosen below (i.e., you can't currently read from
% the server and save to your local computer)

% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile('/Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile('/Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile('/data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'data');

% pick the right dirs.dataroot; note the order of searching
if exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
  %runLocally = 1;
elseif exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
  %runLocally = 1;
elseif exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
  %runLocally = 0;
elseif exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

% Use the FT chan locs file
files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

% figure printing options - see mm_ft_setSaveDirs for other options
files.saveFigs = 1;
%files.figPrintFormat = 'png';
files.figPrintFormat = 'epsc2';

%% Convert the data to FieldTrip structs

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_timelockanalysis';
%ana.artifact.type = {'nsAuto','ftManual'};
ana.artifact.type = {'nsAuto'};

% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'tla';

% any preprocessing?
cfg_pp = [];
% single precision to save space
%cfg_pp.precision = 'single';

% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];

% % do a bandstop filter
% cfg_pp.bsfilter = 'yes';
% % % delta
% % cfg_pp.bsfreq = [1 3];
% % % theta
% % cfg_pp.bsfreq = [3 8];
% % % alpha
% % cfg_pp.bsfreq = [8 12];
% % all
% cfg_pp.bsfreq = [1 12];

% % do a bandpass filter
% cfg_pp.bpfilter = 'yes';
% cfg_pp.bpfreq = [1 40];

% do a lowpass filter
cfg_pp.lpfilter = 'yes';
cfg_pp.lpfreq = 40;
%cfg_pp.lpfreq = 1;

cfg_proc = [];
% do we want to keep the individual trials?
cfg_proc.keeptrials = 'no';

% set the save directories; final argument is prefix of save directory
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);
process_ft_data(ana,cfg_proc,exper,dirs);

%% save the analysis details

saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  fprintf('Done.\n');
else
  error('Not saving! %s already exists.\n',saveFile);
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldTrip format creation ends here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldTrip analysis starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the analysis details

% power
%adFile = '/Volumes/curranlab/Data/TNT/TNT_matt/eeg/-1000_1700/ft_data/NT_TH_eq1/pow_mtmconvol_hanning_pow_-500_980_3_9_avg/analysisDetails.mat';

% % ERP: T/NT/B x 1/2 x R/F
% adFile = '/Volumes/curranlab/Data/TNT/TNT_matt/eeg/-1000_1700/ft_data/B1of2_B2of2_NT1of2For_NT1of2Rec_NT2of2For_NT2of2Rec_T1of2For_T1of2Rec_T2of2For_T2of2Rec_eq0_art_nsAuto/tla_-1000_1700_avg/analysisDetails.mat';

% % ERP: T/NT/B x 1/2
% adFile = '/Volumes/curranlab/Data/TNT/TNT_matt/eeg/-1000_1700/ft_data/B1_B2_NT1_NT2_T1_T2_eq0_art_nsAuto/tla_-1000_1700_avg/analysisDetails.mat';

% ERP: T/NT/B x R/F
adFile = '/Volumes/curranlab/Data/TNT/TNT_matt/eeg/-1000_1700/ft_data/B_NTF_NTR_TF_TR_eq0_art_nsAuto/tla_-1000_1700_avg/analysisDetails.mat';

% % ERP: T/NT/B
% adFile = '/Volumes/curranlab/Data/TNT/TNT_matt/eeg/-1000_1700/ft_data/B_NT_T_eq0_art_nsAuto/tla_-1000_1700_avg/analysisDetails.mat';

[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% list the event values to analyze; specific to each experiment

% this is useful for when there are multiple types of event values, for
% example, hits and CRs in two conditions. You don't have to enter anything
% if you just want all events from exper.eventValues together in a single
% cell because it will get set to {exper.eventValues}, but it needs to be a
% cell containing a cell of eventValue strings

% this is used by mm_ft_checkCondComps to create pairwise combinations
% either within event types {'all_within_types'} or across all event types
% {'all_across_types'}; mm_ft_checkCondComps is called within subsequent
% analysis functions

ana.eventValues = {exper.eventValues};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end

%% load in the subject data

[data_tla] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'tla');

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

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs(exper,ana,15,false);

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = ana.eventValues;
cfg_ana.data_str = 'data_tla';
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

ga_tla = struct;

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sessions)
  for typ = 1:length(ana.eventValues)
    for evVal = 1:length(ana.eventValues{typ})
      %tic
      fprintf('Running ft_timelockgrandaverage on %s...',ana.eventValues{typ}{evVal});
      ga_tla.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
      fprintf('Done.\n');
      %toc
    end
  end
end

%% write out some data

cfg = [];

% cfg.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg.latencies = [0.3 0.5; 0.5 0.8];
% cfg.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg.latencies = [0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];

% cfg.rois = {{'LPS'},{'RPS'}};
% cfg.latencies = [0.5 0.8; 0.5 0.8];
cfg.rois = {{'LPS','RPS'}};
cfg.latencies = [0.5 0.8];
cfg.condByROI = repmat(ana.eventValues,size(cfg.rois));

cfg.parameter = 'avg';

cfg.excludeBadSub = true;

cfg.direction = 'columns';
%cfg.direction = 'rows';

mm_printDataToText(cfg,exper,ana,dirs,data_tla);

%% plot the conditions - simple

cfg_ft = [];
cfg_ft.parameter = 'avg';

cfg_plot = [];

% recognition memory effects
%cfg_ft.xlim = [-1.0 1.7];
cfg_ft.xlim = [-0.5 1.0];
cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
%cfg_plot.ylims = [-4 4; -1 6];
cfg_plot.ylims = [-4 4; -1.5 3.5];
cfg_plot.legendlocs = {'SouthEast','NorthWest'};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-4 4; -4 4; -1 6; -1 6];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','NorthWest','NorthWest'};

% % N170
% cfg_ft.xlim = [-.2 0.4];
% cfg_plot.rois = {{'E64'},{'E95'}};
% cfg_plot.ylims = [-4 4; -1 6];
% cfg_plot.legendlocs = {'SouthEast','NorthWest'};

cfg_plot.is_ga = 1;
% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
%cfg_plot.condByROI = repmat({{'TH','NT','B'}},size(cfg_plot.rois));
% cfg_plot.condByROI = {...
%   {{'TH','NT','B'}},...
%   {{'TH','NT','B'}}};

cfg_plot.condByROI = repmat({{'all'}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'T1of2Rec','T1of2For','NT1of2Rec','NT1of2For'}},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,ga_tla);
end

%% subplots of each subject's ERPs

cfg_plot = [];
cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
%cfg_plot.rois = {{'LAS','RAS'}};
%cfg_plot.roi = {'E124'};
%cfg_plot.roi = {'RAS'};
%cfg_plot.roi = {'LPS','RPS'};
%cfg_plot.roi = {'LPS'};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-.2 1.0];
cfg_plot.ylim = [-8 3];

cfg_plot.parameter = 'avg';

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_plot.condByROI = repmat({{'all'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'TH','NT','B'}},size(cfg_plot.rois));
% cfg_plot.condByROI = {...
%   {{'TH','NT','B'}},...
%   {{'TH','NT','B'}}};

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotER(cfg_plot,ana,exper,data_tla);
end

%% plot the conditions

cfg_ft = [];
cfg_ft.xlim = [-0.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

cfg_plot.ftFxn = 'ft_singleplotER';
%cfg_plot.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
cfg_plot.rois = {{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-4 1; -4 1; -4 1; -1.5 3.5; -1.5 3.5; -1.5 3.5];
% vertical solid lines to plot
cfg_plot.x_bounds = [0.3 0.4; 0.3 0.4; 0.3 0.4; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.plotLegend = 1;
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
%cfg_plot.plotTitle = 1;

% cfg_plot.xlabel = 'Time (s)';
% cfg_plot.ylabel = 'Voltage (\muV)';
% cfg_plot.xlabel = '';
% cfg_plot.ylabel = '';

% cfg_plot.ftFxn = 'ft_topoplotER';
% cfg_plot.plotLegend = 0;
% cfg_plot.plotTitle = 1;
% cfg_plot.ylims = [-6 6];
% %cfg_plot.ylims = 'maxmin';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% %cfg_ft.comment = 'no';
% cfg_plot.rois = {'all'};
% %cfg_ft.xlim = [0 1.2]; % time
% cfg_ft.xlim = [0 0.3]; % time
% cfg_plot.subplot = 1;
% % cfg_plot.rois = {{'LAS'}};
% % cfg_ft.xlim = [0.3 0.5]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [0.5 0.8]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [1.0 1.5]; % time

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_plot.plotLegend = 0;
% cfg_plot.plotTitle = 1;
% cfg_ft.showlabels = 'yes';
% cfg_ft.interactive = 'yes';
% cfg_ft.comment = '';
% cfg_ft.xlim = [-0.1 0.6];
% %cfg_plot.rois = {{'FS'},{'LAS','RAS'},{'LPS','RPS'}};
% %cfg_plot.rois = {{'FS'},{'PS'}};
% %cfg_plot.rois = {'E71'};
% cfg_plot.rois = {'all'};
% cfg_plot.ylims = [-4 4]; % voltage in multiplot
% %cfg_plot.ylims = repmat('maxmin',size(cfg_plot.rois,2),1); % voltage in multiplot

% INSTRUCTIONS: outermost cell holds one cell for each ROI; each ROI cell
% holds one cell for each event type; each event type cell holds strings
% for its conditions

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'B','TR'}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'TH','NT','B'}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'T1of2Rec','T1of2For','NT1of2Rec','NT1of2For'}},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  if strcmp(cfg_plot.ftFxn,'ft_singleplotER')
    cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
  end
  if cfg_plot.plotLegend
    cfg_plot.legendloc = cfg_plot.legendlocs{r};
  end

  mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 1;

% comparisons to make
%cfg_plot.conditions = {{'T','NT'},{'T','B'},{'NT','B'}};
cfg_plot.conditions = {{'TR','TF'},{'NTR','NTF'},{'TR','NTR'},{'TF','NTF'},{'TR','B'},{'TF','B'},{'NTR','B'},{'NTF','B'}};
%cfg_plot.conditions = {{'T1','T2'},{'NT1','NT2'},{'T1','NT1'},{'T2','NT2'},{'T1','B1'},{'T2','B2'},{'NT1','B1'},{'NT2','B2'}};
%cfg_plot.conditions = {{'TR','NTF'}};
%cfg_plot.conditions = {'all'};

cfg_ft = [];
cfg_ft.xlim = [-0.2 1]; % time
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'yes';

cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-1.5 1.5]; % volt
cfg_ft.marker = 'on';
%cfg_ft.markerfontsize = 9;

%cfg_ft.comment = 'no';
cfg_plot.roi = {'LAS','RAS'};
cfg_ft.xlim = [0.3 0.4]; % time
cfg_ft.xlim = [0.2 0.6]; % time
% cfg_plot.roi = {'LPS','RPS'};
% cfg_ft.xlim = [0.5 0.8]; % time

cfg_plot.subplot = 1;
% cfg_ft.fontsize = 6;
% cfg_ft.xlim = [0 1.0]; % time
% %cfg_ft.xlim = (0:0.05:1.0); % time
% %cfg_plot.roi = {'PS'};
cfg_plot.numCols = 3;

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.ylim = [-1.5 1.5]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.xlim = [-0.2 1]; % time
% cfg_ft.ylim = [-1.5 1.5]; % volt

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);

%% descriptive statistics: ttest

cfg_ana = [];

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.4; 0.5 0.8];

%cfg_ana.conditions = {{'T','NT'},{'T','B'},{'NT','B'}};

%cfg_ana.conditions = {{'TR','TF'},{'NTR','NTF'},{'TR','NTR'},{'TF','NTF'},{'TR','B'},{'NTR','B'},{'TF','B'},{'NTF','B'}};
cfg_ana.conditions = {'all'};

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
cfg_plot.ylims = [-4 -1; 2 5];
cfg_plot.plot_order = {'B','TH','NT'};

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla);
end

%% 3-way ANOVA: T/NT vs 1/2 vs R/F

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 0;

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% % IV2: define the conditions tested for each set of ROIs
% %cfg_ana.condByROI = {{'RH','RCR'},{'RCR','RHSC','RHSI'}};
% cfg_ana.condByROI = {{'TH','NT','B'},{'TH','NT','B'}};

cfg_ana.condByROI = {...
  {'NT1of2For', 'NT1of2Rec', 'NT2of2For', 'NT2of2Rec', 'T1of2For', 'T1of2Rec', 'T2of2For', 'T2of2Rec'},...
  {'NT1of2For', 'NT1of2Rec', 'NT2of2For', 'NT2of2Rec', 'T1of2For', 'T1of2Rec', 'T2of2For', 'T2of2Rec'}};

% Define the IVs (type: event, roi, latency)
cfg_ana.IV1.name = 'T/NT';
cfg_ana.IV1.cond = {'T','NT'};
cfg_ana.IV1.type = 'event';

cfg_ana.IV2.name = 'ExperHalf';
cfg_ana.IV2.cond = {'1of2','2of2'};
cfg_ana.IV2.type = 'event';

cfg_ana.IV3.name = 'MemoryAcc';
cfg_ana.IV3.cond = {'Rec','For'};
cfg_ana.IV3.type = 'event';

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  
  mm_ft_rmaov33ER_spec(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Condition x Half

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 0;

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.4; 0.5 0.8];

%cfg_ana.condByROI = repmat({{'NT1of2', 'NT2of2', 'T1of2', 'T2of2', 'B1of2', 'B2of2'}},size(cfg_ana.rois));
cfg_ana.condByROI = repmat({{'NT1', 'NT2', 'T1', 'T2', 'B1', 'B2'}},size(cfg_ana.rois));
% Define the IVs (type: event, roi, latency)
cfg_ana.IV1.name = 'T/NT/B';
cfg_ana.IV1.cond = {'T','NT','B'};
cfg_ana.IV1.type = 'event';
cfg_ana.IV2.name = 'ExperHalf';
cfg_ana.IV2.cond = {'1','2'};
cfg_ana.IV2.type = 'event';

% cfg_ana.condByROI = repmat({{'TR', 'TF', 'NTR', 'NTF'}},size(cfg_ana.rois));
% % Define the IVs (type: event, roi, latency)
% cfg_ana.IV1.name = 'T/NT';
% cfg_ana.IV1.cond = {'T','NT'};
% cfg_ana.IV1.type = 'event';
% cfg_ana.IV2.name = 'R/F';
% cfg_ana.IV2.cond = {'R','F'};
% cfg_ana.IV2.type = 'event';

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  
  mm_ft_rmaov2ER_spec(cfg_ana,exper,ana,data_tla);
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
cfg_ana.condByROI = {{'TH','NT','B'},{'TH','NT','B'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% run the cluster statistics

cfg_ft = [];
cfg_ft.avgoverchan = 'no';
cfg_ft.avgovertime = 'no';

cfg_ft.parameter = 'avg';

cfg_ana = [];
% cfg_ana.roi = 'all';
% cfg_ana.latencies = [0 1.0; 1.0 1.7];
cfg_ana.roi = {'anterior'};
cfg_ana.latencies = [0.2 0.6];

%cfg_ana.conditions = {{'TH','NT'},{'TH','B'},{'NT','B'}};
%cfg_ana.conditions = {'all'};

% cfg_ana.conditions = {...
%   {'T1of2Rec', 'T1of2For'}, {'T2of2Rec', 'T2of2For'}, {'NT1of2Rec', 'NT1of2For'},{'NT2of2Rec', 'NT2of2For'}...
%   {'T1of2Rec', 'NT1of2Rec'},{'T1of2For', 'NT1of2For'}, {'T2of2Rec', 'NT2of2Rec'},{'T2of2For', 'NT2of2For'}...
%   {'T1of2Rec', 'B1of2'},{'T1of2For', 'B1of2'}, {'T2of2Rec', 'B2of2'},{'T2of2For', 'B2of2'}...
%   {'NT1of2Rec', 'B1of2'},{'NT1of2For', 'B1of2'}, {'NT2of2Rec', 'B2of2'},{'NT2of2For', 'B2of2'}...
%   };

%cfg_ana.conditions = {{'T1of2Rec', 'NT1of2For'}, {'T2of2Rec', 'NT2of2For'}};

cfg_ana.conditions = {{'TR','TF'},{'NTR','NTF'},{'TR','NTR'},{'TF','NTF'},{'TR','B'},{'NTR','B'},{'TF','B'},{'NTF','B'}};

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  
  stat_clus = mm_ft_clusterstatER(cfg_ft,cfg_ana,exper,ana,dirs,data_tla);
end

%% plot the cluster statistics

files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = 0.1;

cfg_plot = [];
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.conditions = cfg_ana.conditions;

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  
  mm_ft_clusterplotER(cfg_ft,cfg_plot,ana,files,dirs);
end

%% correlations

cfg_ana = [];

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.dpTypesByROI = {...
  {'Item'},...
  {'Item'}};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds two cells, one for each
% d' type; each d' cell contains strings for its conditions
cfg_ana.condByROI = {...
  {{'TH','NT'}},...
  {{'TH','NT'}}};

% C2 d' values
cfg_ana.d_item =  abs([]);
cfg_ana.d_source =  abs([]);

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data_tla);
end
