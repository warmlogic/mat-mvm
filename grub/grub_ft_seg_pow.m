% Make plots and do analyses for time-frequency EEG data

% See Maris and Oostenveld (2007) for a description of statistics

% if strcmp(dirs.dataroot,dirs.serverLocalDir)
%   dirs.dataroot = strrep(dirs.dataroot,'RAID/','');
%   dirs.saveDirProc = strrep(dirs.saveDirProc,'RAID/','');
%   dirs.saveDirFigs = strrep(dirs.saveDirFigs,'RAID/','');
% end

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'GRUB';

exper.sampleRate = 500;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 2.0];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
exper.eegFileExt = 'egis';
%exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
exper.eventValues = sort({'RCR','RHSC','RHSI'});

% combine the two types of hits into one category
exper.eventValuesExtra.toCombine = {{'RHSC','RHSI'}};
exper.eventValuesExtra.newValue = {{'RH'}};

exper.subjects = {
  'GRUB002';
  'GRUB003';
  'GRUB004';
  'GRUB005';
  'GRUB006';
  'GRUB007';
  'GRUB008';
  'GRUB009';
  'GRUB010';
  'GRUB011';
  'GRUB012';
  'GRUB014';
  'GRUB015';
  'GRUB016';
  'GRUB017';
  'GRUB018';
  'GRUB019';
  'GRUB020';
  'GRUB021';
  'GRUB022';
  'GRUB023';
  'GRUB024';
  'GRUB025';
  'GRUB026';
  };
%  'GRUB001'; % exper was not set up correctly
%  'GRUB013'; % only half the session was recorded

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_0'};

%% set up parameters

% directory where the data to read is located
dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));

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
files.figPrintFormat = 'png';
%files.figPrintFormat = 'epsc2';

%% add NS's artifact information to the event structure
% nsEvFilters.eventValues = exper.eventValues;
% % RCR
% nsEvFilters.RCR.type = 'LURE_PRES';
% nsEvFilters.RCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % RHSC
% nsEvFilters.RHSC.type = 'TARG_PRES';
% nsEvFilters.RHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % RHSI
% nsEvFilters.RHSI.type = 'TARG_PRES';
% nsEvFilters.RHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% 
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sessions)
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,0);
%   end
% end

%% Convert the data to FieldTrip structs

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_freqanalysis';
ana.artifact.type = 'nsAuto';

% any preprocessing?
cfg_pp = [];
% single precision to save space
cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.output = 'pow';
cfg_proc.pad = 'maxperlen';
cfg_proc.keeptrials = 'no';
cfg_proc.keeptapers = 'no';

% % MTM FFT
% cfg_proc.method = 'mtmfft';
% cfg_proc.taper = 'dpss';
% %cfg_proc.foilim = [3 50];
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% %cfg_proc.foi = 3:freqstep:50;
% cfg_proc.foi = 3:freqstep:9;
% cfg_proc.tapsmofrq = 5;
% cfg_proc.toi = -0:0.04:1.0;

% multi-taper method
cfg_proc.method = 'mtmconvol';
cfg_proc.taper = 'hanning';
%cfg_proc.taper = 'dpss';
%cfg_proc.toi = -0.8:0.04:3.0;
cfg_proc.toi = -0.5:0.04:1.0;
freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
cfg_proc.foi = 3:freqstep:40;
%cfg_proc.foi = 3:freqstep:9;
%cfg_proc.foi = 3:1:9;
%cfg_proc.foi = 2:2:30;
cfg_proc.t_ftimwin = 5./cfg_proc.foi;
% tapsmofrq is not used for hanning taper; it is used for dpss
%cfg_proc.tapsmofrq = 0.4*cfg_proc.foi;

% % wavelet
% cfg_proc.method = 'wavelet';
% cfg_proc.width = 5;
% %cfg_proc.toi = -0.8:0.04:3.0;
% cfg_proc.toi = -0.3:0.04:1.0;
% % evenly spaced frequencies, but not as many as foilim makes
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% cfg_proc.foi = 3:freqstep:40;
% %cfg_proc.foi = 3:freqstep:9;
% %cfg_proc.foilim = [3 9];

% set the save directories; final argument is prefix of save directory
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'pow');

% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = cfg_proc.output;

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,cfg_proc,exper,dirs,files);

%% save the analysis details

% overwrite if it already exists
saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
%if ~exist(saveFile,'file')
fprintf('Saving %s...',saveFile);
save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
fprintf('Done.\n');
%else
%  error('Not saving! %s already exists.\n',saveFile);
%end

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

adFile = '/Volumes/curranlab/Data/GRUB/eeg/eppp/-1000_2000/ft_data/RCR_RHSC_RHSI_RH_eq1/pow_mtmconvol_hanning_pow_-500_980_3_40_avg/analysisDetails.mat';
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

%% set up channel groups

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

ana.eventValues = {exper.eventValues};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end

%% load some data

[data_freq] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

%% Test plots to make sure data look ok

cfg_ft = [];
cfg_ft.baseline = [-0.5 -0.1];
cfg_ft.baselinetype = 'absolute';
if strcmp(cfg_ft.baselinetype,'absolute') % data - mean
  %cfg_ft.zlim = [-400 400];
  cfg_ft.zlim = [-2 2];
elseif strcmp(cfg_ft.baselinetype,'relative') % ratio (1 = no change)
  cfg_ft.zlim = [0 2];
elseif strcmp(cfg_ft.baselinetype,'relchange') % (data - mean) / mean
  cfg_ft.zlim = [-1 1];
end
cfg_ft.parameter = 'powspctrm';
%cfg_ft.ylim = [3 9];
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
sub=1;
ses=1;
for i = 1:2
  figure
  ft_multiplotTFR(cfg_ft,data_freq.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

% cfg_ft = [];
% cfg_ft.channel = {'E124'};
% % %cfg_ft.channel = {'E117'};
% cfg_ft.baseline = [-0.3 -0.1];
% cfg_ft.baselinetype = 'absolute';
% % if strcmp(cfg_ft.baselinetype,'absolute')
% %   %cfg_ft.zlim = [-2 2];
% %   cfg_ft.zlim = [-500 500];
% % elseif strcmp(cfg_ft.baselinetype,'relative')
% %   cfg_ft.zlim = [0 1.5];
% % end
% % cfg_ft.showlabels = 'yes';
% % cfg_ft.colorbar = 'yes';
% % cfg_ft.ylim = [4 8];
% figure
% ft_singleplotTFR(cfg_ft,data_freq.(exper.eventValues{1}).sub(1).ses(1).data);

%% Change in freq relative to baseline using absolute power

cfg_fb = [];
%cfg_fb.baseline = [-0.3 -0.1];
cfg_fb.baseline = [-0.5 -0.1];
cfg_fb.baselinetype = 'absolute';

%data_freq_orig = data_freq;

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        fprintf('%s, %s, %s, ',exper.subjects{sub},exper.sessions{ses},ana.eventValues{typ}{evVal});
        data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
      end
    end
  end
end

% % find the time points without NaNs for a particular frequency
% data_freq.(ana.eventValues{1}{1}).sub(1).ses(1).data.time(~isnan(squeeze(data_freq.(ana.eventValues{1}{1}).sub(1).ses(1).data.powspctrm(1,2,:))'))
% ga_freq.(ana.eventValues{1}{1}).time(~isnan(squeeze(ga_freqanaexper.eventValues{1}{1}).powspctrm(1,1,2,:))'))

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {};

% exclude subjects with low event counts
[exper] = mm_threshSubs(exper,ana,20);

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = ana.eventValues;
cfg_ana.data_str = 'data_freq';
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sessions)
  for typ = 1:length(ana.eventValues)
    for evVal = 1:length(ana.eventValues{typ})
      %tic
      fprintf('Running ft_freqgrandaverage on %s...',ana.eventValues{typ}{evVal});
      ga_freq.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
      fprintf('Done.\n');
      %toc
    end
  end
end

%% plot the conditions - simple

cfg_ft = [];
%cfg_ft.baseline = [-0.3 -0.1];
%cfg_ft.baselinetype = 'absolute';
%if strcmp(cfg_ft.baselinetype,'absolute')
cfg_ft.xlim = [0 1];
%cfg_ft.ylim = [3 50];
%cfg_ft.ylim = [3 8];
%cfg_ft.ylim = [8 12];
cfg_ft.zlim = [-2 2];
%cfg_ft.zlim = [-100 100];
%elseif strcmp(cfg_ft.baselinetype,'relative')
%  cfg_ft.zlim = [0 2.0];
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
for typ = 1:length(ana.eventValues)
  for evVal = 1:length(ana.eventValues{typ})
    figure
    ft_multiplotTFR(cfg_ft,ga_freq.(ana.eventValues{typ}{evVal}));
    set(gcf,'Name',sprintf('%s',ana.eventValues{typ}{evVal}))
  end
end

%% subplots of each subject's power spectrum

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.rois = {{'FS'},{'PS'}};
%cfg_plot.roi = {'E124'};
%cfg_plot.roi = {'RAS'};
%cfg_plot.roi = {'LPS','RPS'};
%cfg_plot.roi = {'LPS'};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_plot.condByROI = repmat({{'RCR','RH','RHSC','RHSI'}},size(cfg_plot.rois));

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.zlim = [-2 2];
cfg_ft.parameter = 'powspctrm';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotTFR(cfg_ft,cfg_plot,ana,exper,data_freq);
end

%% make some GA plots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [.5 1.0]; % time
cfg_ft.ylim = [3 8]; % freq
%cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 40]; % freq
cfg_ft.zlim = [-2 2]; % pow

cfg_ft.parameter = 'powspctrm';

cfg_plot = [];
cfg_plot.plotTitle = 1;

%cfg_plot.rois = {{'FS'},{'LAS','RAS'},{'LPS','RPS'}};
%cfg_plot.rois = {{'FS'},{'PS'}};
%cfg_plot.rois = {'E71'};
cfg_plot.rois = {'all'};

cfg_plot.is_ga = 1;
% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_plot.condByROI = repmat({{'RCR','RH','RHSC','RHSI'}},size(cfg_plot.rois));

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

%cfg_plot.ftFxn = 'ft_singleplotTFR';

% cfg_plot.ftFxn = 'ft_topoplotTFR';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.comment = 'no';
% cfg_ft.xlim = [0.5 0.8]; % time
% %cfg_plot.subplot = 1;
% %cfg_ft.xlim = [0 1.0]; % time

cfg_plot.ftFxn = 'ft_multiplotTFR';
cfg_ft.showlabels = 'yes';
cfg_ft.comment = '';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_freq);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 1;

% comparisons to make
cfg_plot.conditions = {'all'};

cfg_ft = [];
%cfg_ft.xlim = [.5 .8]; % time
%cfg_ft.ylim = [3 8]; % freq
%cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.parameter = 'powspctrm';
cfg_ft.zlim = [-1 1]; % pow

cfg_ft.interactive = 'yes';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'yes';

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

%cfg_plot.ftFxn = 'ft_singleplotTFR';
cfg_plot.ftFxn = 'ft_topoplotTFR';
%cfg_ft.marker = 'on';
cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.comment = 'no';
%cfg_ft.xlim = [0.5 0.8]; % time
cfg_plot.subplot = 1;
cfg_ft.xlim = [0 1.0]; % time
%cfg_ft.xlim = (0:0.05:1.0); % time
%cfg_plot.roi = {'PS'};

% cfg_plot.ftFxn = 'ft_multiplotTFR';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';

mm_ft_contrastTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_freq);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [4 8; 4 8];

%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};
cfg_ana.conditions = {'all'};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.avgoverfreq = 'yes';
cfg_ft.parameter = 'powspctrm';
cfg_ft.correctm = 'fdr';

cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
% line plot parameters
cfg_plot.ylims = [-1 1; -1 1];
%cfg_plot.ylims = [-100 100; -100 100];

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_freq);
end

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LAS','RAS'},{'LPS','RPS'}};
% IV2: define the conditions tested for each set of ROIs
cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
%cfg_ana.condByROI = repmat({{'TH','NT'}},size(cfg_ana.rois));
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [3 8; 3 8; 8 12];

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'powspctrm';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.frequency = cfg_ana.frequencies(r,:);
  
  mm_ft_rmaov2TFR(cfg_ana,exper,ana,data_freq);
end

%% cluster statistics

cfg_ft = [];
cfg_ft.avgoverchan = 'no';
cfg_ft.avgovertime = 'no';
cfg_ft.avgoverfreq = 'yes';
%cfg_ft.avgoverfreq = 'no';

cfg_ft.parameter = 'powspctrm';

% debugging
%cfg_ft.numrandomization = 100;

% use .025 for a two-tailed test
cfg_ft.alpha = 0.025;

cfg_ana = [];
cfg_ana.roi = 'all';
%cfg_ana.conditions = {'all'};
cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};

%cfg_ana.frequencies = [3 8; 8 12; 12 28; 28 50; 50 100];
cfg_ana.frequencies = [3 8; 8 12; 12 28; 28 40];
%cfg_ana.frequencies = [3 40];
cfg_ana.latencies = [0 1.0];
%cfg_ana.latencies = [0 0.5; 0.5 1.0];

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  for fr = 1:size(cfg_ana.frequencies,1)
    cfg_ft.frequency = cfg_ana.frequencies(fr,:);
    
    [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_freq);
  end
end

%% plot the cluster statistics

files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = .05;

cfg_plot = [];
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.frequencies = cfg_ana.frequencies;
cfg_plot.latencies = cfg_ana.latencies;

% not averaging over frequencies - only works with ft_multiplotTFR
files.saveFigs = 0;
cfg_ft.avgoverfreq = 'no';
cfg_ft.interactive = 'yes';
cfg_plot.mask = 'yes';
%cfg_ft.maskstyle = 'saturation';
cfg_ft.maskalpha = 0.3;
cfg_plot.ftFxn = 'ft_multiplotTFR';
% http://mailman.science.ru.nl/pipermail/fieldtrip/2009-July/002288.html
% http://mailman.science.ru.nl/pipermail/fieldtrip/2010-November/003312.html

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  for fr = 1:size(cfg_plot.frequencies,1)
    cfg_ft.frequency = cfg_plot.frequencies(fr,:);
    
    mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs);
  end
end

%% correlations

cfg_ana = [];

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
cfg_ana.frequencies = [3 8; 3 8];

cfg_ana.dpTypesByROI = {...
  {'Item','Source'},...
  {'Item','Source'}};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds two cells, one for each
% d' type; each d' cell contains strings for its conditions
cfg_ana.condByROI = {...
  {{'RCR','RH'},{'RHSC','RHSI'}}...
  {{'RCR','RH'},{'RHSC','RHSI'}}};

% d' values
d_item = abs([0.833 2.2676 1.6705 2.9361 2.1089 1.9395 2.2178 2.1909 2.0847 2.4001 1.2743 1.4276 1.1561 2.606 1.7639 2.324 3.1678 1.5034 1.9145 2.0296 2.5468 2.3279 2.4804 1.6196]);
d_source = abs([0.2474 1.1415 1.2173 2.0237 1.55 0.9553 0.755 1.3439 1.1001 1.2752 0.4722 0.6733 0.23 0.9382 0.4262 2.014 2.0981 0.7012 1.2259 1.6901 1.3517 1.3498 1.5733 1.1314]);

cfg_ana.parameter = 'powspctrm';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.frequency = cfg_ana.frequencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data_tla);
end

%% Make contrast plots (with culster stat info) - old function

% set up contrast
cfg_ana = [];
cfg_ana.include_clus_stat = 0;
cfg_ana.timeS = (0:0.05:1.0);
cfg_ana.timeSamp = round(linspace(1,exper.sampleRate,length(cfg_ana.timeS)));

cfg_plot = [];
cfg_plot.minMaxVolt = [-1 1];
cfg_plot.numRows = 4;

cfg_ft = [];
cfg_ft.interactive = 'no';
cfg_ft.elec = ga_freq.elec;
cfg_ft.highlight = 'on';
if cfg_ana.include_clus_stat == 0
  cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS','RPS','LPS'})});
end
%cfg_ft.comment = 'xlim';
cfg_ft.commentpos = 'title';

% create contrast
cont_topo = [];
cont_topo.RHvsRCR = ga_freq.RH;
cont_topo.RHvsRCR.avg = ga_freq.RH.avg - ga_freq.RCR.avg;
cont_topo.RHvsRCR.individual = ga_freq.RH.individual - ga_freq.RCR.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.RHvsRCR.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.RHvsRCR);
end
set(gcf,'Name','H - CR')

% create contrast
cont_topo.RHSCvsRHSI = ga_freq.RHSC;
cont_topo.RHSCvsRHSI.avg = ga_freq.RHSC.avg - ga_freq.RHSI.avg;
cont_topo.RHSCvsRHSI.individual = ga_freq.RHSC.individual - ga_freq.RHSI.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.RHSCvsRHSI.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.RHSCvsRHSI);
end
set(gcf,'Name','HSC - HSI')

% create contrast
cont_topo.RHSCvsRCR = ga_freq.RHSC;
cont_topo.RHSCvsRCR.avg = ga_freq.RHSC.avg - ga_freq.RCR.avg;
cont_topo.RHSCvsRCR.individual = ga_freq.RHSC.individual - ga_freq.RCR.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.RHSCvsRCR.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.RHSCvsRCR);
end
set(gcf,'Name','HSC - CR')

% % mecklinger plot
% cfg_ft.xlim = [1200 1800];
% cfg_ft.zlim = [-2 2];
% cfg_ft.colorbar = 'yes';
% figure
% ft_topoplotER(cfg_ft,cont_topo.RHSCvsRCR);

% create contrast
cont_topo.RHSIvsRCR = ga_freq.RHSI;
cont_topo.RHSIvsRCR.avg = ga_freq.RHSI.avg - ga_freq.RCR.avg;
cont_topo.RHSIvsRCR.individual = ga_freq.RHSI.individual - ga_freq.RCR.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.RHSIvsRCR.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.RHSIvsRCR);
end
set(gcf,'Name','HSI - CR')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % because of a bug (might be fixed now)
% if ~isfield(stat_clus.RHSCvsRHSIvsRCR,'negclusters') && isfield(stat_clus.RHSCvsRHSIvsRCR,'posclusters')
%   fprintf('No neg clusters found\n');
%   stat_clus.RHSCvsRHSIvsRCR.negclusters.prob = .5;
%   stat_clus.RHSCvsRHSIvsRCR.negclusters.clusterstat = 0;
%   stat_clus.RHSCvsRHSIvsRCR.negclusterslabelmat = zeros(size(stat_clus.RHSCvsRHSIvsRCR.posclusterslabelmat));
%   stat_clus.RHSCvsRHSIvsRCR.negdistribution = zeros(size(stat_clus.RHSCvsRHSIvsRCR.posdistribution));
% end
% if ~isfield(stat_clus.RHSCvsRHSIvsRCR,'posclusters') && isfield(stat_clus.RHSCvsRHSIvsRCR,'negclusters')
%   fprintf('No pos clusters found\n');
%   stat_clus.RHSCvsRHSIvsRCR.posclusters.prob = 1;
%   stat_clus.RHSCvsRHSIvsRCR.posclusters.clusterstat = 0;
%   stat_clus.RHSCvsRHSIvsRCR.posclusterslabelmat = zeros(size(stat_clus.RHSCvsRHSIvsRCR.negclusterslabelmat));
%   stat_clus.RHSCvsRHSIvsRCR.posdistribution = zeros(size(stat_clus.RHSCvsRHSIvsRCR.negdistribution));
% end
%
% cfg_ft = [];
% % p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
% cfg_ft.highlightsymbolseries = ['*','*','.','.','.'];
% cfg_ft.layout = ft_prepare_layout(cfg_ft,ga_freq);
% cfg_ft.contournum = 0;
% cfg_ft.emarker = '.';
% cfg_ft.alpha  = 0.05;
% cfg_ft.parameter = 'stat';
% cfg_ft.zlim = [-5 5];
% ft_clusterplot(cfg_ft,stat_clus.RHSCvsRHSIvsRCR);
