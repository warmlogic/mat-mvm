% Make plots and do analyses for time-frequency EEG data

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SOCO';

exper.sampleRate = 250;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 2.0];

% equate the number of trials across event values?
exper.equateTrials = 0;
%exper.eventValuesExtra.equateExtrasSeparately = 1;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
exper.nsFileExt = 'egis';
%exper.nsFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
exper.eventValues = sort({'CR2','HSC2','HSI2','CR6','HSC6','HSI6'});
%exper.eventValues = sort({'F2','F6','N2','N6','RO2','RO6','RS2','RS6'});

% combine some events into higher-level categories
exper.eventValuesExtra.toCombine = {{'HSC2','HSI2'},{'HSC6','HSI6'}};
exper.eventValuesExtra.newValue = {{'H2'},{'H6'}};
%exper.eventValuesExtra.toCombine = {{'F2','F6'},{'N2','N6'},{'RO2','RO6'},{'RS2','RS6'}};
%exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 0;

exper.subjects = {
  'SOCO001';
  'SOCO002';
  'SOCO003';
  'SOCO004';
  'SOCO005';
  'SOCO006';
  'SOCO007';
  'SOCO008';
  'SOCO009';
  'SOCO010';
  'SOCO011';
  'SOCO012';
  'SOCO013';
  'SOCO014';
  'SOCO015';
  'SOCO016';
  'SOCO017';
  'SOCO018';
  'SOCO019';
  'SOCO020';
  'SOCO021';
  'SOCO022';
  'SOCO023';
  'SOCO024';
  'SOCO025';
  'SOCO026';
  'SOCO027';
  'SOCO028';
  'SOCO029';
  'SOCO030';
  };
% SOCO002 ended early by 6(?) trials because of fire alarm

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not the
% session directory names where the FieldTrip data is saved for each
% subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_0'};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));
%dirs.subDir = 'RK';
%dirs.dataDir = fullfile('eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);

% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile('/Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile('/Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile('/data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'data');

% pick the right dirs.dataroot
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
files.figFileExt = 'png';

% %% add NS's artifact information to the event structure
% nsEvFilters.eventValues = exper.eventValues;
% % CR2
% nsEvFilters.CR2.type = 'LURE_PRES';
% nsEvFilters.CR2.testType = 'side';
% nsEvFilters.CR2.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % HSC2
% nsEvFilters.HSC2.type = 'TARG_PRES';
% nsEvFilters.HSC2.testType = 'side';
% nsEvFilters.HSC2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % HSI2
% nsEvFilters.HSI2.type = 'TARG_PRES';
% nsEvFilters.HSI2.testType = 'side';
% nsEvFilters.HSI2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% % CR6
% nsEvFilters.CR6.type = 'LURE_PRES';
% nsEvFilters.CR6.testType = 'task';
% nsEvFilters.CR6.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % HSC6
% nsEvFilters.HSC6.type = 'TARG_PRES';
% nsEvFilters.HSC6.testType = 'task';
% nsEvFilters.HSC6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % HSI6
% nsEvFilters.HSI6.type = 'TARG_PRES';
% nsEvFilters.HSI6.testType = 'task';
% nsEvFilters.HSI6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% 
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sessions)
%     overwriteEvFile = 1;
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,overwriteEvFile);
%   end
% end

%% Convert the data to FieldTrip structs

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_freqanalysis';
ana.artifactType = 'ns';

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
% freqstep = exper.sampleRate/(sum(abs(exper.prepost))*exper.sampleRate)*2;
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
freqstep = exper.sampleRate/(sum(abs(exper.prepost))*exper.sampleRate)*2;
cfg_proc.foi = 2:freqstep:40;
%cfg_proc.foi = 3:freqstep:9;
%cfg_proc.foi = 3:1:9;
%cfg_proc.foi = 2:2:30;
cfg_proc.t_ftimwin = 4./cfg_proc.foi;
% tapsmofrq is not used for hanning taper; it is used for dpss
%cfg_proc.tapsmofrq = 0.4*cfg_proc.foi;

% % wavelet
% cfg_proc.method = 'wavelet';
% cfg_proc.width = 5;
% %cfg_proc.toi = -0.8:0.04:3.0;
% cfg_proc.toi = -0.3:0.04:1.0;
% % evenly spaced frequencies, but not as many as foilim makes
% freqstep = exper.sampleRate/(sum(abs(exper.prepost))*exper.sampleRate)*2;
% %cfg_proc.foi = 3:freqstep:50;
% cfg_proc.foi = 3:freqstep:9;
% %cfg_proc.foilim = [3 9];

% set the save directories; final argument is prefix of save directory
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'pow');

% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = cfg_proc.output;

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,cfg_proc,exper,dirs,files);

%% set up channel groups

% pre-defined in this function
ana = mm_ft_channelgroups(ana);

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

% list the values separated by types: 2Colors, 6Colors
%ana.eventValues = {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}};
ana.eventValues = {{'F','N','RO','RS'}};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end

%% concatenate individual subject files

% fprintf('rsync -av matt@dream.colorado.edu:/data/projects/curranlab/%s/%s/*.mat %s/\n',dirs.dataDir,dirs.saveDirName,dirs.saveDir);

% [data_freq,exper,cfg_proc] = mm_ft_concatSubs_tfr(exper,dirs);

%% save the analysis details
saveFile = fullfile(dirs.saveDir,sprintf('analysisDetails_%s_%s%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','cfg_proc');
  fprintf('Done.\n');
else
  error('Not saving! %s already exists.\n',saveFile);
end

%% if already saved and not yet loaded, load the ft_freqanalysis files

if ~exist('cfg_proc','var')
  savedFiles = dir(fullfile(dirs.saveDir,'analysisDetails*.mat'));
  if length(savedFiles) == 1
    load(fullfile(dirs.saveDir,savedFiles.name));
  elseif length(savedFiles) > 1
    error('Multiple analysisDetails*.mat files found in %s!',dirs.saveDir)
  elseif isempty(savedFiles)
    error('analysisDetails*.mat not found in %s!',dirs.saveDir)
  end
end

if ~exist('data_freq','var')
  if strcmp(cfg_proc.keeptrials,'no')
    savedFiles = dir(fullfile(dirs.saveDir,sprintf('data_%s_%s%s_avg_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000)));
  elseif strcmp(cfg_proc.keeptrials,'yes')
    savedFiles = dir(fullfile(dirs.saveDir,sprintf('data_%s_%s%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000)));
  end
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDir,savedFiles(sf).name));
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

cfg_ft = [];
cfg_ft.baseline = [-0.3 -0.1];
cfg_ft.baselinetype = 'absolute';
if strcmp(cfg_ft.baselinetype,'absolute')
  cfg_ft.zlim = [-400 400];
  %cfg_ft.zlim = [-2 2];
elseif strcmp(cfg_ft.baselinetype,'relative')
  cfg_ft.zlim = [0 2.0];
end
cfg_ft.zparam = 'powspctrm';
%cfg_ft.ylim = [3 9];
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
figure
ft_multiplotTFR(cfg_ft,data_freq.(exper.eventValues{1}).sub(1).ses(1).data);
figure
ft_multiplotTFR(cfg_ft,data_freq.(exper.eventValues{2}).sub(1).ses(1).data);

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
cfg_fb.baseline = [-0.3 -0.1];
cfg_fb.baselinetype = 'absolute';

data_freq_orig = data_freq;

for evVal = 1:length(exper.eventValues)
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      fprintf('%s, %s, %s, ',exper.subjects{sub},exper.sessions{ses},exper.eventValues{evVal});
      data_freq.(exper.eventValues{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_freq.(exper.eventValues{evVal}).sub(sub).ses(ses).data);
    end
  end
end

% % save the structs for loading in later
% if strcmp(cfg_proc.keeptrials,'no')
%   saveFile = fullfile(dirs.saveDir,sprintf('data_%s_blc_%s%s_avg_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
% elseif strcmp(cfg_proc.keeptrials,'yes')
%   saveFile = fullfile(dirs.saveDir,sprintf('data_%s_blc_%s%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
% end
% if ~exist(saveFile,'file')
%   fprintf('Saving %s...',saveFile);
%   save(saveFile,'data_freq');
%   fprintf('Done.\n');
% else
%   error('Not saving! %s already exists.',saveFile);
% end

% % find the time points without NaNs for a particular frequency
% data_freq.(exper.eventValues{1}).sub(1).ses(1).data.time(~isnan(squeeze(data_freq.(exper.eventValues{1}).sub(1).ses(1).data.powspctrm(1,2,:))'))
% ga_freq.(exper.eventValues{1}).time(~isnan(squeeze(ga_freq.(exper.eventValues{1}).powspctrm(1,1,2,:))'))

%% decide who to kick out based on trial counts

numEv = struct;
numEv.thresh = 15;

% Subjects with bad behavior
exper.badBehSub = {};

% exclude subjects with low event counts
[exper,numEv] = mm_threshSubs(exper,ana,numEv,data_freq);

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = exper.eventValues;
cfg_ana.data_str = 'data_freq';
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sessions)
  for evVal = 1:length(exper.eventValues)
    %tic
    fprintf('Running ft_freqgrandaverage on %s...',exper.eventValues{evVal});
    ga_freq.(exper.eventValues{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(exper.eventValues{evVal}){ses}));
    fprintf('Done.\n');
    %toc
  end
end

%% save grand average file

if ~exist('cfg_proc','var')
  savedFiles = dir(fullfile(dirs.saveDir,'analysisDetails*.mat'));
  if length(savedFiles) == 1
    load(fullfile(dirs.saveDir,savedFiles.name),'cfg_proc');
  elseif length(savedFiles) > 1
    error('Multiple analysisDetails*.mat files found in %s!',dirs.saveDir)
  elseif isempty(savedFiles)
    error('analysisDetails*.mat not found in %s!',dirs.saveDir)
  end
end

saveFile = fullfile(dirs.saveDir,sprintf('ga_%s_%s%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(ga_freq.(exper.eventValues{1}).freq(1)),round(ga_freq.(exper.eventValues{1}).freq(end)),ga_freq.(exper.eventValues{1}).time(1)*1000,ga_freq.(exper.eventValues{1}).time(end)*1000));
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'ga_freq');
  fprintf('Done.\n');
else
  error('Not saving! %s already exists.\n',saveFile);
end

%% (re)save the analysis details

saveFile = fullfile(dirs.saveDir,sprintf('analysisDetails_%s_%s%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
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
  % http://www.amirwatad.com/blog/archives/2009/01/31/sending-emails-with-matlab/
  send_to = {'matt.mollison@gmail.com'};
  subject = sprintf('Done with%s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  message = {...
    sprintf('Done with%s %s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    sprintf('%s',saveFile),...
    };
  %attachments = {'picture1.png'};
  attachments = [];
  send_mail(send_to,subject,message,attachments);
end

%% if already saved and not yet loaded, load the ft_freqgrandaverage files

if ~exist('cfg_proc','var')
  savedFiles = dir(fullfile(dirs.saveDir,'analysisDetails*.mat'));
  if length(savedFiles) == 1
    load(fullfile(dirs.saveDir,savedFiles.name));
  elseif length(savedFiles) > 1
    error('Multiple analysisDetails*.mat files found in %s!',dirs.saveDir)
  elseif isempty(savedFiles)
    error('analysisDetails*.mat not found in %s!',dirs.saveDir)
  end
end

if ~exist('ga_freq','var')
  savedFiles = dir(fullfile(dirs.saveDir,sprintf('ga_%s_%s%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000)));
  %savedFiles = dir(fullfile(dirs.saveDir,sprintf('ga_freq_*.mat')));
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDir,savedFiles(sf).name));
    fprintf('Done.\n');
  end
end

%% plot the conditions - simple

cfg_ft = [];
%cfg_ft.baseline = [-0.3 -0.1];
%cfg_ft.baselinetype = 'absolute';
%if strcmp(cfg_ft.baselinetype,'absolute')
%cfg_ft.xlim = [0 1];
%cfg_ft.ylim = [3 50];
cfg_ft.ylim = [3 8];
%cfg_ft.ylim = [8 12];
%cfg_ft.zlim = [-1 1];
cfg_ft.zlim = [-100 100];
%elseif strcmp(cfg_ft.baselinetype,'relative')
%  cfg_ft.zlim = [0 2.0];
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
for evVal = 1:length(exper.eventValues)
  figure
  ft_multiplotTFR(cfg_ft,ga_freq.(exper.eventValues{evVal}));
  set(gcf,'Name',sprintf('%s',exper.eventValues{evVal}))
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
% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% % abbreviations for the condition types
% cfg_plot.typesByROI = {...
%   {'C2','C6'},...
%   {'C2','C6'}};

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.zlim = [-150 150];
cfg_ft.zparam = 'powspctrm';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.types = cfg_plot.typesByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotTFR(cfg_ft,cfg_plot,ana,exper,numEv,data_freq);
end

%% plot the conditions

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
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.zlim = [-150 150]; % pow

cfg_ft.zparam = 'powspctrm';

cfg_plot = [];
cfg_plot.plotTitle = 1;

%cfg_plot.rois = {{'FS'},{'LAS','RAS'},{'LPS','RPS'}};
%cfg_plot.rois = {{'FS'},{'PS'}};
%cfg_plot.rois = {'E71'};
cfg_plot.rois = {'all'};

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
%   };
% 
% cfg_plot.typesByROI = repmat({{'C2','C6'}},size(cfg_plot.condByTypeByROI));

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

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

% cfg_plot.ftFxn = 'ft_multiplotTFR';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.types = cfg_plot.typesByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_freq);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 1;

% comparisons to make
cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';

cfg_ft = [];
%cfg_ft.xlim = [.5 .8]; % time
cfg_ft.ylim = [3 8]; % freq
%cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.zparam = 'powspctrm';
cfg_ft.zlim = [-100 100]; % pow

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
cfg_ft.xlim = [0.5 0.8]; % time
%cfg_plot.subplot = 1;
%cfg_ft.xlim = [0 1.0]; % time
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
cfg_ana.frequencies = [3 8; 3 8];

%cfg_ana.conditions = {{'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},{'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'}};
cfg_ana.conditions = {{'all_within_types'}};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.avgoverfreq = 'yes';
cfg_ft.parameter = 'powspctrm';
cfg_ft.correctm = 'fdr';

% line plot parameters
cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
%cfg_plot.ylims = [-1 1; -1 1; -1 1];
cfg_plot.ylims = [-100 100; -100 100];
%cfg_plot.plot_order = {'CR2','H2','HSC2','HSI2','CR6','H6','HSC6','HSI6'};

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_freq);
end

%% 3-way ANOVA: Hemisphere x Block Type x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [3 8; 3 8];

% IV2: abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% IV3: outermost cell holds one cell for each ROI; each ROI cell holds one
% cell for each event type; each event type cell holds strings for its
% conditions
cfg_ana.condByTypeByROI = {...
  {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
  {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

% For each ROI, what's common among the conditions in each type
cfg_ana.condCommonByROI = {...
  {'CR','H','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Block Type','Condition'};

cfg_ana.parameter = 'powspctrm';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_ana.conditions = cfg_ana.condByTypeByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  
  mm_ft_rmaov33TFR(cfg_ana,exper,ana,data_tla);
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
  {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
  {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [3 8; 8 12];

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'powspctrm';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.frequency = cfg_ana.frequencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov2TFR(cfg_ana,exper,ana,data_freq);
end

%% cluster statistics

cfg_ft = [];
cfg_ft.avgoverchan = 'no';
cfg_ft.avgovertime = 'no';
cfg_ft.avgoverfreq = 'yes';

cfg_ft.parameter = 'powspctrm';

cfg_ana = [];
cfg_ana.roi = 'all';
cfg_ana.latencies = [0 1.0];
cfg_ana.conditions = {'all_within_types'};
% cfg_ana.conditions = {...
%   {'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},...
%   {'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'},...
%   {'CR2','CR6'},{'H2','H6'},{'HSC2','HSC6'},{'HSI2','HSI6'}};

cfg_ana.frequencies = [3 8; 8 12; 12 28; 28 40];
%cfg_ana.frequencies = [3 8; 8 12; 12 28; 28 50; 50 100];
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
cfg_ft.alpha = .1;

cfg_plot = [];
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.frequencies = cfg_ana.frequencies;
cfg_plot.latencies = cfg_ana.latencies;

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  for fr = 1:size(cfg_plot.frequencies,1)
    cfg_ft.frequency = cfg_plot.frequencies(fr,:);
    
    mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs);
  end
end

%% let me know that it's done
emailme = 1;
if emailme
  % http://www.amirwatad.com/blog/archives/2009/01/31/sending-emails-with-matlab/
  send_to = {'matt.mollison@gmail.com'};
  subject = sprintf('Done with %s pow:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  message = {...
    sprintf('Done with %s pow:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    };
  %attachments = {'picture1.png'};
  attachments = [];
  send_mail(send_to,subject,message,attachments);
end

%% correlations

cfg_ana = [];

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [3 8; 8 12];

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

cfg_ana.parameter = 'powspctrm';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.frequency = cfg_ana.frequencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeTFR(cfg_ana,ana,exper,files,dirs,data_tla);
end
