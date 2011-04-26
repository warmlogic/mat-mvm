% Make plots and do analyses for time-frequency EEG data

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SOSI';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
exper.eventValues = sort({'RCR','RHSC','RHSI'});

% combine the two types of hits into one category
exper.eventValuesExtra.toCombine = {{'RHSC','RHSI'}};
exper.eventValuesExtra.newValue = {{'RH'}};

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
exper.nsFileExt = 'egis';
%exper.nsFileExt = 'raw';

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
dirs.saveDirName = fullfile('ft_data',sprintf('pow_%seq%d',evStr,exper.equateTrials));

exper.subjects = {
  'SOSI001';
  'SOSI002';
  'SOSI003';
  'SOSI004';
  'SOSI005';
  'SOSI006';
  'SOSI007';
  'SOSI008';
  'SOSI009';
  'SOSI010';
  'SOSI011';
  'SOSI012';
  'SOSI013';
  'SOSI014';
  'SOSI015';
  'SOSI016';
  'SOSI017';
  'SOSI018';
  'SOSI020';
  'SOSI019';
  'SOSI021';
  'SOSI022';
  'SOSI023';
  'SOSI024';
  'SOSI025';
  'SOSI026';
  'SOSI027';
  'SOSI028';
  'SOSI029';
  'SOSI030';
  };
% original SOSI019 was replaced because the first didn't finish

% the sessions that each subject ran; multi-session support is not yet
% implemented, but for now this cell must contain one string; this will
% (probably/eventually) be the name of the directory containing the EEG
% files
exper.sessions = {'session_0'};

%% set up file and directory handling parameters

dirs.homeDir = getenv('HOME');

% directory where the data to read is located
dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));

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
dirs.saveDir = fullfile(dirs.dataroot,dirs.saveDirName);
if ~exist(dirs.saveDir,'dir')
  mkdir(dirs.saveDir)
end

% Assumes we have chan locs file in ~/Documents/MATLAB/mat_mvm/eeg/; use
% "short" because the Fiduciary points in the FT locs screws up some plots
files.elecfile = fullfile(dirs.homeDir,'Documents/MATLAB/mat_mvm/eeg/GSN_HydroCel_129_short.sfp');
%files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

% figure printing options
files.saveFigs = 0;
files.figFileExt = 'png';
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
dirs.saveDirFigs = fullfile(dirs.saveDir,'figs');
if ~exist(dirs.saveDirFigs,'dir')
  mkdir(dirs.saveDirFigs)
end

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

%% Convert the data to FieldTrip structs - excludes NS artifact trials
ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_freqanalysis';

% any preprocessing?
cfg_pp = [];

cfg_proc = [];
cfg_proc.pad = 'maxperlen';
cfg_proc.output = 'pow';
cfg_proc.keeptrials = 'no';

% wavelet
cfg_proc.method = 'wavelet';
cfg_proc.width = 5;
%cfg_proc.toi = -0.8:0.04:1.5;
cfg_proc.toi = -0.3:0.04:1.0;
% evenly spaced frequencies, but not as many as cfg_proc.foilim makes
freqstep = exper.sampleRate/(sum(abs(exper.prepost))*exper.sampleRate)*2;
%cfg_proc.foi = 3:freqstep:50;
%cfg_proc.foi = 3:freqstep:9;
cfg_proc.foi = 3:freqstep:100;

% % multi-taper method with hanning taper
% cfg_proc.method = 'mtmconvol';
% %cfg_proc.taper = 'hanning';
% cfg_proc.taper = 'dpss';
% %cfg_proc.toi = -0.5:0.04:1.5;
% cfg_proc.toi = -0.8:0.04:3.0;
% freqstep = exper.sampleRate/(sum(abs(exper.prepost))*exper.sampleRate)*2;
% cfg_proc.foi = 3:freqstep:50;
% %cfg_proc.foi = 3:freqstep:9;
% cfg_proc.t_ftimwin = 4./cfg_proc.foi;
% cfg_proc.tapsmofrq = 0.4*cfg_proc.foi;

[data_freq,exper] = create_ft_struct(ana,cfg_pp,cfg_proc,exper,dirs,files);

% [data_freq,exper] = create_ft_struct_peer('seg2ft','ft_freqanalysis',cfg_proc,cfg_pp,dirs,exper,files);

% extra string in filename for identification
if strcmp(cfg_proc.method,'wavelet')
  files.save_str = sprintf('_w%d',cfg_proc.width);
elseif strcmp(cfg_proc.method,'mtmconvol')
  files.save_str = sprintf('_%s',cfg_proc.taper);
elseif strcmp(cfg_proc.method,'mtmfft')
  files.save_str = sprintf('_%s',cfg_proc.taper);
else
  files.save_str = '';
end

% save the structs for loading in later
if strcmp(cfg_proc.keeptrials,'no')
  saveFile = fullfile(dirs.saveDir,sprintf('data_%s_%s%s_avg_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
elseif strcmp(cfg_proc.keeptrials,'yes')
  saveFile = fullfile(dirs.saveDir,sprintf('data_%s_%s%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,files.save_str,round(cfg_proc.foi(1)),round(cfg_proc.foi(end)),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
end
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'data_freq');
  fprintf('Done.\n');
else
  error('Not saving! %s already exists.',saveFile);
end

%% set up channel groups

ana = mm_ft_channelgroups(ana);

%% create fields for analysis and GA plotting; specific to each experiment

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
[exper,numEv] = mm_threshSubs(exper,numEv,data_freq);

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = exper.eventValues;
cfg_ana.data_str = 'data_freq';
cfg_ana.excludeBadSub = 1;
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
cfg_plot.condByROI = repmat({{'RCR','RH','RHSC','RHSI'}},size(cfg_plot.rois));

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.zlim = [-150 150];
cfg_ft.zparam = 'powspctrm';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotTFR(cfg_ft,cfg_plot,ana,exper,numEv,data_freq);
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
% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_plot.condByROI = repmat({{'RCR','RH','RHSC','RHSI'}},size(cfg_plot.rois));

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
cfg_ft.ylim = [8 12]; % freq
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
%cfg_ana.rois = {{'FS'},{'PS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [3 8; 3 8];

%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};
cfg_ana.conditions = {'all'};

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
cfg_ana.condByROI = {{'RCR','RH'},{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
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

cfg_ft.parameter = 'powspctrm';

cfg_ana = [];
cfg_ana.roi = 'all';
%cfg_ana.conditions = {'all'};
cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};

cfg_ana.frequencies = [3 8; 8 12; 12 28; 28 50; 50 100];
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
    
    mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs,stat_clus);
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

% d' values - check on these
cfg_ana.d_item = abs([]);
cfg_ana.d_source = abs([]);

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
% cfg_ft.zparam = 'stat';
% cfg_ft.zlim = [-5 5];
% ft_clusterplot(cfg_ft,stat_clus.RHSCvsRHSIvsRCR);
