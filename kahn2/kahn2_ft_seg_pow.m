% Connectivity analyses on time-frequency data

% See Maris & Oostenveld (2007) for a info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exp.name = 'KAHN2matt';

exp.sampleRate = 250;

% pre- and post-stimulus times to save, in seconds
exp.prepost = [-0.8 1.5];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
exper.eegFileExt = 'set';

% types of events to save; must be the same as the events in the NS files
exp.eventValues = sort({'CR__','CoTa','InTa'});
exp.eventNames = {'Correct Rejections','Correct Source','Incorrect Source'};
%exp.eventValues = sort({'CR__','CoTa','InTa','PP__','PR__','PTa_','RP__','RR__','RTa_'});
%exp.eventValues = sort({'CR__','PP__','PR__','RP__','RR__'});

exp.subjects = {
  'KAHN2 01';
  'KAHN2 03';
  'KAHN2 04';
  'KAHN2 05';
  'KAHN2 07';
  'KAHN2 08';
  'KAHN2 09';
  'KAHN2 10';
  'KAHN2 11';
  'KAHN2 12';
  'KAHN2 13';
  'KAHN2 14';
  'KAHN2 16';
  'KAHN2 17';
  'KAHN2 18';
  'KAHN2 19';
  'KAHN2 20';
  'KAHN2 21';
  'KAHN2 22';
  'KAHN2 23';
  'KAHN2 24';
  'KAHN2 25';
  'KAHN2 26';
  'KAHN2 27';
  'KAHN2 28';
  'KAHN2 29';
  'KAHN2 31';
  'KAHN2 32';
  'KAHN2 34';
  'KAHN2 38';
  'KAHN2 47';
  'KAHN2 94';
  };

% the exp.sessions that each subject ran
exp.sessions = {{'_1', '_2'}};

%% set up parameters

% directory where the data to read is located
dirs.dataDir = fullfile(exp.name,exp.name);

% directory to save the FT data; if undefined, set to dirs.dataDir
dirs.saveDirStem = fullfile('KAHN2matt','eeg','eegpp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));

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
files.figPrintFormat = 'png';
%files.figPrintFormat = 'epsc2';

%% Convert the data to FieldTrip structs

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_freqanalysis';

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
%freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
%cfg_proc.foi = 3:freqstep:50;
%cfg_proc.foi = 3:freqstep:9;
cfg_proc.foi = 3:1:9;
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
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% %cfg_proc.foi = 3:freqstep:50;
% cfg_proc.foi = 3:freqstep:9;
% %cfg_proc.foilim = [3 9];

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'pow');

ana.ftype = cfg_proc.output;

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,cfg_proc,exper,dirs,files);

%% save the analysis details

% overwrite if it already exists
saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails.mat'));
%if ~exist(saveFile,'file')
fprintf('Saving %s...',saveFile);
save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
fprintf('Done.\n');
%else
%  error('Not saving! %s already exists.\n',saveFile);
%end

%% load the analysis details

%[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(filename,1);

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

[data_freq] = mm_ft_loadSubjectData(exper,dirs,'pow');

%% fix the labels in the data becaues the EEGLAB data labels are different

for evVal = 1:length(exp.eventValues)
  for sub = 1:length(exp.subjects)
    for ses = 1:length(exp.sessions)
      data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data.label = data_freq.elec.label;
    end
  end
end

%% Change in freq relative to baseline using absolute power

cfg_fb = [];
cfg_fb.baseline = [-0.8 0];
cfg_fb.baselinetype = 'relative';

% and run ft_freqdescriptives to get the average???
cfg_fd = [];
cfg_fd.keeptrials = 'no';

%data_freq_orig = data_freq;

for evVal = 1:length(exp.eventValues)
  for sub = 1:length(exp.subjects)
    for ses = 1:length(exp.sessions)
      data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data);
      %data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg_fd,data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data);
    end
  end
end

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {};

% exclude subjects with low event counts
[exper] = mm_threshSubs(exper,ana,15);

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

saveFile = fullfile(dirs.saveDirProc,sprintf('ga_freq_%d_%d_%d_%d.mat',ga_freq.(exp.eventValues{1}).freq(1),ga_freq.(exp.eventValues{1}).freq(end),ga_freq.(exp.eventValues{1}).time(1)*1000,ga_freq.(exp.eventValues{1}).time(end)*1000));
if ~exist(saveFile,'file')
  save(saveFile,'ga_freq');
end

%% save the analysis details

saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails_freq_%d_%d_%d_%d.mat',ga_freq.(exp.eventValues{1}).freq(1),ga_freq.(exp.eventValues{1}).freq(end),ga_freq.(exp.eventValues{1}).time(1)*1000,ga_freq.(exp.eventValues{1}).time(end)*1000));
if ~exist(saveFile,'file')
  save(saveFile,'exp','dirs','files','numEv');
end

%% if already saved and not yet loaded, load the ft_freqgrandaverage files

if ~exist('ga_freq','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('ga_freq_*.mat')));
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDirProc,savedFiles(sf).name));
    fprintf('Done.\n');
  end
end

%% set up channel groups

ana.elecGroups = {...
  {'e32','e33','e38','e39','e43','e44','e128'},... % LAI
  {'e1','e114','e115','e120','e121','e122','e125'},... % RAI
  {'e12','e13','e19','e20','e24','e28','e29'},... % LAS
  {'e4','e5','e111','e112','e117','e118','e124'},... % RAS
  {'e37','e42','e52','e53','e54','e60','e61'},... % LPS
  {'e78','e79','e85','e86','e87','e92','e93'},... % RPS
  {'e57','e58','e63','e64','e65','e68','e69'},... %LPI
  {'e89','e90','e94','e95','e96','e99','e100'},... % RPI
  {'e4','e5','e6','e11','e12','e13','e19','e112'},... % Frontocentral
  {'e52','e53','e60','e61','e59','e66','e67'},... % LPS2
  {'e77','e78','e84','e85','e86','e91','e92'},... % RPS2
  {'e75'},... % P1 (central)
  {'e64'},... % N1 (lateral)
  };
ana.elecGroupsStr = {'LAI','RAI','LAS','RAS','LPS','RPS','LPI','RPI','FC','LPS2','RPS2','P1','N1'};

%% plot the conditions - simple

cfg_ft = [];
%cfg_ft.baseline = [-0.3 -0.1];
%cfg_ft.baselinetype = 'absolute';
%if strcmp(cfg_ft.baselinetype,'absolute')
%cfg_ft.zlim = [-400 400];
%elseif strcmp(cfg_ft.baselinetype,'relative')
cfg_ft.zlim = [0 2.0];
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ga_freq);
for evVal = 1:length(exp.eventValues)
  figure
  ft_multiplotTFR(cfg_ft,ga_freq.(exp.eventValues{evVal}));
  set(gcf,'Name',sprintf('%s',exp.eventValues{evVal}))
end

%% plot the conditions - Contrast

cfg_plot = [];
cfg_plot.plotTitle = 1;
cfg_plot.plotColorbar = 1;
if cfg_plot.plotColorbar
  cfg_plot.colorbar_str = '_cb';
else
  cfg_plot.colorbar_str = '';
end

% find all the pairwise event combinations for these contrasts
cfg_plot.comboEv = nchoosek(exp.eventValues,2);
% get the names for plotting
cfg_plot.comboName = nchoosek(exp.eventNames,2);

cfg_ft_m = [];
%cfg_ft_m.zlim = [-400 400];
cfg_ft_m.zlim = [0 1.5];
cfg_ft_m.showlabels = 'yes';
cfg_ft_m.colorbar = 'yes';
cfg_ft_m.interactive = 'yes';
cfg_ft_m.layout = ft_prepare_layout([],ga_freq);

cfg_ft_t = [];
%cfg_ft_t.xlim = [.6 1.2];
cfg_ft_t.xlim = [.5 .8];
%cfg_ft_t.xlim = [.3 .5];
cfg_ft_t.ylim = [4 8];
cfg_ft_t.marker = 'numbers';
cfg_ft_t.markerfontsize = 9;
cfg_ft_t.interactive = 'yes';
cfg_ft_t.colormap = 'hot';
cfg_ft_t.colorbar = 'yes';
cfg_ft_t.layout = ft_prepare_layout([],ga_freq);

% initialize for storing the contrast topoplots
cont_topo = [];

for c = 1:length(cfg_plot.comboEv)
  cfg_plot.cond = cfg_plot.comboEv(c,:);
  cfg_plot.condNames = cfg_plot.comboName(c,:);
  
  % create contrast
  cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})) = ga_freq.(cfg_plot.cond{1});
  if isfield(ga_freq.(cfg_plot.cond{1}),'powspctrm')
    cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})).powspctrm = ga_freq.(cfg_plot.cond{1}).powspctrm - ga_freq.(cfg_plot.cond{2}).powspctrm;
  end
  
  % make a multiplot
  figure
  ft_multiplotTFR(cfg_ft_m,cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})));
  set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
%   if cfg_plot.plotColorbar
%     h = colorbar;
%     set(get(h,'YLabel'),'string','Relative power');
%   end
  if cfg_plot.plotTitle
    title(sprintf('%s - %s',cfg_plot.condNames{1},cfg_plot.condNames{2}));
  end
  publishfig(gcf,0);
  
  % make a topoplot
  figure
  ft_topoplotTFR(cfg_ft_t,cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})));
%   if cfg_plot.plotColorbar
%     h = colorbar;
%     set(get(h,'YLabel'),'string','Relative power');
%   end
  if cfg_plot.plotTitle
    title(sprintf('%s - %s',cfg_plot.condNames{1},cfg_plot.condNames{2}));
  end
  publishfig(gcf,0);
  
  if files.saveFigs
    cfg_plot.figfilename = sprintf('topo_ga_freq_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
    print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
  end
end

%% power mean amplitude line plots (or bar plots)

