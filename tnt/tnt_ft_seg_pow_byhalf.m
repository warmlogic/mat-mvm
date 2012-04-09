% Make plots and do analyses for time-frequency EEG data

% See Maris & Oostenveld (2007) for a info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'TNT';

% types of events to save; these must be the same as the events in the NS
% files
%exper.eventValues = sort({'NT','TH'});
%exper.eventValues = sort({'B','NT','NTFor','NTRec','TFor','TH','TRec'});
exper.eventValues = sort({'B1of2','B2of2','NT1','NT2','TH1','TH2'});

% % combine the two types of hits into one category
% exper.eventValuesExtra.newValue = {{}};
% exper.eventValuesExtra.toCombine = {{}};

if ~isfield(exper,'eventValuesExtra')
  exper.eventValuesExtra = {};
end

% set the names manually to include on contrast plots
%exper.eventNames = {'NT','TH'};
%exper.eventNames = {'Baseline','No Think','NT Forgot','NT Recog','T Forgot','Think','T Recog'};
exper.eventNames = {'B1of2','B2of2','NT1','NT2','TH1','TH2'};

% pre- and post-stimulus times to save, in seconds
exper.prepost = [-0.8 3.0];

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw/sbin must be put in
% dirs.dataroot/ns_raw or egis must be put in dirs.dataroot/ns_egis
exper.eegFileExt = 'egis';
%exper.eegFileExt = 'raw';

% name of the folder to save the FT data in
dirs.saveDirName = 'ft_data';

exper.subjects = {
  'TNT 06';
  'TNT 07';
  'TNT 08';
  'TNT 09';
  'TNT 11';
  'TNT 13';
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
  'TNT 39';
  'TNT 41';
  'TNT 42';
  'TNT 44';
  'TNT 45';
  'TNT 46';
  'TNT 47';
  'TNT 48';
  'TNT 49';
  'TNT 50';
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
exper.sessions = {'session_0'};

%% set up parameters

dirs.homeDir = getenv('HOME');

% post processed
dirs.dataDir = sprintf('TNT_matt/eeg_byhalf/%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000);
%dirs.dataDir = sprintf('TNT_matt/eeg/%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000);
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

% Assumes we have chan locs file in ~/Documents/MATLAB/mat_mvm/eeg/
files.elecfile = fullfile(dirs.homeDir,'Documents/MATLAB/mat_mvm/eeg/GSN_HydroCel_129_short.sfp');
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
cfg_pp = [];

exper.equateTrials = 1;

cfg_proc = [];
cfg_proc.output = 'pow';
cfg_proc.pad = 'maxperlen';
cfg_proc.keeptrials = 'no';

% wavelet
cfg_proc.method = 'wavelet';
cfg_proc.width = 6;
%cfg_proc.toi = -0.5:0.04:1.5;
cfg_proc.toi = -0.8:0.04:3.0;
%cfg_proc.foi = 2:1:9;
cfg_proc.foi = 3:1:50;
%cfg_proc.foilim = [3 9];

% % multi-taper method with hanning taper
% cfg_proc.method = 'mtmconvol';
% cfg_proc.taper = 'hanning';
% %cfg_proc.toi = -0.5:0.04:1.5;
% cfg_proc.toi = -0.8:0.04:3.0;
% cfg_proc.foi = 2:1:50;
% %cfg_proc.foi = 3:1:9;
% %cfg_proc.foi = 4:1:8;
% cfg_proc.t_ftimwin = 5./cfg_proc.foi;
% cfg_proc.tapsmofrq = 0.4*cfg_proc.foi;

[data_freq,exper] = create_ft_struct('seg2ft','ft_freqanalysis',cfg_proc,cfg_pp,dirs,exper,files);

%[data_freq,exper.eventValues] = create_ft_struct_peer('ft_freqanalysis',cfg_proc,cfg_pp,dirs.dataroot,dirs.saveDirName,exper.eegFileExt,exper.subjects,exper.sessions,exper.prepost,files.elecfile,exper.sampleRate,exper.eventValues,exper.eventValuesExtra);
%[data_freq,exper.eventValues] = create_ft_struct_peer('ft_freqanalysis',cfg_proc,cfg_pp,dirs,exper,files);

% save the structs for loading in later
if strcmp(cfg_proc.keeptrials,'no')
  saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s_%s_avg_eq%d_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,exper.equateTrials,cfg_proc.foi(1),cfg_proc.foi(end),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
elseif strcmp(cfg_proc.keeptrials,'yes')
  saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s_%s_eq%d_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,exper.equateTrials,cfg_proc.foi(1),cfg_proc.foi(end),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000));
end
if ~exist(saveFile,'file')
  fprintf('Saving %s\n',saveFile);
  save(saveFile,'data_freq','cfg_proc');
end

%% if already saved and not yet loaded, load the ft_freqanalysis files

if ~exist('data_freq','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('data_%s_%s_avg_eq%d_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.method,exper.equateTrials,cfg_proc.foi(1),cfg_proc.foi(end),cfg_proc.toi(1)*1000,cfg_proc.toi(end)*1000)));
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDirProc,savedFiles(sf).name));
    fprintf('Done.\n');
  end
  % get all the exper.eventValues and exper.eventValuesExtra together; make sure the extra event values aren't in the list
  if ~isempty(exper.eventValuesExtra)
    for nVal = 1:length(exper.eventValuesExtra.newValue)
      if ~ismember(exper.eventValuesExtra.newValue{nVal},exper.eventValues)
        exper.eventValues = cat(2,exper.eventValues,exper.eventValuesExtra.newValue{nVal});
      else
        fprintf('%s is already in the event value list!\n',exper.eventValuesExtra.newValue{nVal}{1});
      end
    end
    exper.eventValues = sort(exper.eventValues);
  end
end

%% Potential bad channels




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
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],data_freq);
figure
ft_multiplotTFR(cfg_ft,data_freq.(exper.eventValues{1}).sub(1).ses(1).data);
% 
% cfg_ft = [];
% cfg_ft.channel = {'e124'};
% %cfg_ft.channel = {'e117'};
% cfg_ft.baseline = [-0.3 -0.1];
% cfg_ft.baselinetype = 'absolute';
% if strcmp(cfg_ft.baselinetype,'absolute')
%   %cfg_ft.zlim = [-2 2];
%   cfg_ft.zlim = [-500 500];
% elseif strcmp(cfg_ft.baselinetype,'relative')
%   cfg_ft.zlim = [0 1.5];
% end
% cfg_ft.showlabels = 'yes';
% cfg_ft.colorbar = 'yes';
% cfg_ft.ylim = [4 8];
% figure
% ft_singleplotTFR(cfg_ft,data_freq.(exper.eventValues{1}).sub(1).ses(1).data);
% 
%figure
% badTrials = [];
% for ev = 1:size(data_freq.(exper.eventValues{1}).sub(1).ses(1).data.powspctrm,1)
%   if data_freq.(exper.eventValues{1}).sub(1).ses(1).data.powspctrm(ev,117,5,20) < -800
%     badTrials = [badTrials ev];
% %     fprintf('%d: Low! %.1f\n',ev,data_freq.(exper.eventValues{1}).sub(1).ses(1).data.powspctrm(ev,117,5,20));
% %     clf
% %     imagesc(squeeze(data_freq.(exper.eventValues{1}).sub(1).ses(1).data.powspctrm(ev,117,:,:)),[-700 700]);
% %     axis xy
% %     title(['Trial ',num2str(ev)]);
% %     xlabel('Time (ms)');
% %     ylabel('Frequency (Hz)');
% %     set(gca,'XTick',1:5:length(data_freq.(exper.eventValues{1}).sub(1).ses(1).data.time));
% %     set(gca,'XTickLabel',data_freq.(exper.eventValues{1}).sub(1).ses(1).data.time(1):.25:data_freq.(exper.eventValues{1}).sub(1).ses(1).data.time(end));
% %     set(gca,'YTickLabel',data_freq.(exper.eventValues{1}).sub(1).ses(1).data.freq);
% %     h = colorbar;
% %     set(get(h,'YLabel'),'string','Power');
% %     keyboard
% %   else
% %     fprintf('%d: %.1f\n',ev,data_freq.(exper.eventValues{1}).sub(1).ses(1).data.powspctrm(ev,117,5,20));
% %     continue
%   end
% end
% 
% cfg_bcr = [];
% cfg_bcr.badchannel = {'e117'};
% cfg_bcr.trials = badTrials;
% [interp] = ft_channelrepair(cfg_bcr,data_freq.(exper.eventValues{1}).sub(1).ses(1).data);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cfg_ft = [];
% cfg_ft.baseline = [-0.5 -0.1];	
% cfg_ft.baselinetype = 'absolute';
% if strcmp(cfg_ft.baselinetype,'absolute')
%   cfg_ft.zlim = [-2 2];
% elseif strcmp(cfg_ft.baselinetype,'relative')
%   cfg_ft.zlim = [0 1.5];
% end
% cfg_ft.xlim = [.5 .8];
% cfg_ft.ylim = [4 8];
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.interactive = 'yes';
% cfg_ft.colormap = 'hot';
% cfg_ft.colorbar = 'yes';
% cfg_ft.layout = ft_prepare_layout([],data_freq);
% figure
% ft_topoplotTFR(cfg_ft,data_freq.(exper.eventValues{1}).sub(1).ses(1).data);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast
% 
% data_freq_RHSCvsRHSI = data_freq.RHSC.sub(1).ses(1).data;
% data_freq_RHSCvsRHSI.powspctrm = data_freq.RHSC.sub(1).ses(1).data.powspctrm - data_freq.RHSI.sub(1).ses(1).data.powspctrm;
% 
% cfg_ft = [];
% cfg_ft.zlim = [-400 400];
% %cfg_ft.zlim = [0 2.0];
% cfg_ft.showlabels = 'yes';
% cfg_ft.colorbar = 'yes';
% cfg_ft.interactive = 'yes';
% cfg_ft.layout = ft_prepare_layout([],data_freq);
% figure
% ft_multiplotTFR(cfg_ft,data_freq_RHSCvsRHSI);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Contrast
% 
% cfg_ft = [];
% cfg_ft.xlim = [.5 .8];
% %cfg_ft.xlim = [.3 .5];
% cfg_ft.ylim = [4 8];
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.interactive = 'yes';
% cfg_ft.colormap = 'hot';
% cfg_ft.colorbar = 'yes';
% cfg_ft.layout = ft_prepare_layout([],data_freq);
% figure
% ft_topoplotTFR(cfg_ft,data_freq_RHSCvsRHSI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Change in freq relative to baseline using absolute power

%data_freq_old = data_freq;

cfg_fb = [];
cfg_fb.baseline = [-0.3 -0.1];
cfg_fb.baselinetype = 'absolute';

% and run ft_freqdescriptives to get the average
cfg_fd = [];
cfg_fd.keeptrials = 'no';

data_freq_orig = data_freq;

for evVal = 1:length(exper.eventValues)
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      data_freq.(exper.eventValues{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_freq.(exper.eventValues{evVal}).sub(sub).ses(ses).data);
      data_freq.(exper.eventValues{evVal}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg_fd,data_freq.(exper.eventValues{evVal}).sub(sub).ses(ses).data);
    end
  end
end

%% decide who to kick out based on trial counts

numEv = struct;
numEv.thresh = 15;

% Subjects with bad behavior
exper.badBehSub = {};

% exclude subjects with low event counts
[exper,numEv] = mm_threshSubs(exper,numEv,data_freq);

%% set up ana.ga_str strings to put in ft_freqgrandaverage

[ana] = mm_ft_ga_str(exper,ana,'data_freq');

%% get the grand average and add the electrode locations again

cfg_ft = [];
%cfg_ft.keeptrials = 'yes';
cfg_ft.keepindividual = 'yes';
for ses = 1:length(exper.sessions)
  for evVal = 1:length(exper.eventValues)
    tic
    fprintf('Running ft_freqgrandaverage on %s...',exper.eventValues{evVal});
    ga_freq.(exper.eventValues{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',ana.ga_str.(exper.eventValues{evVal}){ses}));
    fprintf('Done.\n');
    toc
  end
end
ga_freq.elec = data_freq.elec;

%% save grand average file

saveFile = fullfile(dirs.saveDirProc,sprintf('ga_freq_%d_%d_%d_%d.mat',round(ga_freq.(exper.eventValues{1}).freq(1)),round(ga_freq.(exper.eventValues{1}).freq(end)),ga_freq.(exper.eventValues{1}).time(1)*1000,ga_freq.(exper.eventValues{1}).time(end)*1000));
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'ga_freq');
  fprintf('Done.\n');
end

%% save the analysis details

saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails_freq_%d_%d_%d_%d.mat',round(ga_freq.(exper.eventValues{1}).freq(1)),round(ga_freq.(exper.eventValues{1}).freq(end)),ga_freq.(exper.eventValues{1}).time(1)*1000,ga_freq.(exper.eventValues{1}).time(end)*1000));
if ~exist(saveFile,'file')
  fprintf('Saving %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','numEv');
  fprintf('Done.\n');
end

%% if already saved and not yet loaded, load the ft_freqgrandaverage files

if ~exist('ga_freq','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('ga_freq_%d_%d_%d_%d.mat',ga_freq.(exper.eventValues{1}).freq(1),ga_freq.(exper.eventValues{1}).freq(end),ga_freq.(exper.eventValues{1}).time(1)*1000,ga_freq.(exper.eventValues{1}).time(end)*1000)));
  %savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('ga_freq_*.mat')));
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
  {'e4','e5','e6','e11','e12','e13','e19','e112'},... % Mid-frontal
  {'e52','e53','e60','e61','e59','e66','e67'},... % LPS2
  {'e77','e78','e84','e85','e86','e91','e92'},... % RPS2
  {'e62','e67','e71','e72','e75','e76','e77'},... % PS
  {'e75'},... % P1 (central)
  {'e64'},... % N1 (lateral)
  };
ana.elecGroupsStr = {'LAI','RAI','LAS','RAS','LPS','RPS','LPI','RPI','MF','LPS2','RPS2','PS','P1','N1'};

%% create fields for plotting and analysis; specific to each experiment

% take out any existing fields
if isfield(ana.events,'types')
  ana.events = rmfield(ana.events,'types');
end
if isfield(ana.events,'values')
  ana.events = rmfield(ana.events,'values');
end
if isfield(ana.events,'names')
  ana.events = rmfield(ana.events,'names');
end
if isfield(ana.events,'common')
  ana.events = rmfield(ana.events,'common');
end
%if isfield(ana.ga_str,'ga_plot')
%  ana.events = rmfield(ana.ga_str,'ga_plot');
%end

% divided across session half
ana.events.types{1} = 'First_Half';
ana.events.types{2} = 'Second_Half';
ana.events.values{1} = exper.eventValues(1);
ana.events.values{2} = exper.eventValues(2);
ana.events.names{1} = exper.eventNames(1);
ana.events.names{2} = exper.eventNames(2);
%ana.ga_str.ga_plot{1} = sprintf('ga_freq.%s',exper.eventValues{1});
%ana.ga_str.ga_plot{2} = sprintf('ga_freq.%s',exper.eventValues{2});
for evVal = 1:length(exper.eventValues)
  if ~ismember(exper.eventValues{evVal},ana.events.values{1}) && strcmp(exper.eventValues{evVal},'B1of2') || strcmp(exper.eventValues{evVal},'NT1') || strcmp(exper.eventValues{evVal},'TH1')
    ana.events.values{1} = cat(2,ana.events.values{1},exper.eventValues{evVal});
    ana.events.names{1} = cat(2,ana.events.names{1},exper.eventNames{evVal});
    %ana.ga_str.ga_plot{1} = cat(2,ana.ga_str.ga_plot{1},sprintf(',ga_freq.%s',exper.eventValues{evVal}));
  elseif ~ismember(exper.eventValues{evVal},ana.events.values{2}) && strcmp(exper.eventValues{evVal},'B2of2') || strcmp(exper.eventValues{evVal},'NT2') || strcmp(exper.eventValues{evVal},'TH2')
    ana.events.values{2} = cat(2,ana.events.values{2},exper.eventValues{evVal});
    ana.events.names{2} = cat(2,ana.events.names{2},exper.eventNames{evVal});
    %ana.ga_str.ga_plot{2} = cat(2,ana.ga_str.ga_plot{2},sprintf(',ga_freq.%s',exper.eventValues{evVal}));
  end
end
% what's common between the two types of events
if length(ana.events.types) == 1
  ana.events.common = exper.eventValues;
else
  ana.events.common = {'B','NT','TH'};
end

% % divided across event types
% ana.events.types{1} = 'B';
% ana.events.types{2} = 'NT';
% ana.events.types{3} = 'TH';
% ana.events.values{1} = exper.eventValues(1);
% ana.events.values{2} = exper.eventValues(3);
% ana.events.values{3} = exper.eventValues(5);
% ana.events.names{1} = exper.eventNames(1);
% ana.events.names{2} = exper.eventNames(3);
% ana.events.names{3} = exper.eventNames(5);
% %ana.ga_str.ga_plot{1} = sprintf('ga_freq.%s',exper.eventValues{1});
% %ana.ga_str.ga_plot{2} = sprintf('ga_freq.%s',exper.eventValues{3});
% %ana.ga_str.ga_plot{3} = sprintf('ga_freq.%s',exper.eventValues{5});
% for evVal = 1:length(exper.eventValues)
%   if strcmp(exper.eventValues{evVal},'B2of2')
%     ana.events.values{1} = cat(2,ana.events.values{1},exper.eventValues{evVal});
%     ana.events.names{1} = cat(2,ana.events.names{1},exper.eventNames{evVal});
%     %ana.ga_str.ga_plot{1} = cat(2,ana.ga_str.ga_plot{1},sprintf(',ga_freq.%s',exper.eventValues{evVal}));
%   elseif strcmp(exper.eventValues{evVal},'NT2')
%     ana.events.values{2} = cat(2,ana.events.values{2},exper.eventValues{evVal});
%     ana.events.names{2} = cat(2,ana.events.names{2},exper.eventNames{evVal});
%     %ana.ga_str.ga_plot{2} = cat(2,ana.ga_str.ga_plot{2},sprintf(',ga_freq.%s',exper.eventValues{evVal}));
%   elseif strcmp(exper.eventValues{evVal},'TH2')
%     ana.events.values{3} = cat(2,ana.events.values{3},exper.eventValues{evVal});
%     ana.events.names{3} = cat(2,ana.events.names{3},exper.eventNames{evVal});
%     %ana.ga_str.ga_plot{3} = cat(2,ana.ga_str.ga_plot{3},sprintf(',ga_freq.%s',exper.eventValues{evVal}));
%   end
% end
% % what's common between the two types of events
% if length(ana.events.types) == 1
%   ana.events.common = exper.eventValues;
% else
%   ana.events.common = {'First_Half','Second_Half'};
% end

%% plot the conditions - simple

cfg_ft = [];
%cfg_ft.baseline = [-0.3 -0.1];
%cfg_ft.baselinetype = 'absolute';
%if strcmp(cfg_ft.baselinetype,'absolute')
%cfg_ft.xlim = [0 1];
cfg_ft.ylim = [3 50];
%cfg_ft.ylim = [3 8];
%cfg_ft.zlim = [-1 1];
cfg_ft.zlim = [-100 100];
%elseif strcmp(cfg_ft.baselinetype,'relative')
%  cfg_ft.zlim = [0 2.0];
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ga_freq);
for evVal = 1:length(exper.eventValues)
  figure
  ft_multiplotTFR(cfg_ft,ga_freq.(exper.eventValues{evVal}));
  set(gcf,'Name',sprintf('%s',exper.eventValues{evVal}))
end

%% subplots of each subject's power spectrum

cfg_plot = [];
%cfg_plot.roi = {'LAS','RAS'};
%cfg_plot.roi = {'RAS'};
%cfg_plot.roi = {'LPS','RPS'};
cfg_plot.roi = {'PS'};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
%cfg_plot.zlim = 'maxmin';
%cfg_plot.zlim = [-1 1];
cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.zlim = [-100 100];

mm_ft_subjplotTFR(cfg_ft,cfg_plot,exper,ana,numEv,data_freq);

%% plot the GA topoplots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [.5 1.0]; % time
%cfg_ft.ylim = [3 8]; % freq
cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.zlim = [-100 100]; % pow

cfg_plot = [];

%cfg_plot.ftFxn = 'ft_singleplotTFR';
% %cfg_plot.roi = {'e20'};
%cfg_plot.roi = {'MF'};

cfg_plot.ftFxn = 'ft_topoplotTFR';
cfg_plot.subplot = 1;
cfg_ft.xlim = [0 1.0]; % time
%cfg_plot.roi = {'MF'};
% %cfg_plot.roi = {'LAS','RAS'};
% %cfg_plot.roi = {'RAS'};
% %cfg_plot.roi = {'LPS','RPS'};
% %cfg_plot.roi = {'LPS'};

%cfg_plot.ftFxn = 'ft_multiplotTFR';
% cfg_plot.roi = {'MF'};
% cfg_plot.roi = {'PS'};
% cfg_plot.roi = {'LAS','RAS'};
% cfg_plot.roi = {'RAS'};
% cfg_plot.roi = {'LPS','RPS'};
% cfg_plot.roi = {'LPS'};

mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_freq);

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 1;

cfg_ft = [];
%cfg_ft.xlim = [.5 .8]; % time
cfg_ft.ylim = [3 8]; % freq
%cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.zlim = [-100 100]; % pow
cfg_ft.interactive = 'yes';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'yes';
cfg_ft.layout = ft_prepare_layout([],ga_freq);

cfg_plot.ftFxn = 'ft_topoplotTFR';
cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
%cfg_ft.comment = 'no';
%cfg_ft.xlim = [0 1.0]; % time
cfg_ft.xlim = (0:0.1:1.0); % time
cfg_plot.subplot = 1;

% cfg_plot.ftFxn = 'ft_multiplotTFR';
% cfg_ft.showlabels = 'yes';

cfg_plot.comboEvs = {'all'};

mm_ft_contrastTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_freq);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'},{'MF'},{'PS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8; 0.5 0.8; .5 1.0];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [5 8; 5 8; 5 8; 8 14];

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.avgoverfreq = 'yes';

cfg_ana.individ_plots = 0;
cfg_ana.line_plots = 0;
% line plot parameters
cfg_plot = [];
%cfg_plot.ylims = [-1 1; -1 1; -1 1];
cfg_plot.ylims = [-100 100; -100 100; -100 100; -100 100];

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,ga_freq);
end

%% 3-way ANOVA: Hemisphere x Session Half x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [5 8; 8 14];

% IV2: ana.events.types

% IV3: ana.events.values{typ}

cfg_ana.IV_names = {'ROI','Session Half','Condition'};
%cfg_ana.IV_names = {'ROI','Condition','Session Half'};

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.frequency = cfg_ana.frequencies(r,:);
  
  mm_ft_rmaov33TFR(cfg_ana,ana,ga_freq);
end

%% 2-way ANOVA: Hemisphere x Condition

% will be run once for each 

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'},{'MF','PS'}};
% IV2: define the conditions tested for each set of ROIs; one cell for each
% ROI; each cells must have one cell for each ana.events.types; don't
% include it if you just want to use ana.events.values for each ROI
%cfg_ana.conditions = {ana.events.values,ana.events.values,ana.events.values};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8; 0.5 0.9];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [4 8; 4 8; 8 14];

cfg_ana.IV_names = {'ROI','Condition'};

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  %cfg_ana.cond = cfg_ana.conditions{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.frequency = cfg_ana.frequencies(r,:);
  
  mm_ft_rmaov2TFR(cfg_ana,ana,ga_freq);
end

%% cluster statistics

cfg_ft = [];
%start and stop time of analysis, in sec
%cfg_ft.latency = [1.0 2.0];
cfg_ft.avgovertime = 'no';

%cfg_ft.frequency = [4 8]; % none
%cfg_ft.frequency = [8 14]; % posterior, 600-750: NT > TH (p=.06)
%cfg_ft.frequency = [14 30]; % posterior, 600-750: NT > TH (p=.081)
%cfg_ft.frequency = [33 50]; % .1-3 posterior, 640-880: NT > TH (pos) p = .054
cfg_ft.avgoverfreq = 'yes';

cfg_ft.elec = ga_freq.elec;
cfg_ft.channel = 'all';
%cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS'})});
%cfg_ft.avgoverchan = 'yes';

cfg_ana = [];
cfg_ana.comboEvs = {'all'};

%[stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,ga_freq);

cfg_ana.frequencies = [5 8; 8 14; 14 32; 32 50];
cfg_ana.latencies = [0 1; 1 2; 2 3];

for fr = 1:size(cfg_ana.frequencies,1)
  cfg_ft.frequency = cfg_ana.frequencies(fr,:);
  for lat = 1:size(cfg_ana.latencies,1)
    cfg_ft.latency = cfg_ana.latencies(lat,:);
    
    [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,ga_freq,dirs);
    
    %save(fullfile(dirs.saveDirProc,sprintf('freq_stat_clus_%d_%d_%d_%d.mat',cfg_ft.frequency(1),cfg_ft.frequency(2),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000)),'stat_clus');
  end
end

%load(fullfile(dirs.saveDirProc,sprintf('freq_stat_clus_%d_%d_%d_%d.mat',cfg_ft.frequency(1),cfg_ft.frequency(2),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000)));

% % run the nonparametric cluster statistics
% stat_clus.(sprintf('%svs%svs%s',cfg_ana.cond{1},cfg_ana.cond{2},cfg_ana.cond{3})) = ft_freqstatistics(cfg_ft,ga_tla.(cfg_ana.cond{1}),ga_tla.(cfg_ana.cond{2}),ga_tla.(cfg_ana.cond{3}));

%% plot the cluster statistics

cfg_ft = [];
cfg_ft.alpha = .1;

cfg_plot = [];
cfg_plot.comboEvs = {'all'};

for fr = 1:size(cfg_ana.frequencies,1)
  cfg_ft.frequency = cfg_ana.frequencies(fr,:);
  for lat = 1:size(cfg_ana.latencies,1)
    cfg_ft.latency = cfg_ana.latencies(lat,:);
    
    savedfile = fullfile(dirs.saveDirProc,sprintf('freq_stat_clus_%d_%d_%d_%d.mat',cfg_ft.frequency(1),cfg_ft.frequency(2),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000));
    if exist(savedfile,'file')
      fprintf('%d--%d Hz, %d--%d s\n',cfg_ft.frequency(1),cfg_ft.frequency(2),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000);
      load(savedfile);
    end
    
    mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs,stat_clus,ga_freq);
  end
end

%% Make contrast plots (with culster stat info)

% % set up contrast
% cfg_ana = [];
% cfg_ana.include_clus_stat = 0;
% cfg_ana.timeS = (0:0.05:1.0);
% cfg_ana.timeSamp = round(linspace(1,exper.sampleRate,length(cfg_ana.timeS)));
% 
% cfg_plot = [];
% cfg_plot.zlim = [-1 1];
% cfg_plot.numRows = 4;
% 
% cfg_ft = [];
% cfg_ft.interactive = 'no';
% cfg_ft.elec = ga_freq.elec;
% cfg_ft.highlight = 'on';
% if cfg_ana.include_clus_stat == 0
%   cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS','RPS','LPS'})});
% elseif cfg_ana.include_clus_stat == 1
%   cfg_ft.highlightcolorpos = [0.5 0 1];
%   cfg_ft.highlightcolorneg = [0 0.5 0];
% end
% cfg_ft.comment = 'no';
% %cfg_ft.commentpos = 'title';
% 
% % create contrast
% cont_topo = [];
% cont_topo.RCRvsRH = ga_freq.RH;
% cont_topo.RCRvsRH.powspctrm = ga_freq.RH.powspctrm - ga_freq.RCR.powspctrm;
% if cfg_ana.include_clus_stat == 1
%   pos = stat_clus.RCRvsRH.posclusterslabelmat==1;
% end
% % make a plot
% figure
% for k = 1:length(cfg_ana.timeS)-1
%   subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
%   cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
%   cfg_ft.zlim = [cfg_plot.zlim(1) cfg_plot.zlim(2)];
%   if cfg_ana.include_clus_stat == 1
%     pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
%     cfg_ft.highlightchannel = find(pos_int==1);
%   end
%   ft_topoplotTFR(cfg_ft,cont_topo.RCRvsRH);
% end
% set(gcf,'Name','CR - H')
% 
% % create contrast
% cont_topo.RHSCvsRHSI = ga_freq.RHSC;
% cont_topo.RHSCvsRHSI.powspctrm = ga_freq.RHSC.powspctrm - ga_freq.RHSI.powspctrm;
% if cfg_ana.include_clus_stat == 1
%   pos = stat_clus.RHSCvsRHSI.posclusterslabelmat==1;
% end
% % make a plot
% figure
% for k = 1:length(cfg_ana.timeS)-1
%   subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
%   cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
%   cfg_ft.zlim = [cfg_plot.zlim(1) cfg_plot.zlim(2)];
%   if cfg_ana.include_clus_stat == 1
%     pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
%     cfg_ft.highlightchannel = find(pos_int==1);
%   end
%   ft_topoplotER(cfg_ft,cont_topo.RHSCvsRHSI);
% end
% set(gcf,'Name','HSC - HSI')
% 
% % create contrast
% cont_topo.RHSCvsRCR = ga_freq.RHSC;
% cont_topo.RHSCvsRCR.powspctrm = ga_freq.RHSC.powspctrm - ga_freq.RCR.powspctrm;
% if cfg_ana.include_clus_stat == 1
%   pos = stat_clus.RHSCvsRCR.posclusterslabelmat==1;
% end
% % make a plot
% figure
% for k = 1:length(cfg_ana.timeS)-1
%   subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
%   cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
%   cfg_ft.zlim = [cfg_plot.zlim(1) cfg_plot.zlim(2)];
%   if cfg_ana.include_clus_stat == 1
%     pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
%     cfg_ft.highlightchannel = find(pos_int==1);
%   end
%   ft_topoplotER(cfg_ft,cont_topo.RHSCvsRCR);
% end
% set(gcf,'Name','HSC - CR')
% 
% % % mecklinger plot
% % cfg_ft.xlim = [1200 1800];
% % cfg_ft.zlim = [-2 2];
% % cfg_ft.colorbar = 'yes';
% % figure
% % ft_topoplotER(cfg_ft,cont_topo.RHSCvsRCR);
% 
% % create contrast
% cont_topo.RHSIvsRCR = ga_freq.RHSI;
% cont_topo.RHSIvsRCR.avg = ga_freq.RHSI.avg - ga_freq.RCR.avg;
% cont_topo.RHSIvsRCR.powspctrm = ga_freq.RHSI.powspctrm - ga_freq.RCR.powspctrm;
% if cfg_ana.include_clus_stat == 1
%   pos = stat_clus.RHSIvsRCR.posclusterslabelmat==1;
% end
% % make a plot
% figure
% for k = 1:length(cfg_ana.timeS)-1
%   subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
%   cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
%   cfg_ft.zlim = [cfg_plot.zlim(1) cfg_plot.zlim(2)];
%   if cfg_ana.include_clus_stat == 1
%     pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
%     cfg_ft.highlightchannel = find(pos_int==1);
%   end
%   ft_topoplotER(cfg_ft,cont_topo.RHSIvsRCR);
% end
% set(gcf,'Name','HSI - CR')
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % % because of a bug (might be fixed now)
% % if ~isfield(stat_clus.RHSCvsRHSIvsRCR,'negclusters') && isfield(stat_clus.RHSCvsRHSIvsRCR,'posclusters')
% %   fprintf('No neg clusters found\n');
% %   stat_clus.RHSCvsRHSIvsRCR.negclusters.prob = .5;
% %   stat_clus.RHSCvsRHSIvsRCR.negclusters.clusterstat = 0;
% %   stat_clus.RHSCvsRHSIvsRCR.negclusterslabelmat = zeros(size(stat_clus.RHSCvsRHSIvsRCR.posclusterslabelmat));
% %   stat_clus.RHSCvsRHSIvsRCR.negdistribution = zeros(size(stat_clus.RHSCvsRHSIvsRCR.posdistribution));
% % end
% % if ~isfield(stat_clus.RHSCvsRHSIvsRCR,'posclusters') && isfield(stat_clus.RHSCvsRHSIvsRCR,'negclusters')
% %   fprintf('No pos clusters found\n');
% %   stat_clus.RHSCvsRHSIvsRCR.posclusters.prob = 1;
% %   stat_clus.RHSCvsRHSIvsRCR.posclusters.clusterstat = 0;
% %   stat_clus.RHSCvsRHSIvsRCR.posclusterslabelmat = zeros(size(stat_clus.RHSCvsRHSIvsRCR.negclusterslabelmat));
% %   stat_clus.RHSCvsRHSIvsRCR.posdistribution = zeros(size(stat_clus.RHSCvsRHSIvsRCR.negdistribution));
% % end
% %
% % cfg_ft = [];
% % % p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
% % cfg_ft.highlightsymbolseries = ['*','*','.','.','.'];
% % cfg_ft.layout = ft_prepare_layout(cfg_ft,ga_freq);
% % cfg_ft.contournum = 0;
% % cfg_ft.emarker = '.';
% % cfg_ft.alpha  = 0.05;
% % cfg_ft.parameter = 'stat';
% % cfg_ft.zlim = [-5 5];
% % ft_clusterplot(cfg_ft,stat_clus.RHSCvsRHSIvsRCR);
