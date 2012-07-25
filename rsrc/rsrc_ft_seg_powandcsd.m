% Do a nonparametric statistics clustering analysis for timelocked EEG
% (ERPs)

% See Maris and Oostenveld (2007) for a description of statistics

%% Experiment-specific setup

exp.name = 'RSRC';

% types of events to save; these must be the same as the events in the NS
% files
exp.eventValues = sort({'RCR','RHSC','RHSI'});
%exp.eventValues = sort({'RCR','RH','RHSC','RHSI'});

% combine the two types of hits into one category
exp.eventValuesExtra.newValue = {{'RH'}};
exp.eventValuesExtra.toCombine = {{'RHSC','RHSI'}};

% if ~isfield(exp,'eventValuesExtra')
%   exp.eventValuesExtra = {};
% end

% pre- and post-stimulus times to save, in seconds
exp.prepost = [-1.0 2.0];

exp.sampleRate = 250;

% type of NS file for FieldTrip to read; raw/sbin must be put in
% dirs.dataroot/ns_raw or egis must be put in dirs.dataroot/ns_egis
exp.eegFileExt = 'egis';
%exp.eegFileExt = 'raw';

% name of the folder to save the FT data in
dirs.saveDirName = 'ft_data';

exp.subjects = {
  'RSRC021';
  'RSRC022';
  'RSRC023';
  'RSRC024';
  'RSRC025';
  'RSRC026';
  'RSRC027';
  'RSRC028';
  'RSRC029';
  'RSRC030';
  'RSRC031';
  'RSRC032';
  'RSRC033';
  'RSRC034';
  'RSRC035';
  'RSRC036';
  'RSRC037';
  'RSRC038';
  'RSRC039';
  'RSRC040';
  'RSRC041';
  'RSRC042';
  'RSRC043';
  'RSRC044';
  'RSRC045';
  'RSRC046';
  'RSRC048';
  'RSRC049';
  'RSRC050';
  'RSRC051';
  'RSRC052';
  'RSRC053';
  'RSRC054';
  'RSRC055';
  'RSRC056';
  'RSRC057';
  'RSRC058';
  'RSRC059';
  };
%    'RSRC047'; % bad subject; falling asleep

% the exp.sessions that each subject ran
exp.sessions = {'session_0'};

%% set up parameters

dirs.homeDir = getenv('HOME');

% EP post processed
dirs.dataDir = sprintf('eeg/eppp_nobackup/%d_%d',exp.prepost(1)*1000,exp.prepost(2)*1000);
dirs.serverDir = fullfile('/Volumes/curranlab/Data',exp.name,dirs.dataDir);
dirs.serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data',exp.name,dirs.dataDir);

% pick the right dirs.dataroot
if exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
elseif exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
else
  error('Server data drive not found.');
  %dirs.dataroot = fullfile(dirs.homeDir,'data',exp.name,dirs.dataDir);
end

% directory to save the data
dirs.saveDirProc = fullfile(dirs.dataroot,dirs.saveDirName);
if ~exist(dirs.saveDirProc,'dir')
  mkdir(dirs.saveDirProc)
end

% Assumes we have chan locs file in ~/Documents/MATLAB/mat_mvm/eeg/
files.elecfile = fullfile(dirs.homeDir,'Documents/MATLAB/mat_mvm/eeg/GSN_HydroCel_129_short.sfp');
%locsFormat = 'besa_sfp';
%elec = ft_read_sens(files.elecfile,'fileformat',locsFormat);

% figure printing options - see mm_ft_setSaveDirs for other options
files.saveFigs = 1;
files.figFileExt = 'eps';
if strcmp(files.figFileExt,'eps')
  files.figPrintFormat = '-depsc2';
elseif strcmp(files.figFileExt,'pdf')
  files.figPrintFormat = '-dpdf';
elseif strcmp(files.figFileExt,'png')
  files.figPrintFormat = '-dpng';
end

% directory to save figures
dirs.saveDirFigs = fullfile(dirs.saveDirProc,'figs');
if ~exist(dirs.saveDirFigs,'dir')
  mkdir(dirs.saveDirFigs)
end

%% add NS's artifact information to the event structure
% nsEvFilters.eventValues = exp.eventValues;
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
% for sub = 1:length(exp.subjects)
%   for ses = 1:length(exp.sessions)
%     ns_addArtifactInfo(dirs.dataroot,exp.subjects{sub},exp.sessions{ses},nsEvFilters,0);
%   end
% end

%% Convert the data to FieldTrip structs - excludes NS artifact trials

% cfg_proc = [];
% cfg_proc.method = 'mtmfft';
% cfg_proc.output = 'fourier';
% cfg_proc.foilim = [4 8];
% cfg_proc.tapsmofrq = 5;
% cfg_proc.keeptrials = 'yes';
% cfg_proc.channel = {'e1','e53'};
% cfg_proc.channelcmb = {'e1','e53'};

cfg_proc = [];
cfg_proc.method = 'mtmconvol';
%cfg_proc.output = 'fourier';
cfg_proc.output = 'powandcsd';
cfg_proc.taper = 'dpss';
cfg_proc.foi = 2:1:10;
%cfg_proc.foi = 3:1:9;
cfg_proc.t_ftimwin = 5./cfg_proc.foi;
cfg_proc.tapsmofrq = 0.4*cfg_proc.foi;
cfg_proc.toi = -0.5:0.05:1.5;
%cfg_proc.pad = 'maxperlen';
cfg_proc.keeptrials = 'yes';
cfg_proc.keeptapers = 'no';

% single precision to save space
cfg_proc.precision = 'single';
cfg_proc.channel = {'e1','e53','e122'};
cfg_proc.channelcmb = {'e1','e53';'e1','e122';'e53','e122'};
%cfg_proc.channel = {'e1','e53'};
%cfg_proc.channelcmb = {'e1','e53'};


% cfg_proc = [];
% cfg_proc.method = 'wltconvol';
% cfg_proc.output = 'powandcsd';
% cfg_proc.foi = 4:1:8;
% cfg_proc.toi = -0.5:0.05:1.5;
% cfg.width = 7;
% cfg_proc.keeptrials = 'yes';
% cfg_proc.channel = {'e1','e53'};
% cfg_proc.channelcmb = {'e1','e53'};

cfg_pp = [];

[data_freq,exp.eventValues] = create_ft_struct(cfg_proc,'ft_freqanalysis',cfg_pp,dirs.dataroot,exp.eegFileExt,exp.subjects,exp.sessions,exp.prepost,files.elecfile,exp.sampleRate,exp.eventValues,exp.eventValuesExtra);
%[data_freq,exp.eventValues] = create_ft_struct_parallel(cfg_proc,'ft_freqanalysis',cfg_pp,dirs.dataroot,dirs.saveDirName,exp.eegFileExt,exp.subjects,exp.sessions,exp.prepost,files.elecfile,exp.sampleRate,exp.eventValues,exp.eventValuesExtra);

% % save the structs for loading in later
% saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s_%d_%d_%d_%d.mat',cfg_proc.output,cfg_proc.foi(1),cfg_proc.foi(end),exp.prepost(1)*1000,exp.prepost(2)*1000));
% if ~exist(saveFile,'file')
%   save(saveFile,'data_freq','cfg_proc');
% end


% prestim = -0.5;
% poststim = 1.5;
% timestep = 0.05;
% timewin = .4;
% freqwin = [4 8];
% cfg_proc = [];
% cfg_proc.output ='powandcsd';
% cfg_proc.taper = 'dpss';
% cfg_proc.channel = {'e1','e53'};
% cfg_proc.channelcmb = {'e1','e53'};
% cfg_proc.method          = 'mtmconvol';
% cfg_proc.keeptrials = 'yes';
% cfg_proc.keeptapers = 'no';
% cfg_proc.foi             = (1/timewin):(1/timewin):freqwin(2); % Frequency axis
% cfg_proc.foi             = cfg_proc.foi(cfg_proc.foi>=freqwin(1));
% cfg_proc.t_ftimwin       = zeros(1,length(cfg_proc.foi));
% cfg_proc.t_ftimwin(:)    = timewin; % Time resolution
% cfg_proc.tapsmofrq       = zeros(1,length(cfg_proc.foi)); % Initialize to zero
% cfg_proc.tapsmofrq(:)    = 1/timewin; % Set initial resolution to 1/timewin (i.e. 2.5 Hz) for all frequencis
% cfg_proc.tapsmofrq(cfg_proc.foi>10*(1/timewin))    = 0.1*cfg_proc.foi(cfg_proc.foi>10*(1/timewin));
% cfg_proc.tapsmofrq(cfg_proc.foi>50)                = 5;
% cfg_proc.toi=(prestim+(timewin/2)):timestep:(poststim-(timewin/2)-1/exp.sampleRate); % Time axis


%% Calculate the coherence

cfg_coh = [];
cfg_coh.method = 'coh';
cfg_coh.channelcmb = {'e1','e53';'e1','e122';'e53','e122'};
%cfg_coh.channelcmb = {'e1','e53'};

data_coh = struct();

for evVal = 1:length(exp.eventValues)
  for sub = 1:length(exp.subjects)
    for ses = 1:length(exp.sessions)
      data_coh.(exp.eventValues{evVal}).sub(sub).ses(ses).data = ft_connectivityanalysis(cfg_coh,data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data);
      
      % if the output is powandcsd, the label field doesn't get passed down
      if ~isfield(data_coh.(exp.eventValues{evVal}).sub(sub).ses(ses).data,'label')
        data_coh.(exp.eventValues{evVal}).sub(sub).ses(ses).data.label = data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data.label;
      end
    end
  end
end

cfg_plv = [];
cfg_plv.method = 'plv';
%cfg_plv.channelcmb = {'e1','e53';'e1','e122';'e53','e122'};
cfg_plv.channelcmb = {'e1','e53'};

data_plv = struct();

for evVal = 1:length(exp.eventValues)
  for sub = 1:length(exp.subjects)
    for ses = 1:length(exp.sessions)
      data_plv.(exp.eventValues{evVal}).sub(sub).ses(ses).data = ft_connectivityanalysis(cfg_plv,data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data);
      
      % if the output is powandcsd, the label field doesn't get passed down
      if ~isfield(data_plv.(exp.eventValues{evVal}).sub(sub).ses(ses).data,'label')
        data_plv.(exp.eventValues{evVal}).sub(sub).ses(ses).data.label = data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data.label;
      end
    end
  end
end

%% 2-way ANOVA: Channel x Condition

cfg_ana = [];
cfg_ana.channels = {'e53'};
cfg_ana.cohrefchannel = {'e1'};
cfg_ana.cond = {'RCR','RH'};
%cfg_ana.cond = {'RCR','RHSC','RHSI'};
cfg_ana.latency = [0.3 0.5];
cfg_ana.frequency = [4 8];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;
if length(cfg_ana.channels) > 2 || length(cfg_ana.cond) > 2
  cfg_ana.calcGGHF = 1;
else
  cfg_ana.calcGGHF = 0;
end

fprintf('====================== RMAOV2_mod =========================\n');
fprintf('IV1: Hemisphere (%d;%s), IV2: Condition (%d;%s)\n',length(cfg_ana.channels),sprintf(repmat(' %s',1,length(cfg_ana.channels)),cfg_ana.channels{:}),length(cfg_ana.cond),sprintf(repmat(' %s',1,length(cfg_ana.cond)),cfg_ana.cond{:}));

ses=1;

cfg_ana.anovamat = [];
%for h = 1:length(cfg_ana.channels)
  %cfg_ana.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.channels(h))});
  %cfg_ana.chan = unique(str2double(strrep(cfg_ana.channel,'e','')));
  
  cfg_ana.cohrefidx = find(ismember(data_plv.RCR.sub(1).ses(1).data.label,cfg_ana.cohrefchannel));
  cfg_ana.chanidx = find(ismember(data_plv.RCR.sub(1).ses(1).data.label,cfg_ana.channels));
  
  for c = 1:length(cfg_ana.cond)
    for sub = 1:size(data_plv.(cfg_ana.cond{c}).sub,2)
      cfg_ana.timesel.(cfg_ana.cond{c}) = find(data_plv.(cfg_ana.cond{c}).sub(sub).ses(ses).data.time >= cfg_ana.latency(1) & data_plv.(cfg_ana.cond{c}).sub(sub).ses(ses).data.time <= cfg_ana.latency(2));
      cfg_ana.freqsel.(cfg_ana.cond{c}) = find(data_plv.(cfg_ana.cond{c}).sub(sub).ses(ses).data.freq >= cfg_ana.frequency(1) & data_plv.(cfg_ana.cond{c}).sub(sub).ses(ses).data.freq <= cfg_ana.frequency(2));
      cfg_ana.values.(cfg_ana.cond{c}) = mean(mean(mean(data_plv.(cfg_ana.cond{c}).sub(sub).ses(ses).data.plvspctrm(cfg_ana.cohrefidx,cfg_ana.chanidx,cfg_ana.freqsel.(cfg_ana.cond{c}),cfg_ana.timesel.(cfg_ana.cond{c})),2),3),4);
      %index = h + (c - 1) + (h - 1);
      
      % format for RMAOV2_mod: [data channels cond subNum]
      %cfg_ana.anovamat = cat(1,cfg_ana.anovamat,[cfg_ana.values.(cfg_ana.cond{c}) h c sub]);
      cfg_ana.anovamat = cat(1,cfg_ana.anovamat,[cfg_ana.values.(cfg_ana.cond{c}) c sub]);
    end
  end
%end
%[cfg_ana.p] = RMAOV2_mod(cfg_ana.anovamat,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,cfg_ana.printTable_tex);
[cfg_ana.p] = RMAOV1_mod(cfg_ana.anovamat,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,cfg_ana.printTable_tex);


%% Plot the coherence spectrum

cfg_plot = [];
cfg_plot.xparam = 'time';
cfg_plot.yparam = 'freq';
%cfg_plot.parameter = 'powspctrm';
%cfg_plot.parameter = 'crsspctrm';
cfg_plot.channel = 'e53';
cfg_plot.cohrefchannel = 'e1';
cfg_plot.showlabels       = 'yes';
cfg_plot.colorbar         = 'yes';
figure
cfg_plot.parameter = 'cohspctrm';
ft_singleplotTFR(cfg_plot,data_coh.RCR.sub(1).ses(1).data);

cfg_plot.layout = ft_prepare_layout([],data_freq);
figure
ft_multiplotTFR(cfg_plot,data_coh.RCR.sub(1).ses(1).data);

figure
cfg_plot.parameter = 'plvspctrm';
ft_singleplotTFR(cfg_plot,data_plv.RCR.sub(1).ses(1).data);

%ft_singleplotTFR(cfg_plot,data_freq.RCR.sub(1).ses(1).data);


cfg_plot = [];
cfg_plot.xparam = 'freq';
cfg_plot.xlim = [4 8];

cfg_plot.parameter = 'cohspctrm';
%cfg_plot.xlim = [4 8];
%cfg_plot.ylim = [0 1];
cfg_plot.cohrefchannel = 'e1';
cfg_plot.layout = ft_prepare_layout([],data_freq);
cfg_plot.showlabels = 'yes';
figure;
ft_multiplotER(cfg_plot,data_coh.RCR.sub(1).ses(1).data);

cfg_plot.channel = 'e53';
figure
ft_singleplotER(cfg_plot,data_coh.RCR.sub(1).ses(1).data);
%ft_singleplotTFR(cfg_plot,data_coh.RCR.sub(1).ses(1).data);
title('RCR');
figure
ft_singleplotER(cfg_plot,data_coh.RHSC.sub(1).ses(1).data);
title('RHSI');
figure
ft_singleplotER(cfg_plot,data_coh.RHSI.sub(1).ses(1).data);
title('RHSC');

%% Potential bad channels

% RSRC021 has a screwy channel 117; power in baseline period is really high


%% Test plots to make sure data look ok

% cfg_ft = [];
% cfg_ft.baseline = [-0.3 -0.1];
% cfg_ft.baselinetype = 'absolute';
% if strcmp(cfg_ft.baselinetype,'absolute')
%   cfg_ft.zlim = [-700 700];
% elseif strcmp(cfg_ft.baselinetype,'relative')
%   cfg_ft.zlim = [0 2.0];
% end
% cfg_ft.showlabels = 'yes';
% cfg_ft.colorbar = 'yes';
% cfg_ft.interactive = 'yes';
% cfg_ft.layout = ft_prepare_layout([],data_freq);
% figure
% ft_multiplotTFR(cfg_ft,data_freq.RHSC.sub(1).ses(1).data);
% 
% cfg_ft = [];
% cfg_ft.channel = {'e124'};
% %cfg_ft.channel = {'e117'};
% cfg_ft.baseline = [-0.5 -0.1];
% cfg_ft.baselinetype = 'absolute';
% if strcmp(cfg_ft.baselinetype,'absolute')
%   cfg_ft.zlim = [-700 700];
% elseif strcmp(cfg_ft.baselinetype,'relative')
%   cfg_ft.zlim = [0 1.5];
% end
% cfg_ft.showlabels = 'yes';
% cfg_ft.colorbar = 'yes';
% ft_singleplotTFR(cfg_ft,data_freq.RHSC.sub(1).ses(1).data);
% 
%figure
% badTrials = [];
% for ev = 1:size(data_freq.RHSC.sub(1).ses(1).data.powspctrm,1)
%   if data_freq.RHSC.sub(1).ses(1).data.powspctrm(ev,117,5,20) < -800
%     badTrials = [badTrials ev];
% %     fprintf('%d: Low! %.1f\n',ev,data_freq.RHSC.sub(1).ses(1).data.powspctrm(ev,117,5,20));
% %     clf
% %     imagesc(squeeze(data_freq.RHSC.sub(1).ses(1).data.powspctrm(ev,117,:,:)),[-700 700]);
% %     axis xy
% %     title(['Trial ',num2str(ev)]);
% %     xlabel('Time (ms)');
% %     ylabel('Frequency (Hz)');
% %     set(gca,'XTick',1:5:length(data_freq.RHSC.sub(1).ses(1).data.time));
% %     set(gca,'XTickLabel',data_freq.RHSC.sub(1).ses(1).data.time(1):.25:data_freq.RHSC.sub(1).ses(1).data.time(end));
% %     set(gca,'YTickLabel',data_freq.RHSC.sub(1).ses(1).data.freq);
% %     h = colorbar;
% %     set(get(h,'YLabel'),'string','Power');
% %     keyboard
% %   else
% %     fprintf('%d: %.1f\n',ev,data_freq.RHSC.sub(1).ses(1).data.powspctrm(ev,117,5,20));
% %     continue
%   end
% end
% 
% cfg_bcr = [];
% cfg_bcr.badchannel = {'e117'};
% cfg_bcr.trials = badTrials;
% [interp] = ft_channelrepair(cfg_bcr,data_freq.RHSC.sub(1).ses(1).data);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cfg_ft = [];
% cfg_ft.baseline = [-0.5 -0.1];	
% cfg_ft.baselinetype = 'absolute';
% if strcmp(cfg_ft.baselinetype,'absolute')
%   cfg_ft.zlim = [-700 700];
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
% ft_topoplotTFR(cfg_ft,data_freq.RHSC.sub(1).ses(1).data);
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

%% if already saved and not yet loaded, load the ft_freqanalysis files

if ~exist('data_freq','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('data_freq_%d_%d.mat',exp.prepost(1)*1000,exp.prepost(2)*1000)));
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDirProc,savedFiles(sf).name));
    fprintf('Done.\n');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change in freq relative to baseline using absolute power

%data_freq_old = data_freq;

cfg_fb = [];
cfg_fb.baseline = [-0.5 -0.1];
cfg_fb.baselinetype = 'absolute';

for evVal = 1:length(exp.eventValues)
  for sub = 1:length(exp.subjects)
    for ses = 1:length(exp.sessions)
      data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data);
    end
  end
end

%% decide who to kick out based on trial counts

for evVal = 1:length(exp.eventValues)
  numEv.(exp.eventValues{evVal}) = zeros(length(exp.subjects),length(exp.sessions));
end

if exist('data_freq','var')
  for evVal = 1:length(exp.eventValues)
    for sub = 1:length(exp.subjects)
      for ses = 1:length(exp.sessions)
        numEv.(exp.eventValues{evVal})(sub,ses) = size(ft_findcfg(data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data.cfg,'trl'),1);
      end
    end
  end
elseif exist('ga_freq','var')
  % if we need to recreate the bad subjects list from a grand average file
  goodSub = cell(size(ga_freq.RCR.cfg.previous,2),1);
  for evVal = 1:length(exp.eventValues)
    for sub = 1:size(ga_freq.RCR.cfg.previous,2)
      datafile = ft_findcfg(ga_freq.RCR.cfg.previous{sub},'datafile');
      [~,datafile] = fileparts(datafile);
      dividers = strfind(datafile,'_');
      if isempty(dividers)
        dividers = strfind(datafile,' ');
      end
      goodSub{sub} = datafile(strfind(datafile,exp.name):dividers(1)-1);
      numEv.(exp.eventValues{evVal})(sub,ses) = size(ft_findcfg(ga_freq.(exp.eventValues{evVal}).cfg.previous{sub},'trl'),1);
    end
  end
  exp.badSub = ~ismember(exp.subjects,goodSub);
  [exp.subjects,subInd] = sort(cat(1,goodSub,exp.subjects(exp.badSub)));
  for evVal = 1:length(exp.eventValues)
    numEv.(exp.eventValues{evVal}) = numEv.(exp.eventValues{evVal})(subInd);
  end
  % add the extra event values onto the end of eventValues
  if ~isempty(exp.eventValuesExtra)
    for nVal = 1:length(exp.eventValuesExtra.newValue)
      if ~ismember(exp.eventValuesExtra.newValue{nVal},exp.eventValues)
        exp.eventValues = cat(2,exp.eventValues,exp.eventValuesExtra.newValue{nVal});
      end
    end
    exp.eventValues = sort(exp.eventValues);
  end
end

fprintf('\tRCR\tRHSC\tRHSI\n');
for sub = 1:length(exp.subjects)
  fprintf('%s\t%d\t%d\t%d\n',exp.subjects{sub},numEv.RCR(sub),numEv.RHSC(sub),numEv.RHSI(sub))
end

numEv.thresh = 20;
numEv.lowNum = zeros(length(exp.subjects),length(exp.sessions));
for sub = 1:length(exp.subjects)
  for ses = 1:length(exp.sessions)
    for evVal = 1:length(exp.eventValues)
      if numEv.(exp.eventValues{evVal})(sub,ses) < numEv.thresh
        numEv.lowNum(sub,ses) = 1;
      end
    end
  end
end
numEv.lowNum = logical(numEv.lowNum);

% who we're rejecting
fprintf('Subjects with low trial counts:\n');
exp.subjects(numEv.lowNum)

%exp.badBehSub = {};
% blink > 45%
exp.badBehSub = {'RSRC032';'RSRC040'};
% d' < 0.7
%exp.badBehSub = {'RSRC021','RSRC030','RSRC040','RSRC049'};
% RSRC052 has very positive ERPs in LPS

fprintf('Subjects with bad behavior:\n');
exp.badBehSub

exp.badSub = zeros(length(exp.subjects),length(exp.sessions));
% low event numbers
for ses = 1:length(exp.sessions)
  exp.badSub(:,ses) = sign(sum([numEv.lowNum(:,ses) exp.badSub],2));
end
% bad behavioral subjects
for ses = 1:length(exp.sessions)
  exp.badSub(:,ses) = sign(sum([ismember(exp.subjects,exp.badBehSub) exp.badSub],2));
end
exp.badSub = logical(exp.badSub);

fprintf('Number of events included in EEG analyses (%d subjects; threshold: %d events):\n',sum(~exp.badSub),numEv.thresh);
for evVal = 1:length(exp.eventValues)
  numEv.meanTrials = mean(numEv.(exp.eventValues{evVal})(~exp.badSub));
  numEv.sem = std(numEv.(exp.eventValues{evVal})(~exp.badSub),0,1)/sqrt(sum(~exp.badSub));
  numEv.sd = std(numEv.(exp.eventValues{evVal})(~exp.badSub),0,1);
  fprintf('%s:\tM=%.3f,\tSEM=%.3f,\tSD=%.3f\n',exp.eventValues{evVal},numEv.meanTrials,numEv.sem,numEv.sd);
end

%% set up strings to put in ft_freqgrandaverage

count = 1;
for sub = 1:length(exp.subjects)
  for ses = 1:length(exp.sessions)
    if exp.badSub(sub,ses)
      fprintf('BadSubject: %s\n',exp.subjects{sub});
      continue
    else
      if count == 1
        for evVal = 1:length(exp.eventValues)
          ana.str.(exp.eventValues{evVal}){ses} = sprintf('data_coh.%s.sub(%d).ses(%d).data',exp.eventValues{evVal},sub,ses);
        end
        count = count + 1;
      elseif count > 1
        for evVal = 1:length(exp.eventValues)
          ana.str.(exp.eventValues{evVal}){ses} = sprintf('%s,data_coh.%s.sub(%d).ses(%d).data',ana.str.(exp.eventValues{evVal}){ses},exp.eventValues{evVal},sub,ses);
        end
        count = count + 1;
      end
    end
  end
end

%% get the grand average and add the electrode locations again

cfg_ft = [];
cfg_ft.keepindividual = 'yes';
for ses = 1:length(exp.sessions)
  for evVal = 1:length(exp.eventValues)
    tic
    fprintf('Running ft_freqgrandaverage on %s...',exp.eventValues{evVal});
    ga_freq.(exp.eventValues{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',ana.str.(exp.eventValues{evVal}){ses}));
    fprintf('Done.\n');
    toc
  end
end
ga_freq.elec = data_freq.elec;

%% save grand average file

saveFile = fullfile(dirs.saveDirProc,sprintf('ga_freq_%d_%d_%d_%d.mat',ga_freq.RCR.freq(1),ga_freq.RCR.freq(end),ga_freq.RCR.time(1)*1000,ga_freq.RCR.time(end)*1000));
if ~exist(saveFile,'file')
  save(saveFile,'ga_freq');
end

%% save the analysis details

saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails_freq_%d_%d_%d_%d.mat',ga_freq.RCR.freq(1),ga_freq.RCR.freq(end),ga_freq.RCR.time(1)*1000,ga_freq.RCR.time(end)*1000));
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
cfg_ft.zlim = [-400 400];
%elseif strcmp(cfg_ft.baselinetype,'relative')
%  cfg_ft.zlim = [0 2.0];
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ga_freq);
figure
ft_multiplotTFR(cfg_ft,ga_freq.RHSC);
figure
ft_multiplotTFR(cfg_ft,ga_freq.RHSI);
figure
ft_multiplotTFR(cfg_ft,ga_freq.RCR);

%% subplot with all subjects' ERPs

ses = 1;
cfg_plot = [];
cfg_plot.roi = {'LAS','RAS'};
%cfg_plot.roi = {'LPS','RPS'};
%cfg_plot.roi = {'RAS'};
%cfg_plot.roi = {'LPS'};
cfg_plot.chan_str = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
cfg_plot.chan = str2double(strrep(cfg_plot.chan_str,'e',''));
cfg_plot.numCols = 5;
cfg_plot.excludeBadSub = 0;
if cfg_plot.excludeBadSub
  cfg_plot.numRows = ceil((length(exp.subjects) - sum(exp.badSub)) / cfg_plot.numCols);
else
  cfg_plot.numRows = ceil((length(exp.subjects)) / cfg_plot.numCols);
end
figure
for sub = 1:length(exp.subjects)
  subplot(cfg_plot.numRows,cfg_plot.numCols,sub);
  plot(data_freq.RHSC.sub(sub).ses(ses).data.time,mean(data_freq.RHSC.sub(sub).ses(ses).data.avg(cfg_plot.chan,:),1),'r');
  hold on
  plot(data_freq.RHSI.sub(sub).ses(ses).data.time,mean(data_freq.RHSI.sub(sub).ses(ses).data.avg(cfg_plot.chan,:),1),'b');
  plot(data_freq.RH.sub(sub).ses(ses).data.time,mean(data_freq.RH.sub(sub).ses(ses).data.avg(cfg_plot.chan,:),1),'g');
  plot(data_freq.RCR.sub(sub).ses(ses).data.time,mean(data_freq.RCR.sub(sub).ses(ses).data.avg(cfg_plot.chan,:),1),'k');
  if exp.badSub(sub,ses)
    fprintf('BadSubject: %s\n',exp.subjects{sub});
    title(sprintf('*BAD*: %s; CR:%d, SC:%d, SI:%d',exp.subjects{sub},size(data_freq.RCR.sub(sub).ses(ses).data.powspctrm,1),size(data_freq.RHSC.sub(sub).ses(ses).data.powspctrm,1),size(data_freq.RHSI.sub(sub).ses(ses).data.powspctrm,1)));
  else
    title(sprintf('%s; CR:%d, SC:%d, SI:%d',exp.subjects{sub},size(data_freq.RCR.sub(sub).ses(ses).data.powspctrm,1),size(data_freq.RHSC.sub(sub).ses(ses).data.powspctrm,1),size(data_freq.RHSI.sub(sub).ses(ses).data.powspctrm,1)));
  end
  %ylim([0 1.9e-13]);
  axis([-0.2 1.5 -10 10]);
  hold off
end
subplot(cfg_plot.numRows,cfg_plot.numCols,length(exp.subjects) + 1);
text(0.5,0.9,sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}),'color','k')
text(0.5,0.7,'H-SC','color','r') ;
text(0.5,0.5,'H-SI','color','b')
text(0.5,0.3,'H','color','g')
text(0.5,0.1,'CR','color','k')
axis off

%% plot the conditions

% Contrast

ga_freq_RHSCvsRHSI = ga_freq.RHSC;
ga_freq_RHSCvsRHSI.powspctrm = ga_freq.RHSC.powspctrm - ga_freq.RHSI.powspctrm;

cfg_ft = [];
cfg_ft.zlim = [-400 400];
%cfg_ft.zlim = [0 2.0];
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ga_freq);
figure
ft_multiplotTFR(cfg_ft,ga_freq_RHSCvsRHSI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast

cfg_ft = [];
%cfg_ft.xlim = [.6 1.2];
cfg_ft.xlim = [.5 .8];
%cfg_ft.xlim = [.3 .5];
cfg_ft.ylim = [4 8];
cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.interactive = 'yes';
cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'yes';
cfg_ft.layout = ft_prepare_layout([],ga_freq);
figure
ft_topoplotTFR(cfg_ft,ga_freq_RHSCvsRHSI);



% cfg_ft = [];
% cfg_ft.showlabels = 'yes';
% cfg_ft.fontsize = 9;
% cfg_ft.elec = ga_freq.elec;
% cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% cfg_ft.linestyle = {'-','--','-.'};
% 
% cfg_plot = [];
% %cfg_plot.fillcolor = [.8,.8,.8];
% %cfg_plot.fillalpha = 0.3;
% %cfg_plot.filledge = 'none';
% cfg_plot.minMaxVolt = [-4 3; -1 6];
% cfg_plot.latency = [-0.2 2.0];
% cfg_plot.plotLegend = 1;
% cfg_plot.plot_rois = {{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
% %cfg_plot.plot_rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% 
% for i = 1:length(cfg_plot.plot_rois)
%   cfg_plot.roi = cfg_plot.plot_rois{i};
%   cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
%   figure
%   %hold on
%   %ft_singleplotER(cfg_ft,ga_freq.RHSC,ga_freq.RHSI,ga_freq.RH,ga_freq.RCR);
%   ft_singleplotER(cfg_ft,ga_freq.RHSC,ga_freq.RHSI,ga_freq.RCR);
%   hold on
%   plot([cfg_plot.latency(1) cfg_plot.latency(2)],[0 0],'k--'); % horizontal
%   if ismember('LAS',cfg_plot.roi) || ismember('RAS',cfg_plot.roi)
%     plot([0 0],[cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)],'k--'); % vertical
%     if cfg_plot.plotLegend
%       %legend('H-SC','H-SI','H','CR','Location','SouthEast');
%       legend('H-SC','H-SI','CR','Location','SouthEast');
%     end
%     
%     cfg_plot.vert_latency = [0.3 0.5];
%     %h = fill([cfg_plot.vert_latency(1),cfg_plot.vert_latency(2),cfg_plot.vert_latency(2),cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(1,1),cfg_plot.minMaxVolt(1,1),cfg_plot.minMaxVolt(1,2),cfg_plot.minMaxVolt(1,2)],fillcolor);
%     %set(h,'FaceAlpha',fillalpha);
%     %set(h,'EdgeColor',filledge)
%     plot([cfg_plot.vert_latency(1) cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)],'k'); % vertical
%     plot([cfg_plot.vert_latency(2) cfg_plot.vert_latency(2)],[cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)],'k'); % vertical
%     
%     axis([cfg_plot.latency(1) cfg_plot.latency(2) cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)])
%   elseif ismember('LPS',cfg_plot.roi) || ismember('RPS',cfg_plot.roi)
%     plot([0 0],[cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)],'k--'); % vertical
%     if cfg_plot.plotLegend
%       %legend('H-SC','H-SI','H','CR','Location','NorthEast');
%       legend('H-SC','H-SI','CR','Location','NorthEast');
%     end
%     
%     cfg_plot.vert_latency = [0.5 0.8];
%     %h = fill([cfg_plot.vert_latency(1),cfg_plot.vert_latency(2),cfg_plot.vert_latency(2),cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(2,1),cfg_plot.minMaxVolt(2,1),cfg_plot.minMaxVolt(2,2),cfg_plot.minMaxVolt(2,2)],fillcolor);
%     %set(h,'FaceAlpha',fillalpha);
%     %set(h,'EdgeColor',filledge);
%     plot([cfg_plot.vert_latency(1) cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)],'k'); % vertical
%     plot([cfg_plot.vert_latency(2) cfg_plot.vert_latency(2)],[cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)],'k'); % vertical
%     
%     axis([cfg_plot.latency(1) cfg_plot.latency(2) cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)])
%   end
%   xlabel('Time (s)');
%   ylabel('Voltage (\muV)');
%   set(gcf,'Name',sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}))
%   %title(sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}))
%   %axis ij % negative up
%   publishfig(gcf,1);
%   if cfg_plot.plotLegend
%     cfg_plot.legend_str = '_legend';
%   else
%     cfg_plot.legend_str = '';
%   end
%   if files.saveFigs
%     cfg_plot.figfilename = sprintf('erp_ga_%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_plot.latency(1)*1000,cfg_plot.latency(2)*1000,cfg_plot.legend_str,files.figFileExt);
%     print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
%   end
% end
% 
% % cfg_ft = [];
% % cfg_plot.roi = {'LPS','RPS'};
% % cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
% % cfg_ft.linewidth = 2;
% % cfg_ft.graphcolor = 'rbgk';
% % cfg_ft.linestyle = {'-','--','-','-.'};
% % figure
% % ft_singleplotER(cfg_ft,ga_freq.RHSC,ga_freq.RHSI,ga_freq.RH,ga_freq.RCR);
% % hold on
% % legend('H-SC','H-SI','H','CR','Location','NorthEast');
% % plot([cfg_plot.latency(1) cfg_plot.latency(2)],[0 0],'k--'); % horizontal
% % plot([0 0],[cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)],'k--'); % vertical
% % axis([cfg_plot.latency(1) cfg_plot.latency(2) cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)])
% % xlabel('Time (ms)');
% % ylabel('Voltage (\muV)');
% % %axis ij % negative up
% % publishfig(gcf,1);
% % if files.saveFigs
% %   cfg_plot.figfilename = sprintf('ga_erp_%s%0.1f_%.1f.eps',sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_plot.latency(1),cfg_plot.latency(2));
% %   print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
% % end
% 
% %ft_multiplotER(cfg_ft,ga_freq.RCR,ga_freq.RHSC,ga_freq.RHSI);

%% Plot average time contrast topoplots

cfg_ft = [];
cfg_ft.interactive = 'no';
%cfg_ft.colorbar = 'yes';
cfg_ft.layout = ft_prepare_layout(cfg_ft,ga_freq);
cfg_ft.highlight = 'on';
cfg_ft.highlightsize = 10;
%cfg_ft.comment = 'xlim';
%cfg_ft.commentpos = 'title';
cfg_ft.comment = 'no';
%cfg_ft.colormap = colormap('jet'); % default; blue to red
cfg_ft.colormap = colormap('hot'); % dark to light; better for b&w printers

cfg_plot = [];
cfg_plot.minMaxVolt = [-1 1];
cfg_plot.plotTitle = 1;
cfg_plot.plotColorbar = 1;
if cfg_plot.plotColorbar
  cfg_plot.colorbar_str = '_cb';
else
  cfg_plot.colorbar_str = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAS + RAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cont_topo = [];

% create contrast
cont_topo.RHvsRCR = ga_freq.RH;
cont_topo.RHvsRCR.avg = ga_freq.RH.avg - ga_freq.RCR.avg;
cont_topo.RHvsRCR.individual = ga_freq.RH.individual - ga_freq.RCR.individual;
% make a plot
figure
cfg_plot.roi = {'LAS','RAS'};
cfg_plot.cond = {'H','CR'};
cfg_ft.xlim = [0.3 0.5];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHvsRCR);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Hits - Correct Rejections');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

% create contrast
cont_topo.RHSCvsRHSI = ga_freq.RHSC;
cont_topo.RHSCvsRHSI.avg = ga_freq.RHSC.avg - ga_freq.RHSI.avg;
cont_topo.RHSCvsRHSI.individual = ga_freq.RHSC.individual - ga_freq.RHSI.individual;
% make a plot
figure
cfg_plot.roi = {'LAS','RAS'};
cfg_plot.cond = {'H-SC','H-SI'};
cfg_ft.xlim = [0.3 0.5];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHSCvsRHSI);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Source Correct - Source Incorrect');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

% create contrast
cont_topo.RHSCvsRCR = ga_freq.RHSC;
cont_topo.RHSCvsRCR.avg = ga_freq.RHSC.avg - ga_freq.RCR.avg;
cont_topo.RHSCvsRCR.individual = ga_freq.RHSC.individual - ga_freq.RCR.individual;
% make a plot
figure
cfg_plot.roi = {'LAS','RAS'};
cfg_plot.cond = {'H-SC','CR'};
cfg_ft.xlim = [0.3 0.5];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHSCvsRCR);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Source Correct - Correct Rejections');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

% create contrast
cont_topo.RHSIvsRCR = ga_freq.RHSI;
cont_topo.RHSIvsRCR.avg = ga_freq.RHSI.avg - ga_freq.RCR.avg;
cont_topo.RHSIvsRCR.individual = ga_freq.RHSI.individual - ga_freq.RCR.individual;
% make a plot
figure
cfg_plot.roi = {'LAS','RAS'};
cfg_plot.cond = {'H-SI','CR'};
cfg_ft.xlim = [0.3 0.5];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHSIvsRCR);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Source Incorrect - Correct Rejections');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LPS + RPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cont_topo = [];

% create contrast
cont_topo.RHvsRCR = ga_freq.RH;
cont_topo.RHvsRCR.avg = ga_freq.RH.avg - ga_freq.RCR.avg;
cont_topo.RHvsRCR.individual = ga_freq.RH.individual - ga_freq.RCR.individual;
% make a plot
figure
cfg_plot.roi = {'LPS','RPS'};
cfg_plot.cond = {'H','CR'};
cfg_ft.xlim = [0.5 0.8];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHvsRCR);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Hits - Correct Rejections');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

% create contrast
cont_topo.RHSCvsRHSI = ga_freq.RHSC;
cont_topo.RHSCvsRHSI.avg = ga_freq.RHSC.avg - ga_freq.RHSI.avg;
cont_topo.RHSCvsRHSI.individual = ga_freq.RHSC.individual - ga_freq.RHSI.individual;
% make a plot
figure
cfg_plot.roi = {'LPS','RPS'};
cfg_plot.cond = {'H-SC','H-SI'};
cfg_ft.xlim = [0.5 0.8];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHSCvsRHSI);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Source Correct - Source Incorrect');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

% create contrast
cont_topo.RHSCvsRCR = ga_freq.RHSC;
cont_topo.RHSCvsRCR.avg = ga_freq.RHSC.avg - ga_freq.RCR.avg;
cont_topo.RHSCvsRCR.individual = ga_freq.RHSC.individual - ga_freq.RCR.individual;
% make a plot
figure
cfg_plot.roi = {'LPS','RPS'};
cfg_plot.cond = {'H-SC','CR'};
cfg_ft.xlim = [0.5 0.8];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHSCvsRCR);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Source Correct - Correct Rejections');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

% create contrast
cont_topo.RHSIvsRCR = ga_freq.RHSI;
cont_topo.RHSIvsRCR.avg = ga_freq.RHSI.avg - ga_freq.RCR.avg;
cont_topo.RHSIvsRCR.individual = ga_freq.RHSI.individual - ga_freq.RCR.individual;
% make a plot
figure
cfg_plot.roi = {'LPS','RPS'};
cfg_plot.cond = {'H-SI','CR'};
cfg_ft.xlim = [0.5 0.8];
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
ft_topoplotER(cfg_ft,cont_topo.RHSIvsRCR);
set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
if cfg_plot.plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if cfg_plot.plotTitle
  title('Source Incorrect - Correct Rejections');
end
publishfig(gcf,0);
if files.saveFigs
  cfg_plot.figfilename = sprintf('topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
  print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
end

%% correlations

cfg_ana = [];
% get the electrode numbers
cfg_ana.elecGroupsNum = cell(size(ana.elecGroups));
for eGrp = 1:length(ana.elecGroups)
  cfg_ana.elecGroupsNum{eGrp} = str2double(strrep(ana.elecGroups{eGrp},'e',''));
end

% accuracy
cfg_ana.d_item = abs([0.6806 1.5094 1.7587 1.4347 1.8115 2.2193 1.7682 1.5188 1.5516 0.5009 1.6972 2.0737 1.1155 2.7689 0.9843 1.509 1.7668 1.1046 2.3913 0.1865 2.6979 1.5103 1.3333 2.1747 1.5594 1.0772 1.279 0.4565 0.7866 2.8701 2.012 1.7477 1.1253 1.5221 1.488 0.7433 1.2955 1.7184]);
cfg_ana.d_source = abs([0.9133 1.1439 2.0831 1.683 2.2928 2.4325 1.9979 1.6057 1.3698 1.0691 1.3502 1.8977 0.9357 2.4034 1.1257 1.9048 2.0171 1.168 2.6502 0.0421 2.8527 1.9001 1.2099 1.9649 2.0527 0.9707 1.5696 0.2378 1.3907 2.8484 2.0977 1.637 1.1908 2.2677 1.9213 0.9688 1.1103 1.4801]);

cfg_ana.anaComponent = 'FN400';
%cfg_ana.anaComponent = 'PON';

if strcmp(cfg_ana.anaComponent,'FN400')
  % FN400
  cfg_ana.timeMS = [300 500];
  cfg_ana.roi = {'LAS','RAS'};
  %cfg_ana.roi = {'LAS'};
  %cfg_ana.roi = {'RAS'};
  cfg_ana.chans = cat(2,cfg_ana.elecGroupsNum{ismember(ana.elecGroupsStr,cfg_ana.roi)});
elseif strcmp(cfg_ana.anaComponent,'PON')
  % P O/N
  cfg_ana.timeMS = [500 800];
  cfg_ana.roi = {'LPS','RPS'};
  %cfg_ana.roi = {'LPS'};
  %cfg_ana.roi = {'RPS'};
  cfg_ana.chans = cat(2,cfg_ana.elecGroupsNum{ismember(ana.elecGroupsStr,cfg_ana.roi)});
end

cfg_ana.timeSamp = abs(exp.prepost(1) * exp.sampleRate) + [((cfg_ana.timeMS(1) / (1000/exp.sampleRate)) + 1) ((cfg_ana.timeMS(2) / (1000/exp.sampleRate)) + 1)];

% get the voltage for the correct region
for evVal = 1:length(exp.eventValues)
  cfg_ana.values.(exp.eventValues{evVal}) = [];
end
ses = 1;
for sub = 1:length(exp.subjects)
  if exp.badSub(sub,ses)
    fprintf('BadSubject: %s\n',exp.subjects{sub});
    continue
  else
    for evVal = 1:length(exp.eventValues)
      cfg_ana.values.(exp.eventValues{evVal}) = [cfg_ana.values.(exp.eventValues{evVal}) mean(mean(mean(data_freq.(exp.eventValues{evVal}).sub(sub).ses(ses).data.powspctrm(:,cfg_ana.chans,cfg_ana.timeSamp(1):cfg_ana.timeSamp(2)),3),2),1)];
    end
  end
end

% correlate accuracy with voltage differences

cfg_plot = [];
% set up how the lines will look
cfg_plot.linespec = 'ko';
cfg_plot.linewidth = 1.5;
cfg_plot.marksize = 8;
cfg_plot.markcolor = 'w';

if strcmp(cfg_ana.anaComponent,'FN400')
  % FN400: item accuracy with H - CR
  fprintf('FN400\n');
  
  cfg_plot.testType = 'item';
  cfg_plot.cond = {'H','CR'};
  x1 = cfg_ana.d_item(~exp.badSub);
  y1 = (cfg_ana.values.RH - cfg_ana.values.RCR);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Item d'': %s - %s',cfg_plot.cond{1},cfg_plot.cond{2}));
  xlabel('Item d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cfg_plot.cond{1},cfg_plot.cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  if files.saveFigs
    cfg_plot.figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',cfg_plot.testType,sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ana.timeMS(1),cfg_ana.timeMS(2),files.figFileExt);
    print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
  end
  fprintf('Item d'' with H - CR: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  cfg_plot.testType = 'source';
  cfg_plot.cond = {'H-SC','H-SI'};
  x1 = cfg_ana.d_source(~exp.badSub);
  y1 = (cfg_ana.values.RHSC - cfg_ana.values.RHSI);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Source d'': %s - %s',cfg_plot.cond{1},cfg_plot.cond{2}));
  xlabel('Source d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cfg_plot.cond{1},cfg_plot.cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  if files.saveFigs
    cfg_plot.figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',cfg_plot.testType,sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ana.timeMS(1),cfg_ana.timeMS(2),files.figFileExt);
    print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
  end
  fprintf('Source d'' with HSC - HSI: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
elseif strcmp(cfg_ana.anaComponent,'PON')
  % P O/N: source accuracy with HSC - HSI
  fprintf('P O/N\n');
  
  cfg_plot.testType = 'source';
  cfg_plot.cond = {'H-SC','H-SI'};
  x1 = cfg_ana.d_source(~exp.badSub);
  y1 = (cfg_ana.values.RHSC - cfg_ana.values.RHSI);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Source d'': %s - %s',cfg_plot.cond{1},cfg_plot.cond{2}));
  xlabel('Source d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cfg_plot.cond{1},cfg_plot.cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  if files.saveFigs
    cfg_plot.figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',cfg_plot.testType,sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ana.timeMS(1),cfg_ana.timeMS(2),files.figFileExt);
    print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
  end
  fprintf('Source d'' with HSC - HSI: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
end

%% descriptive statistics: ttest

% fieldtrip ttest

cfg_ana = [];
cfg_ana.roi = {'LAS','RAS'};
%cfg_ana.roi = {'LAS'};
%cfg_ana.roi = {'RAS'};
cfg_ft = [];
cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi)});
cfg_ft.latency = [0.6 1.0];
cfg_ft.frequency = [4 8];

% cfg_ana.roi = {'LPS','RPS'};
% %cfg_ana.roi = {'LPS'};
% %cfg_ana.roi = {'RPS'};
% cfg_ft = [];
% cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi)});
% cfg_ft.latency     = [0.5 0.8];
% cfg_ft.frequency = [4 8];

cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.avgoverfreq = 'yes';
cfg_ft.parameter   = 'powspctrm';
cfg_ft.method      = 'analytic';
cfg_ft.statistic   = 'depsamplesT';
cfg_ft.alpha       = 0.05;
cfg_ft.correctm    = 'no';

cfg_ana.numSub = length(exp.subjects) - sum(exp.badSub);
cfg_ft.design(1,1:2*cfg_ana.numSub) = [ones(1,cfg_ana.numSub) 2*ones(1,cfg_ana.numSub)];
cfg_ft.design(2,1:2*cfg_ana.numSub) = [1:cfg_ana.numSub 1:cfg_ana.numSub];
cfg_ft.ivar = 1; % the 1st row in cfg_ft.design contains the independent variable
cfg_ft.uvar = 2; % the 2nd row in cfg_ft.design contains the subject number

% get times, voltage, SEM

cfg_ana.chan = unique(str2double(strrep(cfg_ft.channel,'e','')));
for evVal = 1:length(exp.eventValues)
  cfg_ana.timesel.(exp.eventValues{evVal}) = find(ga_freq.(exp.eventValues{evVal}).time >= cfg_ft.latency(1) & ga_freq.(exp.eventValues{evVal}).time <= cfg_ft.latency(2));
  cfg_ana.freqsel.(exp.eventValues{evVal}) = find(ga_freq.(exp.eventValues{evVal}).freq >= cfg_ft.frequency(1) & ga_freq.(exp.eventValues{evVal}).freq <= cfg_ft.frequency(2));
  cfg_ana.values.(exp.eventValues{evVal}) = mean(mean(mean(ga_freq.(exp.eventValues{evVal}).powspctrm(:,cfg_ana.chan,cfg_ana.freqsel.(exp.eventValues{evVal}),cfg_ana.timesel.(exp.eventValues{evVal})),2),3),4);
  cfg_ana.sem.(exp.eventValues{evVal}) = std(cfg_ana.values.(exp.eventValues{evVal}))/sqrt(length(cfg_ana.values.(exp.eventValues{evVal})));
end

cfg_ana.goodSub = exp.subjects(~exp.badSub);
% ga
fprintf('Cat\t%s\t%s\t%s\t%s\n','CR','H-SC','H-SI','H');
fprintf('GA\t%.3f\t%.3f\t%.3f\t%.3f\n\n',mean(cfg_ana.values.RCR,1),mean(cfg_ana.values.RHSC,1),mean(cfg_ana.values.RHSI,1),mean(cfg_ana.values.RH,1));
% sub avg
for sub = 1:length(cfg_ana.goodSub)
fprintf('Cat\t%s\t%s\t%s\t%s\n','CR','H-SC','H-SI','H');
  fprintf('%s\t%.3f\t%.3f\t%.3f\t%.3f\n',cfg_ana.goodSub{sub},cfg_ana.values.RCR(sub),cfg_ana.values.RHSC(sub),cfg_ana.values.RHSI(sub),cfg_ana.values.RH(sub));
end

cfg_ana.RHvsRCR = ft_freqstatistics(cfg_ft,ga_freq.RH,ga_freq.RCR);
cfg_ana.RHSCvsRHSI = ft_freqstatistics(cfg_ft,ga_freq.RHSC,ga_freq.RHSI);
cfg_ana.RHSCvsRCR = ft_freqstatistics(cfg_ft,ga_freq.RHSC,ga_freq.RCR);
cfg_ana.RHSIvsRCR = ft_freqstatistics(cfg_ft,ga_freq.RHSI,ga_freq.RCR);
fprintf('H   (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(cfg_ana.values.RH,1),cfg_ana.sem.RH,mean(cfg_ana.values.RCR,1),cfg_ana.sem.RCR,cfg_ana.RHvsRCR.df,cfg_ana.RHvsRCR.stat,cfg_ana.RHvsRCR.prob);
fprintf('HSC (M=%.3f; SEM=%.3f) vs HSI (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(cfg_ana.values.RHSC,1),cfg_ana.sem.RHSC,mean(cfg_ana.values.RHSI,1),cfg_ana.sem.RHSI,cfg_ana.RHSCvsRHSI.df,cfg_ana.RHSCvsRHSI.stat,cfg_ana.RHSCvsRHSI.prob);
fprintf('HSC (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(cfg_ana.values.RHSC,1),cfg_ana.sem.RHSC,mean(cfg_ana.values.RCR,1),cfg_ana.sem.RCR,cfg_ana.RHSCvsRCR.df,cfg_ana.RHSCvsRCR.stat,cfg_ana.RHSCvsRCR.prob);
fprintf('HSI (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(cfg_ana.values.RHSI,1),cfg_ana.sem.RHSI,mean(cfg_ana.values.RCR,1),cfg_ana.sem.RCR,cfg_ana.RHSIvsRCR.df,cfg_ana.RHSIvsRCR.stat,cfg_ana.RHSIvsRCR.prob);

% plot to see the effect in each subject
figure
plot([cfg_ana.values.RHSC cfg_ana.values.RHSI]','o-');
xlim([0.5 2.5])
title(sprintf('%.1fs--%.1fs, %s',cfg_ft.latency(1),cfg_ft.latency(2),sprintf(repmat('%s ',1,length(cfg_ana.roi)),cfg_ana.roi{:})));
ylabel('Microvolts (\muV)');
set(gca,'XTickLabel',{'','H-SC','','H-SI',''})
figure
plot([cfg_ana.values.RH cfg_ana.values.RCR]','o-');
xlim([0.5 2.5])
title(sprintf('%.1fs--%.1fs, %s',cfg_ft.latency(1),cfg_ft.latency(2),sprintf(repmat('%s ',1,length(cfg_ana.roi)),cfg_ana.roi{:})));
ylabel('Microvolts (\muV)');
set(gca,'XTickLabel',{'','H','','CR',''})

% % matlab dependent samples ttest
% RHSCvsRHSI = cfg_ana.values.RHSC - cfg_ana.values.RHSI;
% [h,p,ci,stats] = ttest(RHSCvsRHSI, 0, 0.05); % H0: mean = 0, alpha 0.05

%% mean amplitude line plots

cfg_plot = [];

if ismember('LAS',cfg_ana.roi) || ismember('RAS',cfg_ana.roi)
  cfg_plot.conditions = {{'H','CR'}, {'H-SC','H-SI','CR'}};
  cfg_plot.yminmax = [-5 -1];
elseif ismember('LPS',cfg_ana.roi) || ismember('RPS',cfg_ana.roi)
  cfg_plot.conditions = {{'H-SC','H-SI','CR'}};
  cfg_plot.yminmax = [2 6];
end
%cfg_plot.yminmax = [floor(min([mean(cfg_ana.values.RHSC,1),mean(cfg_ana.values.RHSI,1),mean(cfg_ana.values.RCR,1)])) ceil(max([mean(cfg_ana.values.RHSC,1),mean(cfg_ana.values.RHSI,1),mean(cfg_ana.values.RCR,1)]))];

% set up how the lines will look
cfg_plot.linewidth = 2;
cfg_plot.marksize = 10;
cfg_plot.linespec = 'k--o';
cfg_plot.markcolor = 'w';
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;

for cnd = 1:length(cfg_plot.conditions)
  cfg_plot.cond = cfg_plot.conditions{cnd};
  
  figure
  if length(cfg_plot.cond) == 3
    % lines
    plot([mean(cfg_ana.values.RHSC,1),mean(cfg_ana.values.RHSI,1),mean(cfg_ana.values.RCR,1)],cfg_plot.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    % errorbars
    h = errorbar(1,mean(cfg_ana.values.RHSC,1),cfg_ana.sem.RHSC,cfg_plot.linespec,'LineWidth',cfg_plot.errwidth);
    % remove errorbar ends
    if cfg_plot.removeErrBarEnds
      c = get(h,'Children');
      xdata = get(c(2),'XData');
      ydata = get(c(2),'YData');
      xdata(cfg_plot.errBarEndMarkerInd) = NaN;
      ydata(cfg_plot.errBarEndMarkerInd) = NaN;
      set(c(2),'XData',xdata);
      set(c(2),'YData',ydata);
      set(h,'Children',c);
    end
    
    h = errorbar(2,mean(cfg_ana.values.RHSI,1),cfg_ana.sem.RHSI,cfg_plot.linespec,'LineWidth',cfg_plot.errwidth);
    % remove errorbar ends
    if cfg_plot.removeErrBarEnds
      c = get(h,'Children');
      xdata = get(c(2),'XData');
      ydata = get(c(2),'YData');
      xdata(cfg_plot.errBarEndMarkerInd) = NaN;
      ydata(cfg_plot.errBarEndMarkerInd) = NaN;
      set(c(2),'XData',xdata);
      set(c(2),'YData',ydata);
      set(h,'Children',c);
    end
    
    h = errorbar(3,mean(cfg_ana.values.RCR,1),cfg_ana.sem.RCR,cfg_plot.linespec,'LineWidth',cfg_plot.errwidth);
    % remove errorbar ends
    if cfg_plot.removeErrBarEnds
      c = get(h,'Children');
      xdata = get(c(2),'XData');
      ydata = get(c(2),'YData');
      xdata(cfg_plot.errBarEndMarkerInd) = NaN;
      ydata(cfg_plot.errBarEndMarkerInd) = NaN;
      set(c(2),'XData',xdata);
      set(c(2),'YData',ydata);
      set(h,'Children',c);
    end
    
    % markers
    plot(1,mean(cfg_ana.values.RHSC,1),cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
    plot(2,mean(cfg_ana.values.RHSI,1),cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
    plot(3,mean(cfg_ana.values.RCR,1),cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  elseif length(cfg_plot.cond) == 2
    % lines
    plot([mean(cfg_ana.values.RH,1),mean(cfg_ana.values.RCR,1)],cfg_plot.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    % errorbars
    h = errorbar(1,mean(cfg_ana.values.RH,1),cfg_ana.sem.RH,cfg_plot.linespec,'LineWidth',cfg_plot.errwidth);
    % remove errorbar ends
    if cfg_plot.removeErrBarEnds
      c = get(h,'Children');
      xdata = get(c(2),'XData');
      ydata = get(c(2),'YData');
      xdata(cfg_plot.errBarEndMarkerInd) = NaN;
      ydata(cfg_plot.errBarEndMarkerInd) = NaN;
      set(c(2),'XData',xdata);
      set(c(2),'YData',ydata);
      set(h,'Children',c);
    end
    
    h = errorbar(2,mean(cfg_ana.values.RCR,1),cfg_ana.sem.RCR,cfg_plot.linespec,'LineWidth',cfg_plot.errwidth);
    % remove errorbar ends
    if cfg_plot.removeErrBarEnds
      c = get(h,'Children');
      xdata = get(c(2),'XData');
      ydata = get(c(2),'YData');
      xdata(cfg_plot.errBarEndMarkerInd) = NaN;
      ydata(cfg_plot.errBarEndMarkerInd) = NaN;
      set(c(2),'XData',xdata);
      set(c(2),'YData',ydata);
      set(h,'Children',c);
    end
    
    % markers
    plot(1,mean(cfg_ana.values.RH,1),cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
    plot(2,mean(cfg_ana.values.RCR,1),cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  end
  % make it look good
  axis([.5 (length(cfg_plot.cond) + .5) cfg_plot.yminmax(1) cfg_plot.yminmax(2)])
  xlabel('Condition');
  ylabel('Voltage (\muV)');
  set(gca,'XTick',(1:length(cfg_plot.cond)))
  set(gca,'XTickLabel',cfg_plot.cond)
  set(gca,'YTick',(cfg_plot.yminmax(1):.5:cfg_plot.yminmax(2)))
  axis square
  publishfig(gcf,0);
  if files.saveFigs
    cfg_plot.figfilename = sprintf('line_ga_%s%s%d_%d.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000,files.figFileExt);
    print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
  end
end

%% 2-way ANOVA: Hemisphere x Condition: FN400

cfg_ana = [];
cfg_ana.hemi = {'LAS','RAS'};
cfg_ana.cond = {'RCR','RH'};
%cfg_ana.cond = {'RCR','RHSC','RHSI'};
cfg_ana.latency = [0.3 0.5];
cfg_ana.frequency = [4 8];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;
if length(cfg_ana.hemi) > 2 || length(cfg_ana.cond) > 2
  cfg_ana.calcGGHF = 1;
else
  cfg_ana.calcGGHF = 0;
end

fprintf('====================== RMAOV2_mod =========================\n');
fprintf('IV1: Hemisphere (%d;%s), IV2: Condition (%d;%s)\n',length(cfg_ana.hemi),sprintf(repmat(' %s',1,length(cfg_ana.hemi)),cfg_ana.hemi{:}),length(cfg_ana.cond),sprintf(repmat(' %s',1,length(cfg_ana.cond)),cfg_ana.cond{:}));

cfg_ana.anovamat = [];
for h = 1:length(cfg_ana.hemi)
  cfg_ana.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.hemi(h))});
  cfg_ana.chan = unique(str2double(strrep(cfg_ana.channel,'e','')));
  
  for c = 1:length(cfg_ana.cond)
    cfg_ana.timesel.(cfg_ana.cond{c}) = find(ga_freq.(cfg_ana.cond{c}).time >= cfg_ana.latency(1) & ga_freq.(cfg_ana.cond{c}).time <= cfg_ana.latency(2));
    cfg_ana.freqsel.(cfg_ana.cond{c}) = find(ga_freq.(cfg_ana.cond{c}).freq >= cfg_ft.frequency(1) & ga_freq.(cfg_ana.cond{c}).freq <= cfg_ft.frequency(2));
    cfg_ana.values.(cfg_ana.cond{c}) = mean(mean(mean(ga_freq.(cfg_ana.cond{c}).powspctrm(:,cfg_ana.chan,cfg_ana.freqsel.(cfg_ana.cond{c}),cfg_ana.timesel.(cfg_ana.cond{c})),2),3),4);
    
    %index = h + (c - 1) + (h - 1);
    for sub = 1:size(cfg_ana.values.(cfg_ana.cond{c}),1)
      
      % format for RMAOV2_mod: [data hemi cond subNum]
      cfg_ana.anovamat = cat(1,cfg_ana.anovamat,[cfg_ana.cfg_ana.values.(cfg_ana.cond{c})(sub) h c sub]);
    end
  end
end
[cfg_ana.p] = RMAOV2_mod(cfg_ana.anovamat,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,cfg_ana.printTable_tex);

%%  2-way ANOVA: Hemisphere x Condition: P O/N

cfg_ana = [];
cfg_ana.hemi = {'LPS','RPS'};
cfg_ana.cond = {'RCR','RHSC','RHSI'};
cfg_ana.latency = [0.5 0.8];
cfg_ft.frequency = [4 8];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;
if length(cfg_ana.hemi) > 2 || length(cfg_ana.cond) > 2
  cfg_ana.calcGGHF = 1;
else
  cfg_ana.calcGGHF = 0;
end

fprintf('====================== RMAOV2_mod =========================\n');
fprintf('IV1: Hemisphere (%d;%s), IV2: Condition (%d;%s)\n',length(cfg_ana.hemi),sprintf(repmat(' %s',1,length(cfg_ana.hemi)),cfg_ana.hemi{:}),length(cfg_ana.cond),sprintf(repmat(' %s',1,length(cfg_ana.cond)),cfg_ana.cond{:}));

cfg_ana.anovamat = [];
for h = 1:length(cfg_ana.hemi)
  cfg_ana.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.hemi(h))});
  cfg_ana.chan = unique(str2double(strrep(cfg_ana.channel,'e','')));
  
  for c = 1:length(cfg_ana.cond)
    cfg_ana.timesel.(cfg_ana.cond{c}) = find(ga_freq.(cfg_ana.cond{c}).time >= cfg_ana.latency(1) & ga_freq.(cfg_ana.cond{c}).time <= cfg_ana.latency(2));
    cfg_ana.freqsel.(cfg_ana.cond{c}) = find(ga_freq.(cfg_ana.cond{c}).freq >= cfg_ft.frequency(1) & ga_freq.(cfg_ana.cond{c}).freq <= cfg_ft.frequency(2));
    cfg_ana.values.(cfg_ana.cond{c}) = mean(mean(mean(ga_freq.(cfg_ana.cond{c}).powspctrm(:,cfg_ana.chan,cfg_ana.freqsel.(cfg_ana.cond{c}),cfg_ana.timesel.(cfg_ana.cond{c})),2),3),4);
    
    %index = h + (c - 1) + (h - 1);
    for sub = 1:size(cfg_ana.values.(cfg_ana.cond{c}),1)
      
      % format for RMAOV2_mod: [data hemi cond subNum]
      cfg_ana.anovamat = cat(1,cfg_ana.anovamat,[cfg_ana.cfg_ana.values.(cfg_ana.cond{c})(sub) h c sub]);
    end
  end
end
[cfg_ana.p] = RMAOV2_mod(cfg_ana.anovamat,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,cfg_ana.printTable_tex);

%% cluster statistics

cfg_ft = [];
%start and stop time of analysis, in sec
cfg_ft.latency = [0.0 1.0];
cfg_ft.avgovertime = 'no';

cfg_ft.frequency = [4 8];
cfg_ft.avgoverfreq = 'yes';

cfg_ft.elec = ga_freq.elec;
cfg_ft.channel = 'all';
%cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS'})});
%cfg_ft.avgoverchan = 'yes';

%cfg_ft.parameter = 'powspctrm';
cfg_ft.parameter = 'cohspctrm';

% use the Monte Carlo Method to calculate the significance probability
cfg_ft.method = 'montecarlo';
cfg_ft.correctm = 'cluster';
% alpha level of the sample-specific test statistic that will be used for
% thresholding
cfg_ft.clusteralpha = 0.05;
% test statistic that will be evaluated under the permutation distribution
cfg_ft.clusterstatistic = 'maxsum';
% minimum number of neighborhood channels that is required for a selected
% sample to be included in the clustering algorithm (default = 0)
cfg_ft.minnbchan = 2;
% alpha level of the permutation test
cfg_ft.alpha = 0.025;
% number of draws from the permutation distribution, should be 1000
cfg_ft.numrandomization = 1000;

% set the unit and independent variables
cfg_ft.uvar = 1; % row of design matrix containing subject numbers (DV)
cfg_ft.ivar = 2; % row of design matrix containing conditions (IV)

cfg_ana = [];
cfg_ana.numConds = 2;
% use the dependent samples T-statistic as a measure to evaluate the
% effect at the sample level
if cfg_ana.numConds == 2
  cfg_ft.statistic = 'depsamplesT';
  % test tails: -1 = left tail, 0 = two tail, 1 = right tail?
  cfg_ft.tail = 0;
  cfg_ft.clustertail = 0;
elseif cfg_ana.numConds > 2
  cfg_ft.statistic = 'depsamplesF';
  % test tails: -1 = left tail, 0 = two tail, 1 = right tail?
  cfg_ft.tail = 1;
  cfg_ft.clustertail = 1;
end

% make the design matrix
cfg_ana.numSub = length(exp.subjects) - sum(exp.badSub);
cfg_ft.design = zeros(2,cfg_ana.numSub*cfg_ana.numConds);
for i = 1:cfg_ana.numConds
  for j = 1:cfg_ana.numSub
    cfg_ft.design(1,1+((i - 1)*cfg_ana.numSub) + (j - 1)) = j; % subject #s
    cfg_ft.design(2,1+((i - 1)*cfg_ana.numSub) + (j - 1)) = i; % condition #s
  end
end

% run the nonparametric cluster statistics
stat_clus = [];
stat_clus.RHSCvsRHSI = ft_freqstatistics(cfg_ft,ga_freq.RHSC,ga_freq.RHSI);
stat_clus.RHSCvsRCR = ft_freqstatistics(cfg_ft,ga_freq.RHSC,ga_freq.RCR);
stat_clus.RHSIvsRCR = ft_freqstatistics(cfg_ft,ga_freq.RHSI,ga_freq.RCR);
stat_clus.RHvsRCR = ft_freqstatistics(cfg_ft,ga_freq.RH,ga_freq.RCR);

% save(fullfile(dirs.saveDirProc,'stat_clus.mat'),'stat_clus.RHSCvsRHSI','stat_clus.RHSCvsRCR','stat_clus.RHSCvsRCR','stat_clus.RHvsRCR');
% load(fullfile(dirs.saveDirProc,'stat_clus.mat'));

% % run the nonparametric cluster statistics
% stat.RHSCvsRHSIvsRCR = ft_freqstatistics(cfg_ft,ga_freq.RHSC,ga_freq.RHSI,ga_freq.RCR);

%% plot the cluster statistics

cfg_ft = [];
% p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
cfg_ft.highlightsymbolseries = ['*','+','.','.','.'];
cfg_ft.highlightcolorpos = [0.5 0 1];
cfg_ft.highlightcolorneg = [0 0.5 0];
cfg_ft.elec = ga_freq.elec;
cfg_ft.contournum = 0;
cfg_ft.emarker = '.';
cfg_ft.alpha  = 0.05;
cfg_ft.parameter = 'stat';
cfg_ft.zlim = [-5 5];
ft_clusterplot(cfg_ft,stat_clus.RHSCvsRHSI);
ft_clusterplot(cfg_ft,stat_clus.RHSCvsRCR);
ft_clusterplot(cfg_ft,stat_clus.RHSIvsRCR);
ft_clusterplot(cfg_ft,stat_clus.RHvsRCR);

%% Make contrast plots (with culster stat info)

% set up contrast
cfg_ana = [];
cfg_ana.include_clus_stat = 0;
cfg_ana.timeS = (0:0.05:1.0);
cfg_ana.timeSamp = round(linspace(1,exp.sampleRate,length(cfg_ana.timeS)));

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
