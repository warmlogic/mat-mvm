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
%exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
%exper.eventValues = sort({'CR2','HSC2','HSI2','CR6','HSC6','HSI6'});
%exper.eventValues = sort({'F2','F6','N2','N6','RO2','RO6','RS2','RS6'});
exper.eventValues = sort({'FSC2','FSI2','FSC6','FSI6','NM2','NM6','NS2','NS6','ROSC2','ROSI2','ROSC6','ROSI6','RSSC2','RSSI2','RSSC6','RSSI6'});

% combine some events into higher-level categories
%exper.eventValuesExtra.toCombine = {{'CR2','CR6'},{'HSC2','HSC6'},{'HSI2','HSI6'},{'HSC2','HSI2','HSC6','HSI6'}};
%exper.eventValuesExtra.newValue = {{'RCR'},{'RHSC'},{'RHSI'},{'RH'}};

%exper.eventValuesExtra.toCombine = {{'FSC2','FSC6'},{'FSI2','FSI6'},{'NM2','NM6','NS2','NS6'},{'ROSC2','ROSC6'},{'ROSI2','ROSI6'},{'RSSC2','RSSC6'},{'RSSI2','RSSI6'}};
%exper.eventValuesExtra.newValue = {{'FSC'},{'FSI'},{'N'},{'ROSC'},{'ROSI'},{'RSSC'},{'RSSI'}};

%exper.eventValuesExtra.toCombine = {{'FSC2','FSC6','FSI2','FSI6'},{'NM2','NM6','NS2','NS6'},{'ROSC2','ROSC6','ROSI2','ROSI6'},{'RSSC2','RSSC6','RSSI2','RSSI6'}};
%exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'}};

exper.eventValuesExtra.toCombine = {...
  {'NM2','NM6','NS2','NS6'},...
  {'FSC2','FSC6','ROSC2','ROSC6','RSSC2','RSSC6'},...
  {'FSI2','FSI6','ROSI2','ROSI6','RSSI2','RSSI6'}};
exper.eventValuesExtra.newValue = {{'CR'},{'SC'},{'SI'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 1;

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
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_0'};

%% set up file and directory handling parameters

% directory where the data to read is located
%dirs.subDir = 'RK';
dirs.subDir = '';
dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);

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
%   for ses = 1:length(exper.sesStr)
%     overwriteArtFields = 1;
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,overwriteArtFields);
%   end
% end

%% Convert the data to FieldTrip structs

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_freqanalysis';
ana.artifact.type = {'zeroVar'};
%ana.artifact.type = {'nsAuto'};
%ana.artifact.type = {'nsAuto','preRejManual','ftICA'};

% any preprocessing?
cfg_pp = [];
% single precision to save space
%cfg_pp.precision = 'single';
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];
% cfg_pp.detrend = 'yes';
% cfg_pp.dftfilter = 'yes';
% cfg_pp.dftfreq = [60 120 180];
% % cfg_pp.bsfilter = 'yes';
% % cfg_pp.bsfreq = [59 61; 119 121; 179 181];
% % cfg_pp.bsfreq = [58 62; 118 122; 178 182];
% % cfg_pp.lpfilter = 'yes';
% % cfg_pp.lpfreq = [35];

cfg_proc = [];
cfg_proc.pad = 'maxperlen';
% cfg_proc.output = 'pow';
% % cfg_proc.output = 'powandcsd';
% cfg_proc.keeptrials = 'yes';
% cfg_proc.keeptapers = 'no';

cfg_proc.output = 'fourier';
cfg_proc.keeptrials = 'yes';
cfg_proc.keeptapers = 'yes';

% % MTM FFT
% cfg_proc.method = 'mtmfft';
% cfg_proc.taper = 'dpss';
% %cfg_proc.foilim = [3 50];
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% %cfg_proc.foi = 3:freqstep:50;
% cfg_proc.foi = 3:freqstep:9;
% cfg_proc.tapsmofrq = 5;
% cfg_proc.toi = -0:0.04:1.0;

% % multi-taper method - Usually to up 30 Hz
% cfg_proc.method = 'mtmconvol';
% cfg_proc.taper = 'hanning';
% %cfg_proc.toi = -0.8:0.04:3.0;
% cfg_proc.toi = -0.5:0.04:1.0;
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% cfg_proc.foi = 3:freqstep:40;
% %cfg_proc.foi = 3:freqstep:9;
% %cfg_proc.foi = 3:1:9;
% %cfg_proc.foi = 2:2:30;
% % temporal smoothing
% cfg_proc.t_ftimwin = 5 ./ cfg_proc.foi;
% % frequency smoothing (tapsmofrq) is not used for hanning taper

% % multi-taper method - Usually above 30 Hz
% cfg_proc.method = 'mtmconvol';
% cfg_proc.taper = 'dpss';
% cfg_proc.toi = -0.5:0.04:1.0;
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% cfg_proc.foi = 3:freqstep:40;
% %cfg_proc.foi = 3:freqstep:9;
% %cfg_proc.foi = 3:1:9;
% %cfg_proc.foi = 2:2:30;
% % temporal smoothing
% cfg_proc.t_ftimwin = 5 ./ cfg_proc.foi;
% % frequency smoothing (tapsmofrq) is used for dpss
% cfg_proc.tapsmofrq = 0.3 .* cfg_proc.foi;

% wavelet
cfg_proc.method = 'wavelet';
cfg_proc.width = 6;
%cfg_proc.toi = -0.8:0.04:3.0;
cfg_proc.toi = -0.5:0.04:1.0;
% % evenly spaced frequencies, but not as many as foilim makes
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% % cfg_proc.foi = 3:freqstep:9;
% cfg_proc.foi = 3:freqstep:60;
cfg_proc.foi = 4:1:100;
%cfg_proc.foi = 4:1:60;
%cfg_proc.foilim = [3 9];

% log-spaced freqs
%cfg_proc.foi = (2^(1/8)).^(16:42);

% set the save directories; final argument is prefix of save directory
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,cfg_proc.output);

% ftype: a string used in naming the processed files (data_FTYPE_EVENT.mat)
ana.ftype = cfg_proc.output;

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);
process_ft_data(ana,cfg_proc,exper,dirs);

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

%adFile = '/Volumes/curranlab/Data/SOCO/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0/pow_mtmconvol_hanning_pow_-500_980_2_40_avg/analysisDetails.mat';

% fourier 4-100 Hz
adFile = '/Volumes/curranlab/Data/SOCO/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/fourier_wavelet_w6_fourier_-500_980_4_100/analysisDetails.mat';

% % wavelet
% adFile = '/Volumes/curranlab/Data/SOCO/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/pow_wavelet_w5_pow_-500_980_3_40_avg/analysisDetails.mat';

% % multitaper
% adFile = '/Volumes/curranlab/Data/SOCO/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/pow_mtmconvol_dpss_pow_-500_980_3_40_avg/analysisDetails.mat';

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

% this is only used by mm_ft_checkCondComps to create pairwise combinations
% either within event types {'all_within_types'} or across all event types
% {'all_across_types'}; mm_ft_checkCondComps is called within subsequent
% analysis functions

% list the values separated by types: 2Colors, 6Colors
ana.eventValues = {exper.eventValues};
%ana.eventValues = {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}};
%ana.eventValues = {{'F','N','RO','RS'}};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end


%% load some data

%[data_pow] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

%% new loading workflow - pow

cfg = [];
cfg.keeptrials = 'no';
%cfg.keeptrials = 'yes';
cfg.equatetrials = 'no';
%cfg.equatetrials = 'yes';
cfg.ftype = 'fourier';
cfg.output = 'pow'; % 'pow', 'coh', 'phase'
cfg.normalize = 'log10'; % 'log10', 'log', 'vector', 'dB'
%cfg.baselinetype = 'zscore'; % 'zscore', 'absolute', 'relchange', 'relative', 'condition' (use ft_freqcomparison)
cfg.baselinetype = 'absolute'; % 'zscore', 'absolute', 'relchange', 'relative', 'condition' (use ft_freqcomparison)
cfg.baseline = [-0.4 -0.2];

cfg.saveFile = true;
%cfg.saveFile = false;

cfg.rmevoked = 'yes';
cfg.rmevokedfourier = 'yes';
cfg.rmevokedpow = 'no';
if strcmp(cfg.rmevoked,'yes') && ~exist('data_evoked','var')
  load('/Volumes/curranlab/Data/SOCO/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0_art_zeroVar/tla_-1000_2000_avg/data_evoked.mat');
end

if isfield(cfg,'equatetrials') && strcmp(cfg.equatetrials,'yes')
  eq_str = '_eq';
else
  eq_str = '';
end
if isfield(cfg,'keeptrials') && strcmp(cfg.keeptrials,'yes')
  kt_str = '_trials';
else
  kt_str = '_avg';
end
if isfield(cfg,'rmevoked') && strcmp(cfg.rmevoked,'yes')
  indu_str = '_induced';
else
  indu_str = '';
end
saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s%s%s%s.mat',cfg.output,eq_str,kt_str,indu_str));

if exist(saveFile,'file')
  fprintf('Loading saved file: %s\n',saveFile);
  load(saveFile);
else
  fprintf('Running mm_ft_loadData\n');
  [data_pow,exper] = mm_ft_loadData(cfg,exper,dirs,ana,data_evoked);
  if cfg.saveFile
    fprintf('Saving %s...\n',saveFile);
    save(saveFile,sprintf('data_%s',cfg.output),'exper','cfg');
  end
end
fprintf('Done.\n');

%% new loading workflow - coherence

cfg = [];
cfg.keeptrials = 'no';
cfg.equatetrials = 'no';
%cfg.equatetrials = 'yes';
cfg.ftype = 'fourier';
cfg.output = 'coh'; % 'pow', 'coh', 'phase'
cfg.baselinetype = 'absolute'; % 'absolute', 'relchange', 'relative', 'condition' (use ft_freqcomparison)
cfg.baseline = [-0.4 -0.2];

cfg.saveFile = true;
%cfg.saveFile = false;

if strcmp(cfg.equatetrials,'yes')
  eq_str = '_eq';
elseif strcmp(cfg.equatetrials,'no')
  eq_str = '';
end
if strcmp(cfg.keeptrials,'yes')
  kt_str = '_trials';
elseif strcmp(cfg.keeptrials,'no')
  kt_str = '_avg';
end
saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s%s%s.mat',cfg.output,eq_str,kt_str));

if exist(saveFile,'file')
  fprintf('Loading saved file: %s\n',saveFile);
  load(saveFile);
elseif ~exist(saveFile,'file')
  fprintf('Running mm_ft_loadData\n');
  [data_coh,exper] = mm_ft_loadData(cfg,exper,dirs,ana);
  if cfg.saveFile
    fprintf('Saving %s...\n',saveFile);
    save(saveFile,sprintf('data_%s',cfg.output),'exper','cfg');
  end
end
fprintf('Done.\n');

%% new loading workflow - phase

cfg = [];
cfg.keeptrials = 'yes';
cfg.equatetrials = 'no';
%cfg.equatetrials = 'yes';
cfg.ftype = 'fourier';
cfg.output = 'phase'; % 'pow', 'coh', 'phase'
%cfg.baselinetype = 'absolute'; % 'absolute', 'relchange', 'relative', 'condition' (use ft_freqcomparison)
%cfg.baseline = [-0.4 -0.2];

cfg.phasefreq = [4 8; 8 12; 12 28; 28 50; 50 100];
%cfg.phaseroi = {{'LPS'}};
cfg.phaseroi = {{'E11'},{'E62'}};

%cfg.saveFile = true;
cfg.saveFile = false;

if strcmp(cfg.equatetrials,'yes')
  eq_str = '_eq';
elseif strcmp(cfg.equatetrials,'no')
  eq_str = '';
end
if strcmp(cfg.keeptrials,'yes')
  kt_str = '_trials';
elseif strcmp(cfg.keeptrials,'no')
  kt_str = '_avg';
end
saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s%s%s.mat',cfg.output,eq_str,kt_str));

if exist(saveFile,'file')
  fprintf('Loading saved file: %s\n',saveFile);
  load(saveFile);
elseif ~exist(saveFile,'file')
  fprintf('Running mm_ft_loadData\n');
  [data_phase,exper] = mm_ft_loadData(cfg,exper,dirs,ana);
  if cfg.saveFile
    fprintf('Saving %s...\n',saveFile);
    save(saveFile,sprintf('data_%s',cfg.output),'exper','cfg');
  end
end
fprintf('Done.\n');

% TODO: mm_ft_loadData runs: mm_ft_freqnormalize, mm_ft_freqbaseline

%% plot something

sub=2;
ses=1;

%chan=11; % Fz
%chan=62; % Pz
%chan=20; % LAS middle
%chan=118; % RAS middle
chan=53; % LPS middle
%chan=86; % RPS middle

cond = {'CR','SC','SI'};

%if strcmp(cfg.output,'pow')
param = 'powspctrm';
clim = [-1 1];
for cnd = 1:length(cond)
  if isfield(data_pow.(cond{cnd}).sub(sub).ses(ses).data,param)
    figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(data_pow.(cond{cnd}).sub(sub).ses(ses).data.(param)(chan,:,:)),clim);
    axis xy;colorbar;
    title(sprintf('Z-Power: %s, sub %d, ses %d, chan %d',cond{cnd},sub,ses,chan));
  end
end
%elseif strcmp(cfg.output,'coh')
param = 'powspctrm';
%clim = [0 1];
clim = [0 0.5];
for cnd = 1:length(cond)
  if isfield(data_coh.(cond{cnd}).sub(sub).ses(ses).data,param)
    figure;imagesc(data_coh.(cond{cnd}).sub(sub).ses(ses).data.time,data_coh.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(data_coh.(cond{cnd}).sub(sub).ses(ses).data.(param)(chan,:,:)),clim);
    axis xy;colorbar;
    title(sprintf('ITC - BL: %s, sub %d, ses %d, chan %d',cond{cnd},sub,ses,chan));
  end
end
%end

%% Test plots to make sure data look ok

cfg_ft = [];
cfg_ft.baseline = [-0.3 -0.1];
cfg_ft.baselinetype = 'absolute'; % maybe this
%cfg_ft.baselinetype = 'relative';
%cfg_ft.baselinetype = 'relchange'; % or this
if strcmp(cfg_ft.baselinetype,'absolute')
  cfg_ft.zlim = [-400 400];
  %cfg_ft.zlim = [-2 2];
elseif strcmp(cfg_ft.baselinetype,'relative')
  cfg_ft.zlim = [0 2.0];
elseif strcmp(cfg_ft.baselinetype,'relchange')
  cfg_ft.zlim = [-1.0 1.0];
end
cfg_ft.parameter = 'powspctrm';
%cfg_ft.ylim = [3 9];
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
sub=1;
ses=1;
for i = 1:length(ana.eventValues{1})
  figure
  ft_multiplotTFR(cfg_ft,data_pow.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
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
% ft_singleplotTFR(cfg_ft,data_pow.(exper.eventValues{1}).sub(1).ses(1).data);

%% Change in freq relative to baseline using absolute power

cfg_fb = [];
cfg_fb.baseline = [-0.3 -0.1];
cfg_fb.baselinetype = 'absolute'; % maybe this
%cfg_fb.baselinetype = 'relative';
%cfg_fb.baselinetype = 'relchange'; % or this

%data_pow_orig = data_pow;

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        fprintf('%s, %s, %s, ',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
        data_pow.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_pow.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
      end
    end
  end
end

% % find the time points without NaNs for a particular frequency
% data_pow.(exper.eventValues{1}).sub(1).ses(1).data.time(~isnan(squeeze(data_pow.(exper.eventValues{1}).sub(1).ses(1).data.powspctrm(1,2,:))'))
% ga_pow.(exper.eventValues{1}).time(~isnan(squeeze(ga_pow.(exper.eventValues{1}).powspctrm(1,1,2,:))'))

%% decide who to kick out based on trial counts

% Subjects with bad behavior
%exper.badBehSub = {};
% huge response bias to say "new": 18, 26
exper.badBehSub = {'SOCO018','SOCO026'}; % for publication

% exclude subjects with low event counts
[exper] = mm_threshSubs(exper,ana,15);

%% for evoked TF, get rid of the trial dim

data_evoked_orig = data_evoked;

cfg_fd = [];

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        if isfield(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,param)
          %fprintf('%s, %s, %s, ',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
          
          % % phase - need to deal with in a different way, can't average
          % data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.powspctrm = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.phasespctrm;
          % data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'phasespctrm');
          % data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'fourierspctrm');
          
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg_fd,data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
        end
      end
    end
  end
end

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = ana.eventValues;
cfg_ana.data_str = 'data_pow';
%cfg_ana.data_str = 'data_coh';
%cfg_ana.data_str = 'data_evoked';
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sesStr)
  for typ = 1:length(ana.eventValues)
    for evVal = 1:length(ana.eventValues{typ})
      %tic
      fprintf('Running ft_freqgrandaverage on %s...',ana.eventValues{typ}{evVal});
      if strcmp(cfg_ana.data_str,'data_pow')
        cfg_ft.parameter = 'powspctrm';
        ga_pow.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
      elseif strcmp(cfg_ana.data_str,'data_coh')
        %cfg_ft.parameter = 'plvspctrm';
        cfg_ft.parameter = 'powspctrm';
        ga_coh.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
      elseif strcmp(cfg_ana.data_str,'data_evoked')
        cfg_ft.parameter = 'powspctrm';
        ga_evoked.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
      end
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
%cfg_ft.xlim = [0 1];
%cfg_ft.ylim = [4 8];
cfg_ft.ylim = [4 100];
%cfg_ft.ylim = [8 12];
%cfg_ft.zlim = [-1 1];
%cfg_ft.zlim = [-2 2];
%cfg_ft.zlim = [-400 400];
cfg_ft.zlim = [-0.2 0.2];
%elseif strcmp(cfg_ft.baselinetype,'relative')
%  cfg_ft.zlim = [0 2.0];
%end

% if strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'absolute')
%   cfg_ft.zlim = [-400 400];
% elseif strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'relative')
%   cfg_ft.zlim = [0 2.0];
% elseif strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'relchange')
%   cfg_ft.zlim = [-1.0 1.0];
% end

cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
for typ = 1:length(ana.eventValues)
  for evVal = 1:length(ana.eventValues{typ})
    figure
    ft_multiplotTFR(cfg_ft,ga_pow.(ana.eventValues{typ}{evVal}));
    set(gcf,'Name',sprintf('%s',ana.eventValues{typ}{evVal}))
  end
end

%% subplots of each subject's power spectrum

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.rois = {{'FS'},{'PS'}};
%cfg_plot.rois = {'E17','E62'};
%cfg_plot.roi = {'RAS'};
%cfg_plot.roi = {'LPS','RPS'};
%cfg_plot.roi = {'LPS'};
cfg_plot.excludeBadSub = 1;
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
cfg_ft.parameter = 'powspctrm';

if strcmp(ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype'),'absolute')
  cfg_ft.zlim = [-400 400];
elseif strcmp(ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype'),'relative')
  cfg_ft.zlim = [0 2.0];
elseif strcmp(ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype'),'relchange')
  cfg_ft.zlim = [-1.0 1.0];
end

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotTFR(cfg_ft,cfg_plot,ana,exper,data_pow);
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
%cfg_ft.zlim = [-100 100]; % pow

if strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'absolute')
  cfg_ft.zlim = [-400 400];
elseif strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'relative')
  cfg_ft.zlim = [0 2.0];
elseif strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'relchange')
  cfg_ft.zlim = [-1.0 1.0];
end

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
cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'RCR','RH','RHSC','RHSI'}},size(cfg_plot.rois));

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

%cfg_plot.ftFxn = 'ft_singleplotTFR';

% cfg_plot.ftFxn = 'ft_topoplotTFR';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.comment = 'no';
% %cfg_ft.xlim = [0.5 0.8]; % time
% cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time

cfg_plot.ftFxn = 'ft_multiplotTFR';
cfg_ft.showlabels = 'yes';
cfg_ft.comment = '';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_pow);
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
cfg_ft.parameter = 'powspctrm';

if strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'absolute')
  cfg_ft.zlim = [-400 400];
elseif strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'relative')
  cfg_ft.zlim = [0 2.0];
elseif strcmp(ft_findcfg(ga_pow.(ana.eventValues{1}{1}).cfg,'baselinetype'),'relchange')
  cfg_ft.zlim = [-1.0 1.0];
end

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

mm_ft_contrastTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_pow);

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
%cfg_plot.plot_order = {'CR2','H2','HSC2','HSI2','CR6','H6','HSC6','HSI6'};

%cfg_plot.ylims = repmat([-1 1],size(cfg_ana.rois'));
%cfg_plot.ylims = repmat([-100 100],size(cfg_ana.rois'));
if strcmp(ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype'),'absolute')
  cfg_plot.ylims = repmat([-100 100],size(cfg_ana.rois'));
elseif strcmp(ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype'),'relative')
  cfg_plot.ylims = repmat([0 2.0],size(cfg_ana.rois'));
elseif strcmp(ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype'),'relchange')
  cfg_plot.ylims = repmat([-1.0 1.0],size(cfg_ana.rois'));
else
  cfg_plot.ylims = repmat([-2.0 2.0],size(cfg_ana.rois'));
end

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_pow);
  %mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_coh);
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
  
  mm_ft_rmaov2TFR(cfg_ana,exper,ana,data_pow);
end

%% cluster statistics

cfg_ft = [];
cfg_ft.avgoverchan = 'no';
%cfg_ft.avgovertime = 'no';
cfg_ft.avgovertime = 'yes';
%cfg_ft.avgoverfreq = 'no';
cfg_ft.avgoverfreq = 'yes';

cfg_ft.parameter = 'powspctrm';

% debugging
%cfg_ft.numrandomization = 10;

cfg_ft.numrandomization = 500;
% cfg_ft.clusteralpha = .05;
% cfg_ft.alpha = .025;
cfg_ft.clusteralpha = .1;
cfg_ft.alpha = .05;

cfg_ana = [];
cfg_ana.roi = 'all';
cfg_ana.avgFrq = cfg_ft.avgoverfreq;
cfg_ana.conditions = {'all'};
%cfg_ana.conditions = {'all_within_types'};
% cfg_ana.conditions = {...
%   {'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},...
%   {'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'},...
%   {'CR2','CR6'},{'H2','H6'},{'HSC2','HSC6'},{'HSI2','HSI6'}};

% extra identifier when saving
%thisBL = ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baseline');
%cfg_ana.dirStr = sprintf('_%s_%d_%d',ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype'),thisBL(1)*1000,thisBL(2)*1000);

thisBLtype = 'pow_induced';
%thisBLtype = 'zpow_indu';
%thisBLtype = 'zp_evok';
%thisBLtype = 'coh';
thisBL = [-0.4 -0.2];
cfg_ana.dirStr = sprintf('_%s_%d_%d',thisBLtype,thisBL(1)*1000,thisBL(2)*1000);

if strcmp(cfg_ft.avgovertime,'no')
  cfg_ana.latencies = [0 1.0];
  %cfg_ana.latencies = [0 0.5; 0.5 1.0];
elseif strcmp(cfg_ft.avgovertime,'yes')
  %cfg_ana.latencies = [-0.2 -0.1; -0.1 0; 0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
  %cfg_ana.latencies = [-0.2 0; 0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
  %cfg_ana.latencies = [-0.2:0.1:0.9; -0.1:0.1:1.0]';
  cfg_ana.latencies = [-0.2:0.2:0.8; 0:0.2:1.0]';
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end

if strcmp(cfg_ft.avgoverfreq,'no')
  cfg_ana.frequencies = [4 100];
elseif strcmp(cfg_ft.avgoverfreq,'yes')
  cfg_ana.frequencies = [4 8; 8 12; 12 28; 28 50; 50 100];
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgF'];
end

if strcmp(cfg_ft.avgoverchan,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgC'];
end

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  for fr = 1:size(cfg_ana.frequencies,1)
    cfg_ft.frequency = cfg_ana.frequencies(fr,:);
    
    if ~isempty(strfind(cfg_ana.dirStr,'pow'))
      [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_pow);
    elseif ~isempty(strfind(cfg_ana.dirStr,'coh'))
      [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_coh);
    elseif ~isempty(strfind(cfg_ana.dirStr,'evok'))
      [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_evoked);
    end
  end
end

%% plot the cluster statistics

files.saveFigs = 1;

cfg_ft = [];
%cfg_ft.alpha = .025;
%cfg_ft.alpha = .05;
cfg_ft.alpha = .1;
cfg_ft.avgoverfreq = cfg_ana.avgFrq;

cfg_plot = [];
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.frequencies = cfg_ana.frequencies;
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.dirStr = cfg_ana.dirStr;

if strcmp(cfg_ft.avgoverfreq,'no')
  % not averaging over frequencies - only works with ft_multiplotTFR
  files.saveFigs = 0;
  cfg_ft.interactive = 'yes';
  cfg_plot.mask = 'yes';
  %cfg_ft.maskstyle = 'saturation';
  cfg_ft.maskalpha = 0.3;
  cfg_plot.ftFxn = 'ft_multiplotTFR';
  % http://mailman.science.ru.nl/pipermail/fieldtrip/2009-July/002288.html
  % http://mailman.science.ru.nl/pipermail/fieldtrip/2010-November/003312.html
end

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
  subject = sprintf('Done with %s cluster pow:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with %s cluster pow:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    };
  send_gmail(subject,mail_message);
end

%% line plots

files.saveFigs = 0;

cfg = [];
cfg.parameter = 'powspctrm';

%cfg.times = [-0.2:0.05:0.9; -0.1:0.05:1.0]';
%cfg.times = [-0.2:0.1:0.9; -0.1:0.1:1.0]';
cfg.times = [-0.2:0.2:0.8; 0:0.2:1.0]';

cfg.freqs = [4 8; 8 12; 12 28; 28 50; 50 100];
%cfg.freqs = [4 8];

cfg.rois = {...
  {'LAS'},{'FS'},{'RAS'},...
  {'LPS'},{'PS'},{'RPS'},...
  };

% cfg.rois = {...
%   {'LAI'},{'FI'},{'RAI'},...
%   {'LAS'},{'FS'},{'RAS'},...
%   {'LPS'},{'PS'},{'RPS'},...
%   {'LPI'},{'PI'},{'RPI'},...
%   };

cfg.conditions = ana.eventValues;

cfg.plotTitle = true;
cfg.plotLegend = true;

cfg.plotClusSig = true;
cfg.clusAlpha = 0.1;
%cfg.clusTimes = cfg.times;
cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
cfg.clusLimits = true;

%cfg.ylim = [-0.6 0.6];
%cfg.ylim = [-0.5 0.2];
cfg.nCol = 3;

% % whole power
% cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_pow);

% induced power
cfg.type = 'line_pow_induced';
cfg.clusDirStr = '_pow_induced_-400_-200';
cfg.ylabel = 'Log Pow';
mm_ft_lineTFR(cfg,ana,files,dirs,ga_pow);

% % induced power
% cfg.type = 'line_pow_indu';
% cfg.clusDirStr = '_zpow_indu_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_pow);

% % evoked power
% cfg.type = 'line_pow_evok';
% cfg.clusDirStr = '_zp_evok_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_evoked);

% cfg.type = 'line_coh';
% cfg.clusDirStr = '_coh_-400_-200';
% cfg.ylabel = 'ITC';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_coh);

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
