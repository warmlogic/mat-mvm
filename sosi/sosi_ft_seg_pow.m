% Make plots and do analyses for time-frequency EEG data

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SOSI';

exper.sampleRate = 250;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 2.0];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
%exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
%exper.eventValues = sort({'RCR','RHSC','RHSI'});
exper.eventValues = sort({'FSC','FSI','NM','NS','ROSC','ROSI','RSSC','RSSI'});

% combine the two types of hits into one category
%exper.eventValuesExtra.toCombine = {{'RHSC','RHSI'}};
%exper.eventValuesExtra.newValue = {{'RH'}};

%exper.eventValuesExtra.toCombine = {{'FSC','FSI'},{'NS','NM'},{'ROSC','ROSI'},{'RSSC','RSSI'}};
%exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'}};

exper.eventValuesExtra.toCombine = {{'NS','NM'},{'FSC','ROSC','RSSC'},{'FSI','ROSI','RSSI'}};
exper.eventValuesExtra.newValue = {{'CR'},{'SC'},{'SI'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 1;

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

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_0'};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));
%dirs.dataDir = fullfile(exper.name,'eeg','ftpp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));

% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
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
%   for ses = 1:length(exper.sesStr)
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,0);
%   end
% end

%% Convert the data to FieldTrip structs

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_freqanalysis';
ana.artifact.type = {'zeroVar'};
%ana.artifact.type = {'nsAuto'};
%ana.artifact.type = {'nsAuto','preRejManual','ftICA'};

% ana.otherFxn = {};
% ana.cfg_other = [];
% ana.otherFxn{1} = 'ft_preprocessing';
% ana.cfg_other{1}.demean = 'yes';
% ana.cfg_other{1}.baselinewindow = [-.2 0];
% ana.cfg_other{1}.ftype = 'demean';
% ana.otherFxn{2} = 'ft_preprocessing';
% ana.cfg_other{2}.detrend = 'yes';
% ana.cfg_other{2}.ftype = 'detrend';
% ana.otherFxn{3} = 'ft_preprocessing';
% ana.cfg_other{3}.dftfilter = 'yes';
% ana.cfg_other{3}.dftfreq = [60 120 180];
% ana.cfg_other{3}.ftype = 'dft';
% % ana.otherFxn{4} = 'ft_preprocessing';
% % ana.cfg_other{4}.bsfilter = 'yes';
% % ana.cfg_other{4}.bsfreq = [59 61; 119 121; 179 181];
% % ana.cfg_other{4}.ftype = 'bs';
% % ana.otherFxn{5} = 'ft_preprocessing';
% % ana.cfg_other{5}.lpfilter = 'yes';
% % ana.cfg_other{5}.lpfreq = [35];
% % ana.cfg_other{5}.ftype = 'lp';

% any preprocessing?
cfg_pp = [];
% single precision to save space
%cfg_pp.precision = 'single';
% baseline correct
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
% cfg_proc.foilim = [4 8];
% % freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% % cfg_proc.foi = 3:freqstep:40;
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
%cfg_proc.foi = 4:1:30;
%cfg_proc.foilim = [3 9];

% log-spaced freqs
%cfg_proc.foi = (2^(1/8)).^(16:42);

% set the save directories; final argument is prefix of save directory
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,cfg_proc.output);

% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
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

%adFile = saveFile;

% local, testing 1 subject; fourier
adFile = '/Users/matt/data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/fourier_wavelet_w6_fourier_-500_980_4_30/analysisDetails.mat';

% % fourier 4-100 Hz
% adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/fourier_wavelet_w6_fourier_-500_980_4_100/analysisDetails.mat';

%adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0/pow_mtmconvol_hanning_pow_-500_980_2_40_avg/analysisDetails.mat';

% wavelet6 [3 9]
%adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/pow_wavelet_w6_pow_-500_980_3_9_avg/analysisDetails.mat';

% wavelet5 [3 40]
%adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/pow_wavelet_w5_pow_-500_980_3_40_avg/analysisDetails.mat';

% % multitaper
% adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar/pow_mtmconvol_dpss_pow_-500_980_3_40_avg/analysisDetails.mat';

[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

%% set up channel groups

ana = mm_ft_elecGroups(ana);

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

%% load some data

%[data_pow] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

[data_pow] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,cfg_proc.output);

%% testing ft_freqdescriptives

% can it convert from fourier to power? yes

conds = {'SC','SI','CR'};
sub=1;
ses=1;
chan=55;
%cond = 'SC';

cfg = [];
cfg.variance = 'no';
cfg.jackknife = 'no';
cfg.keeptrials = 'yes';

% calculate power
for cnd = 1:length(conds)
  data_pow.(conds{cnd}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg,data_fourier.(conds{cnd}).sub(sub).ses(ses).data);
end
cnd=1;
figure;
surf(data_pow.(conds{cnd}).sub(sub).ses(ses).data.time,data_pow.(conds{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(data_pow.(conds{cnd}).sub(sub).ses(ses).data.powspctrm(:,chan,:,:),1)));
shading interp;view([0,90]);axis tight;
title(sprintf('%s: straightup power',conds{cnd}));

% log power
data_pow_log = data_pow;
for cnd = 1:length(conds)
  data_pow_log.(conds{cnd}).sub(sub).ses(ses).data.powspctrm = log10(data_pow.(conds{cnd}).sub(sub).ses(ses).data.powspctrm);
end
cnd=1;
figure;
surf(data_pow_log.(conds{cnd}).sub(sub).ses(ses).data.time,data_pow_log.(conds{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(data_pow_log.(conds{cnd}).sub(sub).ses(ses).data.powspctrm(:,chan,:,:),1)));
shading interp;view([0,90]);axis tight;
title(sprintf('%s: straightup log power',conds{cnd}));

% baseline normalize
cfg = [];
cfg.baseline = [-0.3 -0.1];
cfg.baselinetype = 'absolute';
for cnd = 1:length(conds)
  data_pow.(conds{cnd}).sub(sub).ses(ses).data = ft_freqbaseline(cfg,data_pow.(conds{cnd}).sub(sub).ses(ses).data);
end
cnd=1;
figure;
surf(data_pow.(conds{cnd}).sub(sub).ses(ses).data.time,data_pow.(conds{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(data_pow.(conds{cnd}).sub(sub).ses(ses).data.powspctrm(:,chan,:,:),1)));
shading interp;view([0,90]);axis tight;
title(sprintf('%s: baseline normalized power',conds{cnd}));

% cfg = [];
% cfg.channel = 'E55';
% cfg.latency = [0.3 0.5];
% cfg.frequency = [8 12];
% cfg.avgovertime = 'yes';
% cfg.avgoverfreq = 'yes';
% cfg.parameter = 'powspctrm';
% stat = ft_freqstatistics(cfg,data_pow.(conds{1}).sub(sub).ses(ses).data,data_pow.(conds{2}).sub(sub).ses(ses).data);

cfg = [];
cfg.latency = [0.3 0.5];
cfg.frequency = [8 14];
timesel = data_pow.(conds{1}).sub(sub).ses(ses).data.time >= cfg.latency(1) & data_pow.(conds{1}).sub(sub).ses(ses).data.time <= cfg.latency(2);
freqsel = data_pow.(conds{1}).sub(sub).ses(ses).data.freq >= cfg.frequency(1) & data_pow.(conds{1}).sub(sub).ses(ses).data.freq <= cfg.frequency(2);
cond1 = squeeze(mean(mean(data_pow.(conds{1}).sub(sub).ses(ses).data.powspctrm(:,chan,freqsel,timesel),4),3));
cond2 = squeeze(mean(mean(data_pow.(conds{2}).sub(sub).ses(ses).data.powspctrm(:,chan,freqsel,timesel),4),3));
cond3 = squeeze(mean(mean(data_pow.(conds{3}).sub(sub).ses(ses).data.powspctrm(:,chan,freqsel,timesel),4),3));
[h,p,ci,stats] = ttest2(cond1,cond2,0.05,'both');
fprintf('cond1 vs cond2: t(%d)=%.3f, p=%.4f\n',stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest2(cond1,cond3,0.05,'both');
fprintf('cond1 vs cond3: t(%d)=%.3f, p=%.4f\n',stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest2(cond2,cond3,0.05,'both');
fprintf('cond2 vs cond3: t(%d)=%.3f, p=%.4f\n',stats.df,stats.tstat,p);

% log: baseline normalize
cfg = [];
cfg.baseline = [-0.3 -0.1];
cfg.baselinetype = 'absolute';
for cnd = 1:length(conds)
  data_pow_log.(conds{cnd}).sub(sub).ses(ses).data = ft_freqbaseline(cfg,data_pow_log.(conds{cnd}).sub(sub).ses(ses).data);
end
cnd=1;
figure;
surf(data_pow_log.(conds{cnd}).sub(sub).ses(ses).data.time,data_pow_log.(conds{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(data_pow_log.(conds{cnd}).sub(sub).ses(ses).data.powspctrm(:,chan,:,:),1)));
shading interp;view([0,90]);axis tight;
title(sprintf('%s: baseline normalized log power',conds{cnd}));

cfg = [];
cfg.latency = [0.3 0.5];
cfg.frequency = [8 14];
chansel = 55;
timesel = data_pow_log.(conds{1}).sub(sub).ses(ses).data.time >= cfg.latency(1) & data_pow_log.(conds{1}).sub(sub).ses(ses).data.time <= cfg.latency(2);
freqsel = data_pow_log.(conds{1}).sub(sub).ses(ses).data.freq >= cfg.frequency(1) & data_pow_log.(conds{1}).sub(sub).ses(ses).data.freq <= cfg.frequency(2);
cond1 = squeeze(mean(mean(data_pow_log.(conds{1}).sub(sub).ses(ses).data.powspctrm(:,chansel,freqsel,timesel),4),3));
cond2 = squeeze(mean(mean(data_pow_log.(conds{2}).sub(sub).ses(ses).data.powspctrm(:,chansel,freqsel,timesel),4),3));
cond3 = squeeze(mean(mean(data_pow_log.(conds{3}).sub(sub).ses(ses).data.powspctrm(:,chansel,freqsel,timesel),4),3));
[h,p,ci,stats] = ttest2(cond1,cond2,0.05,'both');
fprintf('log cond1 vs cond2: t(%d)=%.3f, p=%.4f\n',stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest2(cond1,cond3,0.05,'both');
fprintf('log cond1 vs cond3: t(%d)=%.3f, p=%.4f\n',stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest2(cond2,cond3,0.05,'both');
fprintf('log cond2 vs cond3: t(%d)=%.3f, p=%.4f\n',stats.df,stats.tstat,p);








% figure;
% surf(data_pow.(cond).sub(sub).ses(ses).data.time,data_pow.(cond).sub(sub).ses(ses).data.freq,squeeze(mean(data_pow.(cond).sub(sub).ses(ses).data.powspctrm(:,chan,:,:),1)));
% shading interp;view([0,90]);axis tight;
% title('straightup power');
% 
% cfg = [];
% cfg.baseline = [-0.3 -0.1];
% cfg.baselinetype = 'absolute';
% data_pow_new_bl = ft_freqbaseline(cfg,data_pow_new);
% figure;
% surf(data_pow.(cond).sub(sub).ses(ses).data.time,data_pow.(cond).sub(sub).ses(ses).data.freq,squeeze(data_pow_new_bl.powspctrm(chan,:,:)));
% shading interp;view([0,90]);axis tight;
% title('baseline corrected');
% 
% pow_norm_11 = qnorm(data_pow_new_bl.powspctrm,[1,1],'vec');
% pow_norm_21 = qnorm(data_pow_new_bl.powspctrm,[2,1],'vec');
% pow_norm_22 = qnorm(data_pow_new_bl.powspctrm,[2,2],'vec');
% pow_norm = pow_norm_22;
% figure;
% surf(data_pow.(cond).sub(sub).ses(ses).data.time,data_pow.(cond).sub(sub).ses(ses).data.freq,squeeze(pow_norm(chan,:,:)));
% shading interp;view([0,90]);axis tight;
% title('norm 22');

%% new loading workflow - pow

cfg = [];
cfg.keeptrials = 'no';
%cfg.keeptrials = 'yes';
cfg.equatetrials = 'no';
%cfg.equatetrials = 'yes';

% type of input (used in the filename to load)
cfg.ftype = 'fourier';

% type of output: 'pow', 'coh', 'phase'
cfg.output = 'pow';

% normalization type: 'log10', 'log', 'vec', 'dB'
cfg.normalize = 'log10'; 
%cfg.normalize = 'dB';
%cfg.normalize = 'vec';

% baseline type
% % 'zscore', 'absolute', 'relchange', 'relative', 'condition' (use ft_freqcomparison)
%cfg.baselinetype = 'zscore';
cfg.baselinetype = 'absolute';
% cfg.baselinetype = 'relchange';
% cfg.baselinetype = 'relative';

% baseline period
cfg.baseline = [-0.4 -0.2];
%cfg.baseline = [-0.2 0];

% at what data stage should it be baseline corrected?
cfg.baselinedata = 'mod';
% cfg.baselinedata = 'pow';


%cfg.saveFile = true;
cfg.saveFile = false;

% only keep induced data by removing evoked?
cfg.rmevoked = 'no';
cfg.rmevokedfourier = 'no';
cfg.rmevokedpow = 'no';
% cfg.rmevoked = 'yes';
% cfg.rmevokedfourier = 'yes';
% cfg.rmevokedpow = 'no';
if strcmp(cfg.rmevoked,'yes') && ~exist('data_evoked','var')
  load('/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0_art_zeroVar/tla_-1000_2000_avg/data_evoked.mat');
  
  % local testing
  %load('/Users/matt/data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar_badChanManual_badChanEP/tla_-1000_2000_avg/data_evoked.mat');
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
  indu_str = '_whole';
end
saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s%s%s%s.mat',cfg.output,eq_str,kt_str,indu_str));

if exist(saveFile,'file')
  fprintf('Loading saved file: %s\n',saveFile);
  load(saveFile);
else
  fprintf('Running mm_ft_loadData\n');
  if exist('data_evoked','var')
    [data_pow,exper] = mm_ft_loadData(cfg,exper,dirs,ana,data_evoked);
  else
    [data_pow,exper] = mm_ft_loadData(cfg,exper,dirs,ana);
  end
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

% %% rename plv to pow
% 
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sesStr)
%     for typ = 1:length(ana.eventValues)
%       for evVal = 1:length(ana.eventValues{typ})
%         if isfield(data_phase.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'plvspctrm')
%           fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
%           data_phase.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.powspctrm = data_phase.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.plvspctrm;
%           data_phase.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_phase.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'plvspctrm');
%         end
%       end
%     end
%   end
% end

%% plot something

sub=1;
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
clim = [0 1];
for cnd = 1:length(cond)
  if isfield(data_coh.(cond{cnd}).sub(sub).ses(ses).data,param)
    figure;imagesc(data_coh.(cond{cnd}).sub(sub).ses(ses).data.time,data_coh.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(data_coh.(cond{cnd}).sub(sub).ses(ses).data.(param)(chan,:,:)),clim);
    axis xy;colorbar;
    title(sprintf('ITC - BL: %s, sub %d, ses %d, chan %d',cond{cnd},sub,ses,chan));
  end
end
%end

cond = 'SI';
figure;
surf(data_pow.(cond).sub(1).ses(1).data.time,1:size(data_pow.(cond).sub(1).ses(1).data.powspctrm,1),squeeze(data_pow.(cond).sub(1).ses(1).data.powspctrm(:,55,3,:)));
shading interp;view([0,90]);axis tight;

%% playing with ft_connectivityanalysis

sub=1;
ses=1;
cond = {'SC'};
cnd = 1;

cfg = [];
%cfg.method = 'plv';
cfg.method = 'wpli_debiased';
cfg.channel = {'E11','E55'};

% this tries to run univariate2bivariate
data_conn = ft_connectivityanalysis(cfg,data_pow.(cond{cnd}).sub(sub).ses(ses).data);
param = sprintf('%sspctrm',cfg.method);
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(data_conn.(param)(1,2,:,:)));axis xy;colorbar;
title(sprintf('%s: %s-%s',cfg.method,cfg.channel{1},cfg.channel{2}));


%% playing with log/baseline correcting

% samples to use for baseline
blt = 6:11;

sub=1;
ses=1;
cond = {'SC'};
cnd = 1;

chan = 11;

% get power
%p = data_pow.(cond{cnd}).sub(sub).ses(ses).data.powspctrm;
% p(p == 0) = eps(0);

% get complex fourier spectra
f = data_pow.(cond{cnd}).sub(sub).ses(ses).data.fourierspctrm;

% turn fourier into phase-locking values (inter-trial coherence?)
pl = f ./ abs(f);
pl_itc = abs(squeeze(mean(pl,1)));
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(pl_itc(chan,:,:)));axis xy;colorbar;
title('raw phase');

% phase: subtract baseline
%
% Doesn't matter which one we use, results are the same
%
% average phase over trials, then over time
pl_blm = repmat(abs(squeeze(mean(mean(pl(:,:,:,blt),1),4))),[1,1,38]);
% average phase over time, then over trials
%pl_blm = repmat(abs(squeeze(mean(mean(pl(:,:,:,blt),4),1))),[1,1,38]);
pl_itc_abs = pl_itc - pl_blm;
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(pl_itc_abs(chan,:,:)));axis xy;colorbar;
title('raw phase - bl');

% turn fourier into power
p = abs(f).^2;
p(p == 0) = eps(0);

% raw power
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(p(:,chan,:,:),1)));axis xy;colorbar;
title('raw power');

% log10 normalized power
plog = log10(p);
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(plog(:,chan,:,:),1)));axis xy;colorbar;
title('log10(raw power)');

% z-transform power relative to baseline period
p(p == 0) = eps(0);
plog = log10(p);
plog_blm = repmat(mean(plog(:,:,:,blt),4),[1,1,1,38]);
plog_blstd = repmat(mean(std(plog(:,:,:,blt),0,1),4),[85,1,1,38]);
plog_z = (plog - plog_blm) ./ plog_blstd;
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(plog_z(:,chan,:,:),1)),[-1 1]);axis xy;colorbar;
title('z(log10(raw power) ~blm)');

% absolute: subtract baseline
plog_abs = plog - plog_blm;
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(plog_abs(:,chan,:,:),1)),[-1 1]);axis xy;colorbar;
title('absolute: log10(raw power) - blm');

% relative: divide by baseline - NOT WORKING
plog_rel = plog ./ plog_blm;
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(plog_rel(:,chan,:,:),1)));axis xy;colorbar;
title('relative: log10(raw power) ./ blm');

% relchange: subtract and by baseline (center at 0) - NOT WORKING
plog_relchg = (plog - plog_blm) ./ plog_blm;
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(plog_relchg(:,chan,:,:),1)));axis xy;colorbar;
title('relchange: (log10(raw power) - blm) ./ blm');

% % robert (get rid of trials) - no different from z-transforming
% plog_r_blm = repmat(squeeze(mean(mean(plog(:,:,:,blt),4),1)),[1,1,38]);
% plog_r_blstd = repmat(squeeze(mean(std(plog(:,:,:,blt),0,1),4)),[1,1,38]);
% plog_r_z = (squeeze(mean(plog,1)) - plog_r_blm) ./ plog_r_blstd;
% figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(plog_r_z(chan,:,:)),[-1 1]);axis xy;colorbar;
% title('robert z(log10(output)~blm)');

% dB(raw)
p(p == 0) = eps(0);
pdb = 10*log10(p);
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(pdb(:,chan,:,:),1)));axis xy;colorbar;
title('dB = 10*log10(raw power)');

% dB(raw ./ bl)
p(p == 0) = eps(0);
p_blm = repmat(mean(p(:,:,:,blt),4),[1,1,1,38]);
pdb_ers = 10*log10(p ./ p_blm);
figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(pdb_ers(:,chan,:,:),1)),[-4 4]);axis xy;colorbar;
title('10*log10(raw power ./ mean(bl))');

% % db(raw) ./ db(bl) - NOT WORKING
% p(p == 0) = eps(0);
% p_blm = repmat(mean(p(:,:,:,blt),4),[1,1,1,38]);
% pdb_ers_sep = (10*log10(p)) ./ (10*log10(p_blm));
% figure;imagesc(data_pow.(cond{cnd}).sub(sub).ses(ses).data.time,data_pow.(cond{cnd}).sub(sub).ses(ses).data.freq,squeeze(mean(pdb_ers_sep(:,chan,:,:),1)),[-4 4]);axis xy;colorbar;
% title('10*log10(raw power) ./ 10*log10(mean(bl))');




%% Test plots to make sure data look ok

cfg_ft = [];
% cfg_ft.baseline = [-0.3 -0.1];
% cfg_ft.baselinetype = 'absolute'; % maybe this
% %cfg_ft.baselinetype = 'relative';
% %cfg_ft.baselinetype = 'relchange'; % or this
% if strcmp(cfg_ft.baselinetype,'absolute')
%   cfg_ft.zlim = [-400 400];
%   %cfg_ft.zlim = [-2 2];
% elseif strcmp(cfg_ft.baselinetype,'relative')
%   cfg_ft.zlim = [0 2.0];
% elseif strcmp(cfg_ft.baselinetype,'relchange')
%   cfg_ft.zlim = [-1.0 1.0];
% end
cfg_ft.zlim = [-1.0 1.0];

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
  if isfield(cfg_ft,'baselinetype') && ~isempty(cfg_ft.baselinetype)
    title(sprintf('%s, baseline: %.1f, %.1f (%s)',ana.eventValues{1}{i},cfg_ft.baseline(1),cfg_ft.baseline(2),cfg_ft.baselinetype));
  else
    title(sprintf('%s',ana.eventValues{1}{i}));
  end
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
exper.badBehSub = {'SOSI011','SOSI030','SOSI007'}; % for ERP publication; also 001 has low trial counts

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
%cfg_ft.keepindividual = 'yes';
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

%% evoked plots

evVal=1;
chan=20;
clim=[-1.2 1.2];
figure;
imagesc(ga_evoked.(ana.eventValues{typ}{evVal}).time,ga_evoked.(ana.eventValues{typ}{evVal}).freq,squeeze(ga_evoked.(ana.eventValues{typ}{evVal}).powspctrm(chan,:,:)),clim);
axis xy;
colorbar;
title(sprintf('GA, %s, chan %d',ana.eventValues{typ}{evVal},chan));

%% plot the conditions - simple

cfg_ft = [];
%cfg_ft.baseline = [-0.3 -0.1];
%cfg_ft.baselinetype = 'absolute';
%if strcmp(cfg_ft.baselinetype,'absolute')
%cfg_ft.xlim = [0 1];
cfg_ft.ylim = [4 100];
%cfg_ft.ylim = [4 8];
%cfg_ft.ylim = [8 12];
%cfg_ft.zlim = [-0.4 0.4];
%cfg_ft.zlim = [-1.2 1.2];
%cfg_ft.zlim = [-1 1];
%cfg_ft.zlim = [-2 2];
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
    ft_multiplotTFR(cfg_ft,ga_pow.(ana.eventValues{typ}{evVal}));
    %ft_multiplotTFR(cfg_ft,ga_evoked.(ana.eventValues{typ}{evVal}));
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
cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'RCR','RH','RHSC','RHSI'}},size(cfg_plot.rois));

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.parameter = 'powspctrm';

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
cfg_plot.conditions = {'all'};

cfg_ft = [];
%cfg_ft.xlim = [.5 .8]; % time
%cfg_ft.ylim = [3 8]; % freq
cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.parameter = 'powspctrm';

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

mm_ft_contrastTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_pow);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
%cfg_ana.rois = {{'PS'},{'FS'}};
cfg_ana.rois = {{'LAS'},{'LPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.21 0.33; 0.6 1.0];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [35 80; 4 8];

%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};
%cfg_ana.conditions = {{'SC','SI'},{'SC','SI'}};
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
end

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LAS','RAS'},{'LPS','RPS'}};
% IV2: define the conditions tested for each set of ROIs
cfg_ana.condByROI = {{'CR','SC','SI'},{'CR','SC','SI'}};
%cfg_ana.condByROI = repmat({{'TH','NT'}},size(cfg_ana.rois));
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [3 8; 8 12];

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'powspctrm';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.frequency = cfg_ana.frequencies(r,:);
  
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
% cfg_ft.clusteralpha = 0.05;
cfg_ft.clusteralpha = 0.1;
cfg_ft.alpha = 0.05;

cfg_ana = [];
cfg_ana.roi = 'all';
cfg_ana.avgFrq = cfg_ft.avgoverfreq;
cfg_ana.conditions = {'all'};
%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};

% extra identifier when saving
%thisBLtype = ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baselinetype');
%thisBL = ft_findcfg(data_pow.(ana.eventValues{1}{1}).sub(1).ses(1).data.cfg,'baseline');

%thisBLtype = 'zpow';
%thisBLtype = 'zpow_indu';
thisBLtype = 'pow_induced';
%thisBLtype = 'zp_evok';
%thisBLtype = 'coh';
thisBL = [-0.4 -0.2];
cfg_ana.dirStr = sprintf('_%s_%d_%d',thisBLtype,thisBL(1)*1000,thisBL(2)*1000);

if strcmp(cfg_ft.avgovertime,'no')
  cfg_ana.latencies = [0 1.0];
  %cfg_ana.latencies = [0 0.5; 0.5 1.0];
elseif strcmp(cfg_ft.avgovertime,'yes')
  %cfg_ana.latencies = [-0.2:0.1:0.9; -0.1:0.1:1.0]';
  cfg_ana.latencies = [-0.2:0.2:0.8; 0:0.2:1.0]';
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end

if strcmp(cfg_ft.avgoverfreq,'no')
  cfg_ana.frequencies = [4 100];
elseif strcmp(cfg_ft.avgoverfreq,'yes')
  %cfg_ana.frequencies = [4 8; 8 12; 12 28; 28 40];
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

files.saveFigs = 1;

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
cfg_ft.elec = ga_pow.elec;
cfg_ft.highlight = 'on';
if cfg_ana.include_clus_stat == 0
  cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS','RPS','LPS'})});
end
%cfg_ft.comment = 'xlim';
cfg_ft.commentpos = 'title';

% create contrast
cont_topo = [];
cont_topo.RHvsRCR = ga_pow.RH;
cont_topo.RHvsRCR.avg = ga_pow.RH.avg - ga_pow.RCR.avg;
cont_topo.RHvsRCR.individual = ga_pow.RH.individual - ga_pow.RCR.individual;
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
cont_topo.RHSCvsRHSI = ga_pow.RHSC;
cont_topo.RHSCvsRHSI.avg = ga_pow.RHSC.avg - ga_pow.RHSI.avg;
cont_topo.RHSCvsRHSI.individual = ga_pow.RHSC.individual - ga_pow.RHSI.individual;
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
cont_topo.RHSCvsRCR = ga_pow.RHSC;
cont_topo.RHSCvsRCR.avg = ga_pow.RHSC.avg - ga_pow.RCR.avg;
cont_topo.RHSCvsRCR.individual = ga_pow.RHSC.individual - ga_pow.RCR.individual;
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
cont_topo.RHSIvsRCR = ga_pow.RHSI;
cont_topo.RHSIvsRCR.avg = ga_pow.RHSI.avg - ga_pow.RCR.avg;
cont_topo.RHSIvsRCR.individual = ga_pow.RHSI.individual - ga_pow.RCR.individual;
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
% cfg_ft.layout = ft_prepare_layout(cfg_ft,ga_pow);
% cfg_ft.contournum = 0;
% cfg_ft.emarker = '.';
% cfg_ft.alpha  = 0.05;
% cfg_ft.parameter = 'stat';
% cfg_ft.zlim = [-5 5];
% ft_clusterplot(cfg_ft,stat_clus.RHSCvsRHSIvsRCR);
