% Do a nonparametric statistics clustering analysis for timelocked EEG
% (ERPs)

% See Maris and Oostenveld (2007) for a description of statistics

%% Experiment-specific setup

expName = 'SLORK';

% pre- and post-stimulus times to save, in seconds
prepost = [-1.0 2.0];

sampleRate = 250;

%eventValues = sort({'CR','HSC','HSI','TCR','THSC','THSI'});

%% list subjects
subjects = {
    'SLORK002'...
    'SLORK003'...
    'SLORK004'...
    'SLORK005'...
    'SLORK006'...
    'SLORK007'...
    'SLORK008'...
    'SLORK009'...
    'SLORK010'...
    'SLORK011'...
    'SLORK012'...
    'SLORK013'...
    'SLORK014'...
    'SLORK015'...
    'SLORK016'...
    'SLORK017'...
    'SLORK018'...
    'SLORK019'...
    'SLORK020'...
    'SLORK022'...
    'SLORK023'...
    'SLORK024'...
    'SLORK025'...
    'SLORK026'...
    'SLORK027'...
    'SLORK029'...
    'SLORK030'...
    'SLORK031'...
    'SLORK032'...
    'SLORK033'...
    };
  % 'SLORK001'... % exp was not set up correctly
  % 'SLORK021'... % didn't record
  % 'SLORK028'... % did't record

%% set up parameters

homeDir = getenv('HOME');

% EP post processed
serverDir = fullfile('/Volumes/curranlab/Data',expName,'eeg/eppp_nobackup/-1000_2000');
serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data',expName,'eeg/eppp_nobackup/-1000_2000');

% NS post processed
%serverDir = fullfile('/Volumes/curranlab/Data/',expName,'/eeg/nspp_nobackup/-1.0_2.0');
%serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data/',expName,'/eeg/nspp_nobackup/-1.0_2.0');

% pick the right dataroot
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  dataroot = fullfile(homeDir,'data',expName,'eeg');
end

% directory to save the data
saveDir = fullfile(dataroot,'ft_data');
if ~exist(saveDir,'dir')
  mkdir(saveDir);
end

% directory to save figures
saveDirFigs = fullfile(saveDir,'figs');
if ~exist(saveDirFigs,'dir')
  mkdir(saveDirFigs)
end

% Assumes we have chan locs file in ~/eeg/
elecfile = fullfile(homeDir,'eeg/GSN_HydroCel_129_short.sfp');
locsFormat = 'besa_sfp';
elec = read_sens(elecfile,'fileformat',locsFormat);

% figure printing options - see mm_ft_setSaveDirs for other options
saveFigs = 1;
figFileExt = 'eps';
if strcmp(figFileExt,'eps')
  figPrintFormat = '-depsc2';
elseif strcmp(figFileExt,'pdf')
  figPrintFormat = '-dpdf';
elseif strcmp(figFileExt,'png')
  figPrintFormat = '-dpng';
end

%% Save the data in FieldTrip structs

for sub = 1:length(subjects)
  fprintf('Working on %s...\n',subjects{sub});
  
  filename = sprintf('%s_data_preproc_%d_%d.mat',subjects{sub},prepost(1)*1000,prepost(2)*1000);
  dataSaveFile = fullfile(saveDir,filename);
  
  if ~exist(dataSaveFile,'file')
    fprintf('Creating %s FieldTrip struct...\n',subjects{sub});
    data_preproc = seg2ft_nsmat(dataroot,subjects{sub},prepost,sampleRate,elecfile);
    
    fprintf('Saving %s FieldTrip struct...\n',subjects{sub});
    save(dataSaveFile,'data_preproc','prepost','elec');
  end
  
  % add NS's artifact information to the event structure
  ns_addArtifactInfo_slork(dataroot,subjects{sub});

  fprintf('Done with %s.\n',subjects{sub});
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Testing
% 
% cfg = [];
% cfg.channel = {'all'};
% %cfg.channel = {'all','-EOGH','-EOGV_L','-EOGV_R'};
% %cfg.latency = [.5 .8];
% data_avg = ft_timelockanalysis(cfg,data_preproc(1));
% 
% cfg = [];
% cfg.showlabels = 'yes'; 
% cfg.fontsize = 9; 
% cfg.layout = data_avg.elec;
% %cfg.ylim = [-3e-13 3e-13];
% ft_multiplotER(cfg,data_avg); 
% 
% cfg = [];
% cfg.channel = {'e42'};
% ft_singleplotER(cfg,data_avg);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cfg = [];
% cfg.channel = {'all'};
% %cfg.channel = {'all','-EOGH','-EOGV_L','-EOGV_R'};
% cfg.latency = [.5 .8];
% data_topo = timelockanalysis(cfg,data_preproc);
% cfg = [];
% cfg.fontsize = 9; 
% ft_topoplotER(cfg,data_topo);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cfg = [];
% cfg.channel = {'all'};
% %cfg.channel = {'all','-EOGH','-EOGV_L','-EOGV_R'};
% %cfg.latency = [.5 .8];
% data_avg = ft_timelockanalysis(cfg,data_HSC(1));
% 
% cfg = [];
% cfg.showlabels = 'yes'; 
% cfg.fontsize = 9; 
% cfg.layout = data_avg.elec;
% %cfg.ylim = [-3e-13 3e-13];
% figure
% ft_multiplotER(cfg,data_avg); 

%% Run ft_timelockanalysis; make sure events go in alphabetical order

%clear data_CRM data_CRS data_SFAF data_SFAO data_SFAS data_SMM data_SMS data_SSCF data_SSCO data_SSCS data_SSIF data_SSIO data_SSIS data_TCRM data_TCRS data_TFAF data_TFAO data_TFAS data_TMM data_TMS data_TSCF data_TSCO data_TSCS data_TSIF data_TSIO data_TSIS

%clear data_CR data_HSC data_HSI data_TCR data_THSC data_THSI

clear data_* data

cfg_pp = [];
% % baseline correction is already done by EP Toolkit
% cfg_pp.blc = 'yes';
% cfg_pp.blcwindow  = [-0.2 0];

% apply lowpass filter at 35 Hz
%cfg_pp.lpfilter = 'yes';
%cfg_pp.lpfreq = 35;

cfg_proc = [];
cfg_proc.keeptrials = 'yes';

% cfg_tlb = [];
% cfg_tlb.baseline = [-0.2 0];
% cfg_tlb.channel = 'all';

% for keeping track of the number of trials that each subject has
num_CR = zeros(length(subjects),1);
num_HSC = zeros(length(subjects),1);
num_H = zeros(length(subjects),1);
num_HSI = zeros(length(subjects),1);
num_TCR = zeros(length(subjects),1);
num_TH = zeros(length(subjects),1);
num_THSC = zeros(length(subjects),1);
num_THSI = zeros(length(subjects),1);

%matlabpool open

%parfor sub = 1:length(subjects)
for sub = 1:length(subjects)
  filename = sprintf('%s_data_preproc_%d_%d.mat',subjects{sub},prepost(1)*1000,prepost(2)*1000);
  dataSaveFile = fullfile(saveDir,filename);
  fprintf('Loading %s...\n',subjects{sub});
  data = load(dataSaveFile);
  
  for i = 1:length(data.data_preproc)
    eventVal = data.data_preproc(i).cfg.trialdef.eventvalue;
    
    if data.data_preproc(i).hdr.nTrials > 0
      
      % do any preprocessing
      if ~isempty(cfg_pp)
        pp_str = sprintf('data_%s_pp = ft_preprocessing(cfg_pp,data.data_preproc(i));',eventVal);
      else
        pp_str = sprintf('data_%s_pp = data.data_preproc(i);',eventVal);
      end
      eval(pp_str);
      
      % run ft_timelockanalysis
      tla_str = sprintf('data_%s(sub) = ft_timelockanalysis(cfg_proc,data_%s_pp);',eventVal,eventVal);
      eval(tla_str);
      
      %%%%%%%%
      
      %tla_str = sprintf('data_%s_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(i));',eventVal);
      %eval(tla_str);
      
      %tlb_str = sprintf('data_%s(sub) = ft_timelockbaseline(cfg_tlb,data_%s_raw);',eventVal,eventVal);
      %eval(tlb_str);
      
      % get the var back
      %tla_str = sprintf('data_%s_var(sub) = ft_timelockanalysis(cfg_proc,data_%s(sub));',eventVal,eventVal);
      %eval(tla_str);
      
    else
      eval(sprintf('data_%s(sub).trial = [];',eventVal));
      eval(sprintf('data_%s(sub).label = {};',eventVal));
      eval(sprintf('data_%s(sub).time = [];',eventVal));
      eval(sprintf('data_%s(sub).hdr = [];',eventVal));
      eval(sprintf('data_%s(sub).fsample = [];',eventVal));
      eval(sprintf('data_%s(sub).elec = [];',eventVal));
      eval(sprintf('data_%s(sub).cfg = [];',eventVal));
      eval(sprintf('data_%s(sub).dimord = [];',eventVal));
      
%       eval(sprintf('data_%s_pp(sub).trial = [];',eventVal));
%       eval(sprintf('data_%s_pp(sub).label = {};',eventVal));
%       eval(sprintf('data_%s_pp(sub).time = [];',eventVal));
%       eval(sprintf('data_%s_pp(sub).hdr = [];',eventVal));
%       eval(sprintf('data_%s_pp(sub).fsample = [];',eventVal));
%       eval(sprintf('data_%s_pp(sub).elec = [];',eventVal));
%       eval(sprintf('data_%s_pp(sub).cfg = [];',eventVal));
%       eval(sprintf('data_%s_pp(sub).dimord = [];',eventVal));
    end
    
    % keep track of the number of trials for each subject
    num_str = sprintf('num_%s(sub) = size(data_%s(sub).trial,1);',eventVal,eventVal);
    eval(num_str);
    
  end
  
  if ~isempty(data_HSC(sub).trial) && ~isempty(data_HSI(sub).trial)
    data_H_pp = ft_appenddata([],data_HSC_pp,data_HSI_pp);
  elseif ~isempty(data_HSC(sub).trial) && isempty(data_HSI(sub).trial)
    data_H_pp = data_HSC_pp;
  elseif isempty(data_HSC(sub).trial) && ~isempty(data_HSI(sub).trial)
    data_H_pp = data_HSI_pp;
  else
    data_H_pp.trial = [];
  end
  if ~isempty(data_H_pp.trial)
    data_H(sub) = ft_timelockanalysis(cfg_proc,data_H_pp);
  else
    data_H(sub).trial = [];
  end
  
  if ~isempty(data_THSC(sub).trial) && ~isempty(data_THSI(sub).trial)
    data_TH_pp = ft_appenddata([],data_THSC_pp,data_THSI_pp);
  elseif ~isempty(data_THSC(sub).trial) && isempty(data_THSI(sub).trial)
    data_TH_pp = data_THSC_pp;
  elseif isempty(data_THSC(sub).trial) && ~isempty(data_THSI(sub).trial)
    data_TH_pp = data_THSI_pp;
  else
    data_TH_pp.trial = [];
  end
  if ~isempty(data_TH_pp.trial)
    data_TH(sub) = ft_timelockanalysis(cfg_proc,data_TH_pp);
  else
    data_TH(sub).trial = [];
  end
  
  num_H(sub) = num_HSC(sub) + num_HSI(sub);
  num_TH(sub) = num_THSC(sub) + num_THSI(sub);
  
  %eval(sprintf('clear data_*_raw'));
  clear data_*_pp
  
%   % Just the 3 event categories
%   if data.data_preproc(1).hdr.nTrials > 0
%     data_CR_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(1));
%     data_CR(sub) = ft_timelockbaseline(cfg_tlb,data_CR_raw);
%     % get the var field back
%     %data_CR_var(sub) = ft_timelockanalysis(cfg_proc,data_CR(sub));
%   else
%     data_CR(sub).trial = [];
%   end
%   if data.data_preproc(2).hdr.nTrials > 0
%     data_HSC_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(2));
%     data_HSC(sub) = ft_timelockbaseline(cfg_tlb,data_HSC_raw);
%     % get the var field back
%     %data_HSC_var(sub) = ft_timelockanalysis(cfg_proc,data_HSC(sub));
%   else
%     data_HSC(sub).trial = [];
%   end
%   if data.data_preproc(3).hdr.nTrials > 0
%     data_HSI_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(3));
%     data_HSI(sub) = ft_timelockbaseline(cfg_tlb,data_HSI_raw);
%     % get the var field back
%     %data_HSI_var(sub) = ft_timelockanalysis(cfg_proc,data_HSI(sub));
%   else
%     data_HSI(sub).trial = [];
%   end
%   
%   if data.data_preproc(4).hdr.nTrials > 0
%     data_TCR_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(4));
%     data_TCR(sub) = ft_timelockbaseline(cfg_tlb,data_TCR_raw);
%     % get the var field back
%     %data_TCR_var(sub) = ft_timelockanalysis(cfg_proc,data_TCR(sub));
%   else
%     data_TCR(sub).trial = [];
%   end
%   if data.data_preproc(5).hdr.nTrials > 0
%     data_THSC_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(5));
%     data_THSC(sub) = ft_timelockbaseline(cfg_tlb,data_THSC_raw);
%     % get the var field back
%     %data_THSC_var(sub) = ft_timelockanalysis(cfg_proc,data_THSC(sub));
%   else
%     data_THSC(sub).trial = [];
%   end
%   if data.data_preproc(6).hdr.nTrials > 0
%     data_THSI_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(6));
%     data_THSI(sub) = ft_timelockbaseline(cfg_tlb,data_THSI_raw);
%     % get the var field back
%     %data_THSI_var(sub) = ft_timelockanalysis(cfg_proc,data_THSI(sub));
%   else
%     data_THSI(sub).trial = [];
%   end
  
%   data_CRM(sub) = ft_timelockanalysis(cfg,data.data_preproc(1));
%   data_CRS(sub) = ft_timelockanalysis(cfg,data.data_preproc(2));
%   data_SFAF(sub) = ft_timelockanalysis(cfg,data.data_preproc(3));
%   data_SFAO(sub) = ft_timelockanalysis(cfg,data.data_preproc(4));
%   data_SFAS(sub) = ft_timelockanalysis(cfg,data.data_preproc(5));
%   data_SMM(sub) = ft_timelockanalysis(cfg,data.data_preproc(6));
%   data_SMS(sub) = ft_timelockanalysis(cfg,data.data_preproc(7));
%   data_SSCF(sub) = ft_timelockanalysis(cfg,data.data_preproc(8));
%   data_SSCO(sub) = ft_timelockanalysis(cfg,data.data_preproc(9));
%   data_SSCS(sub) = ft_timelockanalysis(cfg,data.data_preproc(10));
%   data_SSIF(sub) = ft_timelockanalysis(cfg,data.data_preproc(11));
%   data_SSIO(sub) = ft_timelockanalysis(cfg,data.data_preproc(12));
%   data_SSIS(sub) = ft_timelockanalysis(cfg,data.data_preproc(13));
%   data_TCRM(sub) = ft_timelockanalysis(cfg,data.data_preproc(14));
%   data_TCRS(sub) = ft_timelockanalysis(cfg,data.data_preproc(15));
%   data_TFAF(sub) = ft_timelockanalysis(cfg,data.data_preproc(16));
%   data_TFAO(sub) = ft_timelockanalysis(cfg,data.data_preproc(17));
%   data_TFAS(sub) = ft_timelockanalysis(cfg,data.data_preproc(18));
%   data_TMM(sub) = ft_timelockanalysis(cfg,data.data_preproc(19));
%   data_TMS(sub) = ft_timelockanalysis(cfg,data.data_preproc(20));
%   data_TSCF(sub) = ft_timelockanalysis(cfg,data.data_preproc(21));
%   data_TSCO(sub) = ft_timelockanalysis(cfg,data.data_preproc(22));
%   data_TSCS(sub) = ft_timelockanalysis(cfg,data.data_preproc(23));
%   data_TSIF(sub) = ft_timelockanalysis(cfg,data.data_preproc(24));
%   data_TSIO(sub) = ft_timelockanalysis(cfg,data.data_preproc(25));
%   data_TSIS(sub) = ft_timelockanalysis(cfg,data.data_preproc(26));
  
  fprintf('Done.\n');
end

%matlabpool close

% % ft_appenddata must be done after preprocessing, before ft_timelockanalysis
% cfg = [];
% for sub = 1:length(subjects)
%   data_HSC_pp(sub) = ft_appenddata(cfg,data_SSCS_pp(sub),data_SSCO_pp(sub),data_SSCF_pp(sub));
%   data_HSI_pp(sub) = ft_appenddata(cfg,data_SSIS_pp(sub),data_SSIO_pp(sub),data_SSIF_pp(sub));
%   data_CR_pp(sub) = ft_appenddata(cfg,data_CRS_pp(sub),data_CRM_pp(sub));
%   data_THSC_pp(sub) = ft_appenddata(cfg,data_TSCS_pp(sub),data_TSCO_pp(sub),data_TSCF_pp(sub));
%   data_THSI_pp(sub) = ft_appenddata(cfg,data_TSIS_pp(sub),data_TSIO_pp(sub),data_TSIF_pp(sub));
%   data_TCR_pp(sub) = ft_appenddata(cfg,data_TCRS_pp(sub),data_TCRM_pp(sub));
% end

%% save the ft_timelockanalysis files for each event

tlaFilename = sprintf('data_CR_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_CR','num_CR');

tlaFilename = sprintf('data_H_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_H','num_H');

tlaFilename = sprintf('data_HSC_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_HSC','num_HSC');

tlaFilename = sprintf('data_HSI_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_HSI','num_HSI');

%% if already saved, load the ft_timelockanalysis files

savedFiles = dir(fullfile(saveDir,'data_*_tla_*.mat'));
for i = 1:length(savedFiles)
  fprintf('Loading %s...',savedFiles(i).name);
  load(fullfile(saveDir,savedFiles(i).name));
  fprintf('Done.\n');
end

cfg = [];
for sub = 1:length(subjects)
  data_CR_temp(sub) = ft_appenddata(cfg,data_SCR(sub),data_TCR(sub));
  data_H_temp(sub) = ft_appenddata(cfg,data_SH(sub),data_TH(sub));
  data_HSC_temp(sub) = ft_appenddata(cfg,data_SHSC(sub),data_THSC(sub));
  data_HSI_temp(sub) = ft_appenddata(cfg,data_SHSI(sub),data_THSI(sub));
  
  data_CR(sub) = ft_timelockanalysis(cfg,data_CR_temp(sub));
  data_H(sub) = ft_timelockanalysis(cfg,data_H_temp(sub));
  data_HSC(sub) = ft_timelockanalysis(cfg,data_HSC_temp(sub));
  data_HSI(sub) = ft_timelockanalysis(cfg,data_HSI_temp(sub));
  
  num_CR(sub) = num_SCR(sub) + num_TCR(sub);
  num_H(sub) = num_SH(sub) + num_TH(sub);
  num_HSC(sub) = num_SHSC(sub) + num_THSC(sub);
  num_HSI(sub) = num_SHSI(sub) + num_THSI(sub);
end
num_CR = num_CR';
num_H = num_H';
num_HSC = num_HSC';
num_HSI = num_HSI';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove the degrees of freedom? field. I'm not sure what it does.
data_CR = rmfield(data_CR,'dof');
data_H = rmfield(data_H,'dof');
data_HSC = rmfield(data_HSC,'dof');
data_HSI = rmfield(data_HSI,'dof');

%% decide who to kick out based on trial counts
[num_CR num_HSC num_HSI]

lowNum = [];
thresh = 15;
% get the number of trials
for sub = 1:length(subjects)
  if num_CR(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_HSC(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_HSI(sub) < thresh
    lowNum = [lowNum sub];
  end
end

% who we're rejecting
subjects(unique(lowNum))

% should be: 'SLORK009'    'SLORK015'    'SLORK016'    'SLORK024'
% 'SLORK026'    'SLORK032'    'SLORK033'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

badBehSub = {};

% item d' < 1.0
%badBehSub = {'SLORK010','SLORK011','SLORK030'};
%badBehSub = {'SLORK011','SLORK017'};

% SLORK012 and SLORK030's EEG is really noisy in LAS and RAS
%badBehSub = {'SLORK012','SLORK030'};

badSub = unique([badBehSub subjects(lowNum)]);

fprintf('Number of events included in EEG analyses (%d subjects; threshold: %d events):\n',sum(~ismember(subjects,badSub)),thresh);
cr_mean = mean(num_CR(~ismember(subjects,badSub)));
cr_sem = std(num_CR(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('CR: M = %.3f, SEM = %.3f\n',cr_mean,cr_sem);
hsc_mean = mean(num_HSC(~ismember(subjects,badSub)));
hsc_sem = std(num_HSC(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('HSC: M = %.3f, SEM = %.3f\n',hsc_mean,hsc_sem);
hsi_mean = mean(num_HSI(~ismember(subjects,badSub)));
hsi_sem = std(num_HSI(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('HSI: M = %.3f, SEM = %.3f\n',hsi_mean,hsi_sem);
h_mean = mean(num_H(~ismember(subjects,badSub)));
h_sem = std(num_H(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('H: M = %.3f, SEM = %.3f\n',h_mean,h_sem);

%% set up analysis strings

str_CR = [];
str_H = [];
str_HSC = [];
str_HSI = [];

for sub = 1:length(subjects)
  if ismember(subjects{sub},badSub)
    continue
  else
    if sub == 1
      str_CR = sprintf('data_CR(%d)',sub);
      str_HSC = sprintf('data_HSC(%d)',sub);
      str_HSI = sprintf('data_HSI(%d)',sub);
      str_H = sprintf('data_H(%d)',sub);
    elseif sub > 1
      str_CR = sprintf('%s,data_CR(%d)',str_CR,sub);
      str_HSC = sprintf('%s,data_HSC(%d)',str_HSC,sub);
      str_HSI = sprintf('%s,data_HSI(%d)',str_HSI,sub);
      str_H = sprintf('%s,data_H(%d)',str_H,sub);
    end
  end
end
if strcmp(str_HSC(1),',')
  str_CR = str_CR(2:end);
  str_HSC = str_HSC(2:end);
  str_HSI = str_HSI(2:end);
  str_H = str_H(2:end);
end

%% get the grand average and add the electrode locations again

cfg = [];
cfg.keepindividual = 'yes';
eval(sprintf('ga_CR = ft_timelockgrandaverage(cfg,%s);',str_CR));
eval(sprintf('ga_HSC = ft_timelockgrandaverage(cfg,%s);',str_HSC));
eval(sprintf('ga_HSI = ft_timelockgrandaverage(cfg,%s);',str_HSI));
eval(sprintf('ga_H = ft_timelockgrandaverage(cfg,%s);',str_H));

ga_CR.elec = elec;
ga_HSC.elec = elec;
ga_HSI.elec = elec;
ga_H.elec = elec;

%% set up channel groups

elecGroups = {...
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
elecGroupsStr = {'LAI','RAI','LAS','RAS','LPS','RPS','LPI','RPI','FC','LPS2','RPS2','P1','N1'};

% %% plot the conditions
% cfg = [];
% cfg.showlabels = 'yes'; 
% cfg.fontsize = 9; 
% %cfg.ylim = [-3e-13 3e-13];
% 
% cfg.layout = ft_prepare_layout(cfg,ga_CR);
% figure
% ft_multiplotER(cfg,ga_CR,ga_HSC,ga_HSI);
% figure
% ft_multiplotER(cfg,ga_TCR,ga_THSC,ga_THSI);

%cfg.interplimits = 'electrodes';
%ft_topoplotER(cfg,ga_CR);

%% plot the conditions - simple

cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 9;
cfg.xlim = [-.2 2.0];
cfg.linewidth = 2;
cfg.graphcolor = 'krb';
cfg.linestyle = {'-','--','-.'};
% FN400
cfg.ylim = [-4 3];
cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'LAS','RAS'})});
figure
ft_singleplotER(cfg,ga_CR,ga_HSC,ga_HSI);
legend('CR','HSC','HSI','Location','SouthEast');
% P O/N
cfg.ylim = [-2 5];
cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'LPS','RPS'})});
figure
ft_singleplotER(cfg,ga_CR,ga_HSC,ga_HSI);
legend('CR','HSC','HSI','Location','NorthEast');

%% subplot with all subjects' ERPs

roi = {'LAS','RAS'};
%roi = {'LPS','RPS'};
%roi = {'RAS'};
%roi = {'LPS'};
chan_str = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
chan = str2double(strrep(chan_str,'e',''));
numCols = 5;
excludeBadSub = 0;
if excludeBadSub
  numRows = ceil((length(subjects) - length(badSub)) / numCols);
else
  numRows = ceil((length(subjects)) / numCols) + 1;
end

% side
testType = 'both';
figure
for sub = 1:length(subjects)
  subplot(numRows,numCols,sub);
  plot(data_CR(sub).time,mean(data_CR(sub).avg(chan,:),1),'b');
  hold on
  plot(data_HSC(sub).time,mean(data_HSC(sub).avg(chan,:),1),'r');
  plot(data_HSI(sub).time,mean(data_HSI(sub).avg(chan,:),1),'g');
  plot(data_H(sub).time,mean(data_H(sub).avg(chan,:),1),'k');
  if ismember(subjects{sub},badSub)
    title(sprintf('*BAD*: %s; CR:%d, SC:%d, SI:%d',subjects{sub},size(data_CR(sub).trial,1),size(data_HSC(sub).trial,1),size(data_HSI(sub).trial,1)));
  else
    title(sprintf('%s; CR:%d, SC:%d, SI:%d',subjects{sub},size(data_CR(sub).trial,1),size(data_HSC(sub).trial,1),size(data_HSI(sub).trial,1)));
  end
  %ylim([0 1.9e-13]);
  axis([-0.2 1.5 -10 10]);
  hold off
end
set(gcf,'Name',[testType,': ',sprintf(repmat('%s ',1,length(roi)),roi{:})])
subplot(numRows,numCols,length(subjects) + 1);
text(0.5,0.9,[testType,', ',sprintf(repmat('%s ',1,length(roi)),roi{:})],'color','k')
text(0.5,0.7,'CR','color','b')
text(0.5,0.5,'HSC','color','r') ;
text(0.5,0.3,'HSI','color','g')
text(0.5,0.1,'H','color','k')
axis off

%% plot the conditions

cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 9;
cfg.layout = ft_prepare_layout(cfg,ga_CR);
%fillcolor = [.8,.8,.8];
%fillalpha = 0.3;
%filledge = 'none';
minMaxVolt = [-4 3; -2 5];
latency = [-0.2 2.0];
plotLegend = 1;

cfg.linewidth = 2;
cfg.graphcolor = 'krb';
cfg.linestyle = {'-','--','-.'};

plot_rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'},{'LAS','RAS'},{'LPS','RPS'}};
%plot_rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};

% both
testType = 'both';
for i = 1:length(plot_rois)
  roi = plot_rois{i};
  cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
  figure
  %hold on
  %ft_singleplotER(cfg,ga_CR,ga_HSC,ga_HSI,ga_H);
  ft_singleplotER(cfg,ga_CR,ga_HSC,ga_HSI);
  hold on
  plot([latency(1) latency(2)],[0 0],'k--'); % horizontal
  if ismember('LAS',roi) || ismember('RAS',roi)
    plot([0 0],[minMaxVolt(1,1) minMaxVolt(1,2)],'k--'); % vertical
    if plotLegend
      %legend('CR','HSC','HSI','H','Location','SouthEast');
      legend('CR','HSC','HSI','Location','SouthEast');
    end
    
    vertlatency = [0.3 0.5];
    %h = fill([vertlatency(1),vertlatency(2),vertlatency(2),vertlatency(1)],[minMaxVolt(1,1),minMaxVolt(1,1),minMaxVolt(1,2),minMaxVolt(1,2)],fillcolor);
    %set(h,'FaceAlpha',fillalpha);
    %set(h,'EdgeColor',filledge)
    plot([vertlatency(1) vertlatency(1)],[minMaxVolt(1,1) minMaxVolt(1,2)],'k'); % vertical
    plot([vertlatency(2) vertlatency(2)],[minMaxVolt(1,1) minMaxVolt(1,2)],'k'); % vertical
    
    axis([latency(1) latency(2) minMaxVolt(1,1) minMaxVolt(1,2)])
  elseif ismember('LPS',roi) || ismember('RPS',roi)
    plot([0 0],[minMaxVolt(2,1) minMaxVolt(2,2)],'k--'); % vertical
    if plotLegend
      %legend('CR','HSC','HSI','H','Location','NorthEast');
      legend('CR','HSC','HSI','Location','NorthEast');
    end
    
    vertlatency = [0.5 0.8];
    %h = fill([vertlatency(1),vertlatency(2),vertlatency(2),vertlatency(1)],[minMaxVolt(2,1),minMaxVolt(2,1),minMaxVolt(2,2),minMaxVolt(2,2)],fillcolor);
    %set(h,'FaceAlpha',fillalpha);
    %set(h,'EdgeColor',filledge);
    plot([vertlatency(1) vertlatency(1)],[minMaxVolt(2,1) minMaxVolt(2,2)],'k'); % vertical
    plot([vertlatency(2) vertlatency(2)],[minMaxVolt(2,1) minMaxVolt(2,2)],'k'); % vertical
    
    axis([latency(1) latency(2) minMaxVolt(2,1) minMaxVolt(2,2)])
  end
  xlabel('Time (s)');
  ylabel('Voltage (\muV)');
  set(gcf,'Name',[sprintf(repmat('%s ',1,length(roi)),roi{:})])
  %title(sprintf(repmat('%s ',1,length(roi)),roi{:}))
  %axis ij % negative up
  publishfig(gcf,1);
  if plotLegend
    legend_str = '_legend';
  else
    legend_str = '';
  end
  figfilename = sprintf('erp_ga_%s_%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(roi)),roi{:}),latency(1)*1000,latency(2)*1000,legend_str,figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
end

%% Plot average time contrasts - both

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Side: LAS + RAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.interactive = 'no';
minMaxVolt = [-1 1];
%cfg.colorbar = 'yes';
cfg.layout = ft_prepare_layout(cfg,ga_HSC);
cfg.highlight = 'on';
cfg.highlightsize = 10;
%cfg.comment = 'xlim';
%cfg.commentpos = 'title';
cfg.comment = 'no';
plotTitle = 1;
plotColorbar = 1;
if plotColorbar
  colorbar_str = '_cb';
else
  colorbar_str = '';
end

cfg.colormap = colormap('jet'); % default; blue to red
%cfg.colormap = colormap('hot'); % dark to light; better for b&w printers

testType = 'both';

% create contrast
ga_HvsCR = ga_H;
ga_HvsCR.avg = ga_H.avg - ga_CR.avg;
ga_HvsCR.individual = ga_H.individual - ga_CR.individual;
% make a plot
figure
roi = {'LAS','RAS'};
cond = {'H','CR'};
cfg.xlim = [0.3 0.5];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Hits - Correct Rejections',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_HSCvsCR = ga_HSC;
ga_HSCvsCR.avg = ga_HSC.avg - ga_CR.avg;
ga_HSCvsCR.individual = ga_HSC.individual - ga_CR.individual;
% make a plot
figure
roi = {'LAS','RAS'};
cond = {'HSC','CR'};
cfg.xlim = [0.3 0.5];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HSCvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Source Correct - Correct Rejections',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_HSIvsCR = ga_HSI;
ga_HSIvsCR.avg = ga_HSI.avg - ga_CR.avg;
ga_HSIvsCR.individual = ga_HSI.individual - ga_CR.individual;
% make a plot
figure
roi = {'LAS','RAS'};
cond = {'HSI','CR'};
cfg.xlim = [0.3 0.5];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HSIvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Source Incorrect - Correct Rejections',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_HSCvsHSI = ga_HSC;
ga_HSCvsHSI.avg = ga_HSC.avg - ga_HSI.avg;
ga_HSCvsHSI.individual = ga_HSC.individual - ga_HSI.individual;
% make a plot
figure
roi = {'LAS','RAS'};
cond = {'HSC','HSI'};
cfg.xlim = [0.3 0.5];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HSCvsHSI);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Source Correct - Source Incorrect',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Side: LPS + RPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create contrast
ga_HvsCR = ga_H;
ga_HvsCR.avg = ga_H.avg - ga_CR.avg;
ga_HvsCR.individual = ga_H.individual - ga_CR.individual;
% make a plot
figure
roi = {'LPS','RPS'};
cond = {'H','CR'};
cfg.xlim = [0.5 0.8];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Hits - Correct Rejections',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_HSCvsCR = ga_HSC;
ga_HSCvsCR.avg = ga_HSC.avg - ga_CR.avg;
ga_HSCvsCR.individual = ga_HSC.individual - ga_CR.individual;
% make a plot
figure
roi = {'LPS','RPS'};
cond = {'HSC','CR'};
cfg.xlim = [0.5 0.8];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HSCvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Source Correct - Correct Rejections',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_HSIvsCR = ga_HSI;
ga_HSIvsCR.avg = ga_HSI.avg - ga_CR.avg;
ga_HSIvsCR.individual = ga_HSI.individual - ga_CR.individual;
% make a plot
figure
roi = {'LPS','RPS'};
cond = {'HSI','CR'};
cfg.xlim = [0.5 0.8];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HSIvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Source Incorrect - Correct Rejections',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_HSCvsHSI = ga_HSC;
ga_HSCvsHSI.avg = ga_HSC.avg - ga_HSI.avg;
ga_HSCvsHSI.individual = ga_HSC.individual - ga_HSI.individual;
% make a plot
figure
roi = {'LPS','RPS'};
cond = {'HSC','HSI'};
cfg.xlim = [0.5 0.8];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_HSCvsHSI);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('%s, Source Correct - Source Incorrect',testType);
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

%% correlations

% get the electrode numbers
elecGroupsNum = cell(size(elecGroups));
for i = 1:length(elecGroups)
  elecGroupsNum{i} = str2double(strrep(elecGroups{i},'e',''));
end

% accuracy
d_item_side = abs([1.61 2.08 2.70 1.80 1.98 1.10 2.09 2.09 0.66 0.64 1.68 1.79 1.60 2.45 2.30 1.73 1.57 1.25 1.08 3.03 1.96 1.50 0.88 2.78 1.55 1.35 0.85 0.84 0.71 1.21]);
d_source_side = abs([1.55 1.75 1.95 1.75 1.69 0.08 1.29 2.00 0.48 0.66 1.13 1.88 1.15 3.05 1.44 0.80 1.71 0.88 1.19 1.74 1.71 2.33 0.10 3.28 1.23 0.94 -0.61 0.70 0.20 0.62]);
d_item_q = abs([1.87 2.01 3.01 1.16 1.62 1.19 1.67 1.78 0.99 1.06 1.54 2.11 1.49 2.68 2.32 1.74 1.76 1.54 1.44 2.41 2.22 1.71 1.39 2.89 1.30 1.01 1.03 1.26 1.60 1.54]);
d_source_q = abs([0.74 1.06 1.54 0.74 0.80 -0.32 1.07 1.32 0.07 0.52 1.36 1.23 0.74 1.62 2.17 0.52 1.47 1.14 0.36 0.81 1.07 1.89 -0.25 1.36 1.21 0.40 -0.11 1.22 0.51 0.86]);

%anaComponent = 'FN400';
anaComponent = 'PON';

if strcmp(anaComponent,'FN400')
  % FN400
  timeMS = [300 500];
  roi = {'LAS','RAS'};
  %roi = {'LAS'};
  %roi = {'RAS'};
  chans = cat(2,elecGroupsNum{ismember(elecGroupsStr,roi)});
elseif strcmp(anaComponent,'PON')
  % P O/N
  timeMS = [500 800];
  %roi = {'LPS','RPS'};
  %roi = {'LPS'};
  roi = {'RPS'};
  chans = cat(2,elecGroupsNum{ismember(elecGroupsStr,roi)});
end

timeSamp = abs(prepost(1) * sampleRate) + [((timeMS(1) / (1000/sampleRate)) + 1) ((timeMS(2) / (1000/sampleRate)) + 1)];

% get the voltage for the correct region
volt_CR = [];
volt_H = [];
volt_HSC = [];
volt_HSI = [];
volt_TCR = [];
volt_TH = [];
volt_THSC = [];
volt_THSI = [];
for sub = 1:length(subjects)
  if ismember(subjects{sub},badSub)
    continue
  else
    volt_CR = [volt_CR mean(mean(mean(data_CR(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_HSC = [volt_HSC mean(mean(mean(data_HSC(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_HSI = [volt_HSI mean(mean(mean(data_HSI(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_H = [volt_H mean(mean(mean(data_H(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TCR = [volt_TCR mean(mean(mean(data_TCR(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_THSC = [volt_THSC mean(mean(mean(data_THSC(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_THSI = [volt_THSI mean(mean(mean(data_THSI(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TH = [volt_TH mean(mean(mean(data_TH(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
  end
end

% correlate accuracy with voltage differences

cfg_plot = [];
% set up how the lines will look
cfg_plot.linespec = 'ko';
cfg_plot.linewidth = 1.5;
cfg_plot.marksize = 8;
cfg_plot.markcolor = 'w';

if strcmp(anaComponent,'FN400')
  % FN400: item accuracy with H - CR
  fprintf('FN400\n');
  
  % Item
  
  % side
  testType = 'side_item';
  cond = {'H','CR'};
  x1 = d_item_side(~ismember(subjects,badSub));
  y1 = (volt_H - volt_CR);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Item d'': %s - %s',cond{1},cond{2}));
  xlabel('Item d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     Item d'' with H - CR: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_item';
  cond = {'QH','QCR'};
  x1 = d_item_q(~ismember(subjects,badSub));
  y1 = (volt_TH - volt_TCR);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Item d'': %s - %s',cond{1},cond{2}));
  xlabel('Item d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question: Item d'' with H - CR: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % Source
  
  % side
  testType = 'side_source';
  cond = {'HSC','HSI'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_HSC - volt_HSI);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Source d'': %s - %s',cond{1},cond{2}));
  xlabel('Source d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     Source d'' with HSC - HSI: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_source';
  cond = {'QHSC','QHSI'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_THSC - volt_THSI);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Source d'': %s - %s',cond{1},cond{2}));
  xlabel('Source d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question: Source d'' with HSC - HSI: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);

elseif strcmp(anaComponent,'PON')
  % P O/N: source accuracy with HSC - HSI
  fprintf('P O/N\n');
  
  % Source
  % side
  testType = 'side_source';
  cond = {'HSC','HSI'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_HSC - volt_HSI);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Source d'': %s - %s',cond{1},cond{2}));
  xlabel('Source d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     Source d'' with HSC - HSI: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_source';
  cond = {'QHSC','QHSI'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_THSC - volt_THSI);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('Source d'': %s - %s',cond{1},cond{2}));
  xlabel('Source d''');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question: Source d'' with HSC - HSI: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',rho(2,1),rho(2,1)^2,p(2,1),m,b);
end

%% get the trial counts for each subject to use in subsampling

% nTrial_CR = num_CR(~ismember(subjects,badSub));
% nTrial_HSC = num_HSC(~ismember(subjects,badSub));
% nTrial_HSI = num_HSI(~ismember(subjects,badSub));
% nTrial_H = num_H(~ismember(subjects,badSub));
% nTrial_TCR = num_TCR(~ismember(subjects,badSub));
% nTrial_THSC = num_THSC(~ismember(subjects,badSub));
% nTrial_THSI = num_THSI(~ismember(subjects,badSub));
% nTrial_TH = num_TH(~ismember(subjects,badSub));
% 
% %save('~/Desktop/trialCounts.mat','nTrial_*');
% 
% % nTrial_CR = num_CR(~ismember(subjects,badSub));
% % nTrial_HSC = num_HSC(~ismember(subjects,badSub));
% % nTrial_HSI = num_HSI(~ismember(subjects,badSub));
% % nTrial_H = num_H(~ismember(subjects,badSub));
% 
% % exp 1
% load('~/Desktop/trialCounts.mat');
% 
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% nExp3Sub = length(nTrial_CR);
% nExp1Sub = sum(~ismember(subjects,badSub));
% subsRand = randperm(min([nExp3Sub,nExp1Sub]));
% 
% ind_CR = nan(length(nTrial_CR),1);
% ind_H = nan(length(nTrial_H),1);
% ind_HSC = nan(length(nTrial_HSC),1);
% ind_HSI = nan(length(nTrial_HSI),1);
% counter = 1;
% for sub = 1:length(subsRand)
%   ind_CR(sub) = randperm(max([nTrial_CR(subsRand(sub)),num_CR(counter)]));
%   ind_CR(sub) = ind_CR(sub,1:min([nTrial_CR(subsRand(sub)),num_CR(counter)]));
%   
%   ind_H(sub) = randperm(max([nTrial_H(subsRand(sub)),num_H(counter)]));
%   ind_H(sub) = ind_H(sub,1:min([nTrial_H(subsRand(sub)),num_H(counter)]));
%   
%   ind_HSC(sub) = randperm(max([nTrial_HSC(subsRand(sub)),num_HSC(counter)]));
%   ind_HSC(sub) = ind_HSC(sub,1:min([nTrial_HSC(subsRand(sub)),num_HSC(counter)]));
%   
%   ind_HSI(sub) = randperm(max([nTrial_HSI(subsRand(sub)),num_HSI(counter)]));
%   ind_HSI(sub) = ind_HSI(sub,1:min([nTrial_HSI(subsRand(sub)),num_HSI(counter)]));
%   
%   counter = counter + 1;
% end
% 
% % exp 2
% load('~/Desktop/trialCounts.mat');
% 
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% nExp3Sub = length(nTrial_CR);
% nExp2Sub = sum(~ismember(subjects,badSub));
% subs = randperm(min([nExp3Sub,nExp1Sub]));
% 
% ind_CR = nan(length(nTrial_TCR),1);
% ind_H = nan(length(nTrial_TH),1);
% ind_HSC = nan(length(nTrial_THSC),1);
% ind_HSI = nan(length(nTrial_THSI),1);
% counter = 1;
% for sub = 1:length(subsRand)
%   ind_CR(sub) = randperm(max([nTrial_TCR(subsRand(sub)),num_CR(counter)]));
%   ind_CR(sub) = ind_CR(sub,1:min([nTrial_TCR(subsRand(sub)),num_CR(counter)]));
%   
%   ind_H(sub) = randperm(max([nTrial_TH(subsRand(sub)),num_H(counter)]));
%   ind_H(sub) = ind_H(sub,1:min([nTrial_TH(subsRand(sub)),num_H(counter)]));
%   
%   ind_HSC(sub) = randperm(max([nTrial_THSC(subsRand(sub)),num_HSC(counter)]));
%   ind_HSC(sub) = ind_HSC(sub,1:min([nTrial_THSC(subsRand(sub)),num_HSC(counter)]));
%   
%   ind_HSI(sub) = randperm(max([nTrial_THSI(subsRand(sub)),num_HSI(counter)]));
%   ind_HSI(sub) = ind_HSI(sub,1:min([nTrial_THSI(subsRand(sub)),num_HSI(counter)]));
%   
%   counter = counter + 1;
% end
% 
% % both
% for sub = 1:length(subjects)
%   cfg = [];
%   cfg.trials = ind_CR;
%   data_CR_subsample(sub) = ft_preprocessing(cfg,data_CR(sub));
%   cfg = [];
%   cfg.trials = ind_H;
%   data_H_subsample(sub) = ft_preprocessing(cfg,data_H(sub));
%   cfg = [];
%   cfg.trials = ind_HSC;
%   data_HSC_subsample(sub) = ft_preprocessing(cfg,data_HSC(sub));
%   cfg = [];
%   cfg.trials = ind_HSI;
%   data_HSI_subsample(sub) = ft_preprocessing(cfg,data_HSI(sub));
% end

%% descriptive statistics: ttest

% fieldtrip ttest

cfg = [];
roi = {'LAS','RAS'};
%roi = {'LAS'};
%roi = {'RAS'};
cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
cfg.latency     = [0.3 0.5];

% cfg = [];
% roi = {'LPS','RPS'};
% %roi = {'LPS'};
% %roi = {'RPS'};
% cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
% cfg.latency     = [0.5 0.8];

cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
cfg.parameter   = 'individual';
cfg.method      = 'analytic';
cfg.statistic   = 'depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';

numSub = length(subjects) - length(badSub);
cfg.design(1,1:2*numSub) = [ones(1,numSub) 2*ones(1,numSub)];
cfg.design(2,1:2*numSub) = [1:numSub 1:numSub];
cfg.ivar = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar = 2; % the 2nd row in cfg.design contains the subject number

chan = unique(str2double(strrep(cfg.channel,'e','')));

% side
timesel_CR = find(ga_CR.time >= cfg.latency(1) & ga_CR.time <= cfg.latency(2));
timesel_HSC = find(ga_HSC.time >= cfg.latency(1) & ga_HSC.time <= cfg.latency(2));
timesel_HSI = find(ga_HSI.time >= cfg.latency(1) & ga_HSI.time <= cfg.latency(2));
timesel_H = find(ga_H.time >= cfg.latency(1) & ga_H.time <= cfg.latency(2));
values_CR = mean(mean(ga_CR.individual(:,chan,timesel_CR),2),3);
values_HSC = mean(mean(ga_HSC.individual(:,chan,timesel_HSC),2),3);
values_HSI = mean(mean(ga_HSI.individual(:,chan,timesel_HSI),2),3);
values_H = mean(mean(ga_H.individual(:,chan,timesel_H),2),3);
sem_CR = std(values_CR)/sqrt(length(values_CR));
sem_HSC = std(values_HSC)/sqrt(length(values_HSC));
sem_HSI =  std(values_HSI)/sqrt(length(values_HSI));
sem_H =  std(values_H)/sqrt(length(values_H));

fprintf('\t%s\t%s\t%s\t%s\n','CR','HSC','HSI','H');
% ga
[mean(values_CR,1) mean(values_HSC,1) mean(values_HSI,1) mean(values_H,1)]
% sub avg
%[values_CR values_HSC values_HSI values_H]

% both
stat_ttest_HvsCR = ft_timelockstatistics(cfg,ga_H,ga_CR);
stat_ttest_HSCvsHSI = ft_timelockstatistics(cfg,ga_HSC,ga_HSI);
stat_ttest_HSCvsCR = ft_timelockstatistics(cfg,ga_HSC,ga_CR);
stat_ttest_HSIvsCR = ft_timelockstatistics(cfg,ga_HSI,ga_CR);

% both
fprintf('H   (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(values_H,1),sem_H,mean(values_CR,1),sem_CR,stat_ttest_HvsCR.df,stat_ttest_HvsCR.stat,stat_ttest_HvsCR.prob);
fprintf('HSC (M=%.3f; SEM=%.3f) vs HSI (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(values_HSC,1),sem_HSC,mean(values_HSI,1),sem_HSI,stat_ttest_HSCvsHSI.df,stat_ttest_HSCvsHSI.stat,stat_ttest_HSCvsHSI.prob);
fprintf('HSC (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(values_HSC,1),sem_HSC,mean(values_CR,1),sem_CR,stat_ttest_HSCvsCR.df,stat_ttest_HSCvsCR.stat,stat_ttest_HSCvsCR.prob);
fprintf('HSI (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(values_HSI,1),sem_HSI,mean(values_CR,1),sem_CR,stat_ttest_HSIvsCR.df,stat_ttest_HSIvsCR.stat,stat_ttest_HSIvsCR.prob);
fprintf('\n');

% plot to see the effect in each subject
figure
plot([values_HSC values_HSI]','o-');
xlim([0.5 2.5])
title(sprintf('%.1fs--%.1fs, %s',cfg.latency(1),cfg.latency(2),sprintf(repmat('%s ',1,length(roi)),roi{:})));
ylabel('Microvolts (\muV)');
set(gca,'XTickLabel',{'','HSC','','HSI',''})
figure
plot([values_H values_CR]','o-');
xlim([0.5 2.5])
title(sprintf('%.1fs--%.1fs, %s',cfg.latency(1),cfg.latency(2),sprintf(repmat('%s ',1,length(roi)),roi{:})));
ylabel('Microvolts (\muV)');
set(gca,'XTickLabel',{'','H','','CR',''})

% % matlab dependent samples ttest
% HSCminHSI = values_HSC - values_HSI;
% [h,p,ci,stats] = ttest(HSCminHSI, 0, 0.05); % H0: mean = 0, alpha 0.05

%% mean amplitude line plots

cfg_plot = [];

if ismember('LAS',roi) || ismember('RAS',roi)
  conditions = {{'H','CR'}, {'H-SC','H-SI','CR'}};
  cfg_plot.yminmax = [-5 -1];
elseif ismember('LPS',roi) || ismember('RPS',roi)
  conditions = {{'H-SC','H-SI','CR'}};
  cfg_plot.yminmax = [1 5];
end
%cfg_plot.yminmax = [floor(min([mean(values_HSC,1),mean(values_HSI,1),mean(values_CR,1)])) ceil(max([mean(values_HSC,1),mean(values_HSI,1),mean(values_CR,1)]))];

% set up how the lines will look
cfg_plot.linewidth = 2;
cfg_plot.marksize = 10;
%cfg_plot.linespec = 'k--o';
%cfg_plot.markcolor = 'w';
cfg_plot.side.linespec = 'k--o';
cfg_plot.task.linespec = 'k-o';
cfg_plot.side.markcolor = 'w';
cfg_plot.task.markcolor = 'k';
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;

for i = 1:length(conditions)
  cond = conditions{i};
  
  %%%%%%%%%%%%%%%%%%%
  % Side
  %%%%%%%%%%%%%%%%%%%
  
  testType = 'side';
  figure
  if length(cond) == 3
    % lines
    plot([mean(values_HSC,1),mean(values_HSI,1),mean(values_CR,1)],cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    % errorbars
    h = errorbar(1,mean(values_HSC,1),sem_HSC,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_HSI,1),sem_HSI,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(3,mean(values_CR,1),sem_CR,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    plot(1,mean(values_HSC,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(2,mean(values_HSI,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(3,mean(values_CR,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
  elseif length(cond) == 2
    % lines
    plot([mean(values_H,1),mean(values_CR,1)],cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    % errorbars
    h = errorbar(1,mean(values_H,1),sem_H,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_CR,1),sem_CR,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    plot(1,mean(values_H,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(2,mean(values_CR,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
  end
  % make it look good
  axis([.5 (length(cond) + .5) cfg_plot.yminmax(1) cfg_plot.yminmax(2)])
  xlabel('Condition');
  ylabel('Voltage (\muV)');
  set(gca,'XTick',[1:length(cond)])
  set(gca,'XTickLabel',cond)
  set(gca,'YTick',[cfg_plot.yminmax(1):.5:cfg_plot.yminmax(2)])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('line_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.latency(1)*1000,cfg.latency(2)*1000,figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  
  %%%%%%%%%%%%%%%%%%%
  % Task
  %%%%%%%%%%%%%%%%%%%
  
  testType = 'task';
  figure
  if length(cond) == 3
    % lines
    plot([mean(values_THSC,1),mean(values_THSI,1),mean(values_TCR,1)],cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    % errorbars
    h = errorbar(1,mean(values_THSC,1),sem_THSC,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_THSI,1),sem_THSI,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(3,mean(values_TCR,1),sem_TCR,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    plot(1,mean(values_THSC,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    plot(2,mean(values_THSI,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    plot(3,mean(values_TCR,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
  elseif length(cond) == 2
    % lines
    plot([mean(values_TH,1),mean(values_TCR,1)],cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    % errorbars
    h = errorbar(1,mean(values_TH,1),sem_TH,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_TCR,1),sem_TCR,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    plot(1,mean(values_TH,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    plot(2,mean(values_TCR,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
  end
  % make it look good
  axis([.5 (length(cond) + .5) cfg_plot.yminmax(1) cfg_plot.yminmax(2)])
  xlabel('Condition');
  ylabel('Voltage (\muV)');
  set(gca,'XTick',[1:length(cond)])
  set(gca,'XTickLabel',cond)
  set(gca,'YTick',[cfg_plot.yminmax(1):.5:cfg_plot.yminmax(2)])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('line_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.latency(1)*1000,cfg.latency(2)*1000,figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  
  %%%%%%%%%%%%%%%%%%%
  % Side+Task
  %%%%%%%%%%%%%%%%%%%
  
  figure
  if length(cond) == 3
    % lines
    plot([mean(values_HSC,1),mean(values_HSI,1),mean(values_CR,1)],cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    plot([mean(values_THSC,1),mean(values_THSI,1),mean(values_TCR,1)],cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth);
    % errorbars
    h = errorbar(1,mean(values_HSC,1),sem_HSC,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_HSI,1),sem_HSI,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(3,mean(values_CR,1),sem_CR,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(1,mean(values_THSC,1),sem_THSC,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_THSI,1),sem_THSI,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(3,mean(values_TCR,1),sem_TCR,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    h1 = plot(1,mean(values_HSC,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(2,mean(values_HSI,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(3,mean(values_CR,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    h2 = plot(1,mean(values_THSC,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    plot(2,mean(values_THSI,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    plot(3,mean(values_TCR,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    % legend
    legend([h1 h2],'Side','Question','Location','NorthEast');
  elseif length(cond) == 2
    % lines
    plot([mean(values_H,1),mean(values_CR,1)],cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    plot([mean(values_TH,1),mean(values_TCR,1)],cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth);
    % errorbars
    h = errorbar(1,mean(values_H,1),sem_H,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_CR,1),sem_CR,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(1,mean(values_TH,1),sem_TH,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_TCR,1),sem_TCR,cfg_plot.task.linespec,'LineWidth',cfg_plot.errwidth);
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
    h1 = plot(1,mean(values_H,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(2,mean(values_CR,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    h2 = plot(1,mean(values_TH,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    plot(2,mean(values_TCR,1),cfg_plot.task.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.task.markcolor);
    % legend
    legend([h1 h2],'Side','Question','Location','NorthEast');
  end
  % make it look good
  axis([.5 (length(cond) + .5) cfg_plot.yminmax(1) cfg_plot.yminmax(2)])
  xlabel('Condition');
  ylabel('Voltage (\muV)');
  set(gca,'XTick',[1:length(cond)])
  set(gca,'XTickLabel',cond)
  set(gca,'YTick',[cfg_plot.yminmax(1):.5:cfg_plot.yminmax(2)])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('line_ga_both_%s%s%d_%d.%s',sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.latency(1)*1000,cfg.latency(2)*1000,figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
end

%% descriptive statistics: ANOVA

%% 3-way ANOVA: Hemisphere x Condition x Test Type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FN400
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.hemi = {'LAS','RAS'};
%cfg.cond = {'CR','H'};
cfg.cond = {'CR','HSC','HSI'};
cfg.test = {'S','T'};
cfg.latency = [0.3 0.5];
cfg.alpha = 0.05;
cfg.showtable = 1;
cfg.printTable_tex = 1;
if length(cfg.hemi) > 2 || length(cfg.cond) > 2 || length(cfg.test) > 2
  cfg.calcGGHF = 1;
else
  cfg.calcGGHF = 0;
end

fprintf('====================== RMAOV33_mod =========================\n');
fprintf('IV1: Hemisphere (%d;%s), IV2: Condition (%d;%s), IV3: Test Type (%d;%s)\n',...
  length(cfg.hemi),sprintf(repmat(' %s',1,length(cfg.hemi)),cfg.hemi{:}),...
  length(cfg.cond),sprintf(repmat(' %s',1,length(cfg.cond)),cfg.cond{:}),...
  length(cfg.test),sprintf(repmat(' %s',1,length(cfg.test)),cfg.test{:}));

anovamat = [];
for h = 1:length(cfg.hemi)
  cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{cfg.hemi{h}})});
  chan = unique(str2double(strrep(cfg.channel,'e','')));
  
  for c = 1:length(cfg.cond)
    for t = 1:length(cfg.test)
      timesel.(cfg.test{t}).(cfg.cond{c}) = eval(sprintf('find(ga_%s%s.time >= cfg.latency(1) & ga_%s%s.time <= cfg.latency(2));',cfg.test{t},cfg.cond{c},cfg.test{t},cfg.cond{c}));
      values.(cfg.test{t}).(cfg.cond{c}) = eval(sprintf('mean(mean(ga_%s%s.individual(:,chan,timesel.%s.%s),2),3);',cfg.test{t},cfg.cond{c},cfg.test{t},cfg.cond{c}));
      
      for sub = 1:size(values.(cfg.test{t}).(cfg.cond{c}),1)
        % format for RMAOV33_mod: [data hemi cond test subNum]
        anovamat = [anovamat; values.(cfg.test{t}).(cfg.cond{c})(sub) h c t sub];
      end
    end
  end
end
p_vals = RMAOV33_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parietal Old/New
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.hemi = {'LPS','RPS'};
cfg.cond = {'CR','HSC','HSI'};
cfg.test = {'S','T'};
cfg.latency = [0.5 0.8];
cfg.alpha = 0.05;
cfg.showtable = 1;
cfg.printTable_tex = 1;
if length(cfg.hemi) > 2 || length(cfg.cond) > 2 || length(cfg.test) > 2
  cfg.calcGGHF = 1;
else
  cfg.calcGGHF = 0;
end

fprintf('====================== RMAOV33_mod =========================\n');
fprintf('IV1: Hemisphere (%d;%s), IV2: Condition (%d;%s), IV3: Test Type (%d;%s)\n',...
  length(cfg.hemi),sprintf(repmat(' %s',1,length(cfg.hemi)),cfg.hemi{:}),...
  length(cfg.cond),sprintf(repmat(' %s',1,length(cfg.cond)),cfg.cond{:}),...
  length(cfg.test),sprintf(repmat(' %s',1,length(cfg.test)),cfg.test{:}));

anovamat = [];
for h = 1:length(cfg.hemi)
  cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{cfg.hemi{h}})});
  chan = unique(str2double(strrep(cfg.channel,'e','')));
  
  for c = 1:length(cfg.cond)
    for t = 1:length(cfg.test)
      timesel.(cfg.test{t}).(cfg.cond{c}) = eval(sprintf('find(ga_%s%s.time >= cfg.latency(1) & ga_%s%s.time <= cfg.latency(2));',cfg.test{t},cfg.cond{c},cfg.test{t},cfg.cond{c}));
      values.(cfg.test{t}).(cfg.cond{c}) = eval(sprintf('mean(mean(ga_%s%s.individual(:,chan,timesel.%s.%s),2),3);',cfg.test{t},cfg.cond{c},cfg.test{t},cfg.cond{c}));
      
      for sub = 1:size(values.(cfg.test{t}).(cfg.cond{c}),1)
        % format for RMAOV33_mod: [data hemi cond test subNum]
        anovamat = [anovamat; values.(cfg.test{t}).(cfg.cond{c})(sub) h c t sub];
      end
    end
  end
end
p_vals = RMAOV33_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

%% 2-way ANOVA: Hemisphere x Condition for each trial type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FN400
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.hemi = {'LAS','RAS'};
%cfg.cond = {'CR','H'};
cfg.cond = {'CR','HSC','HSI'};
cfg.latency = [0.3 0.5];
cfg.alpha = 0.05;
cfg.showtable = 1;
cfg.printTable_tex = 1;
if length(cfg.hemi) > 2 || length(cfg.cond) > 2
  cfg.calcGGHF = 1;
else
  cfg.calcGGHF = 0;
end

fprintf('====================== RMAOV2_mod =========================\n');
fprintf('IV1: Hemisphere (%d;%s), IV2: Condition (%d;%s)\n',length(cfg.hemi),sprintf(repmat(' %s',1,length(cfg.hemi)),cfg.hemi{:}),length(cfg.cond),sprintf(repmat(' %s',1,length(cfg.cond)),cfg.cond{:}));

anovamat = [];
for h = 1:length(cfg.hemi)
  cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{cfg.hemi{h}})});
  chan = unique(str2double(strrep(cfg.channel,'e','')));
  
  for c = 1:length(cfg.cond)
    timesel.(cfg.cond{c}) = eval(sprintf('find(ga_%s.time >= cfg.latency(1) & ga_%s.time <= cfg.latency(2));',cfg.cond{c},cfg.cond{c}));
    values.(cfg.cond{c}) = eval(sprintf('mean(mean(ga_%s.individual(:,chan,timesel.%s),2),3);',cfg.cond{c},cfg.cond{c}));
    
    %index = h + (c - 1) + (h - 1);
    for sub = 1:size(values.(cfg.cond{c}),1)
      
      % format for RMAOV2_mod: [data hemi cond subNum]
      anovamat = [anovamat; values.(cfg.cond{c})(sub) h c sub];
    end
  end
end
p_vals = RMAOV2_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P O/N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.hemi = {'LPS','RPS'};
cfg.cond = {'CR','HSC','HSI'};
cfg.latency = [0.5 0.8];
cfg.alpha = 0.05;
cfg.showtable = 1;
cfg.printTable_tex = 1;
if length(cfg.hemi) > 2 || length(cfg.cond) > 2
  cfg.calcGGHF = 1;
else
  cfg.calcGGHF = 0;
end

fprintf('====================== RMAOV2_mod =========================\n');
fprintf('IV1: Hemisphere (%d;%s), IV2: Condition (%d;%s)\n',length(cfg.hemi),sprintf(repmat(' %s',1,length(cfg.hemi)),cfg.hemi{:}),length(cfg.cond),sprintf(repmat(' %s',1,length(cfg.cond)),cfg.cond{:}));

anovamat = [];
for h = 1:length(cfg.hemi)
  cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{cfg.hemi{h}})});
  chan = unique(str2double(strrep(cfg.channel,'e','')));
  
  for c = 1:length(cfg.cond)
    timesel.(cfg.cond{c}) = eval(sprintf('find(ga_%s.time >= cfg.latency(1) & ga_%s.time <= cfg.latency(2));',cfg.cond{c},cfg.cond{c}));
    values.(cfg.cond{c}) = eval(sprintf('mean(mean(ga_%s.individual(:,chan,timesel.%s),2),3);',cfg.cond{c},cfg.cond{c}));
    
    %index = h + (c - 1) + (h - 1);
    for sub = 1:size(values.(cfg.cond{c}),1)
      
      % format for RMAOV2_mod: [data hemi cond subNum]
      anovamat = [anovamat; values.(cfg.cond{c})(sub) h c sub];
    end
  end
end
p_vals = RMAOV2_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

%% cluster statistics

cfg = [];
%start and stop time of analysis, in sec
cfg.latency = [0 1.0];
cfg.avgovertime = 'no';

cfg.channel = 'all';
%cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'LAS','RAS'})});
%cfg.avgoverchan = 'yes';

cfg.parameter = 'individual';

% use the Monte Carlo Method to calculate the significance probability
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
% alpha level of the sample-specific test statistic that will be used for
% thresholding
cfg.clusteralpha = 0.05;
% test statistic that will be evaluated under the permutation distribution
cfg.clusterstatistic = 'maxsum';
% minimum number of neighborhood channels that is required for a selected
% sample to be included in the clustering algorithm (default = 0)
cfg.minnbchan = 2;
% alpha level of the permutation test
cfg.alpha = 0.025;
% number of draws from the permutation distribution, should be 1000
cfg.numrandomization = 1000;

% set the unit and independent variables
cfg.uvar = 1; % row of design matrix containing subject numbers (DV)
cfg.ivar = 2; % row of design matrix containing conditions (IV)

numConds = 2;
% use the dependent samples T-statistic as a measure to evaluate the
% effect at the sample level
if numConds == 2
  cfg.statistic = 'depsamplesT';
  % test tails: -1 = left tail, 0 = two tail, 1 = right tail?
  cfg.tail = 0;
  cfg.clustertail = 0;
elseif numConds > 2
  cfg.statistic = 'depsamplesF';
  % test tails: -1 = left tail, 0 = two tail, 1 = right tail?
  cfg.tail = 1;
  cfg.clustertail = 1;
end

% make the design matrix
numSub = length(subjects) - length(badSub);
design = zeros(2,numSub*numConds);
for i = 1:numConds
  for j = 1:numSub
    design(1,1+((i - 1)*numSub) + (j - 1)) = j; % subject #s
    design(2,1+((i - 1)*numSub) + (j - 1)) = i; % condition #s
  end
end
cfg.design = design;

% side
stat_clus_HSCvsHSI = ft_timelockstatistics(cfg,ga_HSC,ga_HSI);
stat_clus_HSCvsCR = ft_timelockstatistics(cfg,ga_HSC,ga_CR);
stat_clus_HSIvsCR = ft_timelockstatistics(cfg,ga_HSI,ga_CR);

% task
stat_clus_THSCvsTHSI = ft_timelockstatistics(cfg,ga_THSC,ga_THSI);
stat_clus_THSCvsTCR = ft_timelockstatistics(cfg,ga_THSC,ga_TCR);
stat_clus_THSIvsTCR = ft_timelockstatistics(cfg,ga_THSI,ga_TCR);

% side vs task
stat_clus_HSCvsTHSC = ft_timelockstatistics(cfg,ga_HSC,ga_THSC);
stat_clus_HSIvsTHSI = ft_timelockstatistics(cfg,ga_HSI,ga_THSI);
stat_clus_CRvsTCR = ft_timelockstatistics(cfg,ga_CR,ga_TCR);

%% plot the cluster statistics

cfg = [];
% p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
cfg.highlightsymbolseries = ['*','+','.','.','.'];
cfg.highlightcolorpos = [0.5 0 1];
cfg.highlightcolorneg = [0 0.5 0];
cfg.layout = ft_prepare_layout(cfg,ga_HSC);
cfg.contournum = 0;
cfg.emarker = '.';
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.zlim = [-5 5];
%cfg.interplimits = 'electrodes';
ft_clusterplot(cfg,stat_clus_HSCvsHSI); % late FN400? P O/N
ft_clusterplot(cfg,stat_clus_HSCvsCR); % P O/N, late pos front
ft_clusterplot(cfg,stat_clus_HSIvsCR); % late P O/N dipole? LPS not sig

ft_clusterplot(cfg,stat_clus_THSCvsTHSI); % P O/N
ft_clusterplot(cfg,stat_clus_THSCvsTCR); % P O/N, late pos front
ft_clusterplot(cfg,stat_clus_THSIvsTCR); % late pos front dipole?

ft_clusterplot(cfg,stat_clus_HSCvsTHSC); % none
ft_clusterplot(cfg,stat_clus_HSIvsTHSI); % none
ft_clusterplot(cfg,stat_clus_CRvsTCR); % none

%% older stuff

% create contrast
ga_HSCvsHSI = ga_HSC;
ga_HSCvsHSI.avg = ga_HSC.avg - ga_HSI.avg;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
pos = stat_clus_HSCvsHSI.posclusterslabelmat==1;
for k = 1:length(j)-1;
  subplot(4,5,k);
  cfg = [];
  cfg.xlim=[j(k) j(k+1)];
  cfg.zlim = [-5 5];
  pos_int = mean(pos(:,m(k):m(k+1))')'; % mean(pos(:,m(k):m(k+1)),2);
  cfg.highlight = 'on';
  cfg.highlightchannel = find(pos_int==1);
  cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  ft_topoplotER(cfg,ga_HSCvsHSI);
end

% create contrast
ga_HSCvsCR = ga_HSC;
ga_HSCvsCR.avg = ga_HSC.avg - ga_CR.avg;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
pos = stat_clus_HSCvsCR.posclusterslabelmat==1;
for k = 1:length(j)-1;
  subplot(4,5,k);
  cfg = [];
  cfg.xlim=[j(k) j(k+1)];
  cfg.zlim = [-5 5];
  pos_int = mean(pos(:,m(k):m(k+1))')';
  cfg.highlight = 'on';
  cfg.highlightchannel = find(pos_int==1);
  cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  ft_topoplotER(cfg,ga_HSCvsCR);
end

% create contrast
ga_HSIvsCR = ga_HSI;
ga_HSIvsCR.avg = ga_HSI.avg - ga_CR.avg;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
pos = stat_clus_HSIvsCR.posclusterslabelmat==1;
for k = 1:length(j)-1;
  subplot(4,5,k);
  cfg = [];
  cfg.xlim=[j(k) j(k+1)];
  cfg.zlim = [-5 5];
  pos_int = mean(pos(:,m(k):m(k+1))')';
  cfg.highlight = 'on';
  cfg.highlightchannel = find(pos_int==1);
  cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  ft_topoplotER(cfg,ga_HSIvsCR);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% because of a bug (might be fixed now)
if ~isfield(stat_HSCvsHSIvsCR,'negclusters') && isfield(stat_HSCvsHSIvsCR,'posclusters')
  fprintf('No neg clusters found\n');
  stat_HSCvsHSIvsCR.negclusters.prob = .5;
  stat_HSCvsHSIvsCR.negclusters.clusterstat = 0;
  stat_HSCvsHSIvsCR.negclusterslabelmat = zeros(size(stat_HSCvsHSIvsCR.posclusterslabelmat));
  stat_HSCvsHSIvsCR.negdistribution = zeros(size(stat_HSCvsHSIvsCR.posdistribution));
end
if ~isfield(stat_HSCvsHSIvsCR,'posclusters') && isfield(stat_HSCvsHSIvsCR,'negclusters')
  fprintf('No pos clusters found\n');
  stat_HSCvsHSIvsCR.posclusters.prob = 1;
  stat_HSCvsHSIvsCR.posclusters.clusterstat = 0;
  stat_HSCvsHSIvsCR.posclusterslabelmat = zeros(size(stat_HSCvsHSIvsCR.negclusterslabelmat));
  stat_HSCvsHSIvsCR.posdistribution = zeros(size(stat_HSCvsHSIvsCR.negdistribution));
end

cfg = [];
% p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
cfg.highlightsymbolseries = ['*','*','.','.','.'];
cfg.layout = ft_prepare_layout(cfg,ga_HSC);
cfg.contournum = 0;
cfg.emarker = '.';
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.zlim = [-5 5];
ft_clusterplot(cfg,stat_HSCvsHSIvsCR);
