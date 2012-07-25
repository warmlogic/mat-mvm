% Do a nonparametric statistics clustering analysis for timelocked EEG
% (ERPs)

% See Maris and Oostenveld (2007) for a description of statistics

%% Experiment-specific setup

expName = 'SLORK';

% pre- and post-stimulus times to save, in seconds
prepost = [-1.0 2.0];

sampleRate = 250;

%eventValues = sort({'CR','KNOW','RO','RS'});

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
saveDir = fullfile(dataroot,'ft_data_st');
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
    data_preproc = seg2ft_nsmat(fullfile(dataroot,'ns_mat_st'),subjects{sub},prepost,sampleRate,elecfile);
    
    fprintf('Saving %s FieldTrip struct...\n',subjects{sub});
    save(dataSaveFile,'data_preproc','prepost','elec');
  end
  
  % add NS's artifact information to the event structure
  %ns_addArtifactInfo_slork(dataroot,subjects{sub});

  fprintf('Done with %s.\n',subjects{sub});
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% data_avg = ft_timelockanalysis(cfg,data_F(1));
% 
% cfg = [];
% cfg.showlabels = 'yes'; 
% cfg.fontsize = 9; 
% cfg.layout = data_avg.elec;
% %cfg.ylim = [-3e-13 3e-13];
% figure
% ft_multiplotER(cfg,data_avg); 

%% Run ft_timelockanalysis; make sure events go in alphabetical order

%clear data_N data_F data_RS

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
num_SN = zeros(length(subjects),1);
num_SF = zeros(length(subjects),1);
num_SRO = zeros(length(subjects),1);
num_SRS = zeros(length(subjects),1);
num_TN = zeros(length(subjects),1);
num_TF = zeros(length(subjects),1);
num_TRO = zeros(length(subjects),1);
num_TRS = zeros(length(subjects),1);

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
  
  %eval(sprintf('clear data_*_raw'));
  clear data_*_pp
  
  fprintf('Done.\n');
end

%matlabpool close

%% save the ft_timelockanalysis files for each event

tlaFilename = sprintf('data_SN_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_SN','num_SN');

tlaFilename = sprintf('data_SF_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_SF','num_SF');

tlaFilename = sprintf('data_SRO_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_SRO','num_SRO');

tlaFilename = sprintf('data_SRS_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_SRS','num_SRS');

tlaFilename = sprintf('data_TN_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_TN','num_TN');

tlaFilename = sprintf('data_TF_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_TF','num_TF');

tlaFilename = sprintf('data_TRO_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_TRO','num_TRO');

tlaFilename = sprintf('data_TRS_tla_%d_%d.mat',prepost(1)*1000,prepost(2)*1000);
tlaSaveFile = fullfile(saveDir,tlaFilename);
save(tlaSaveFile,'data_TRS','num_TRS');

%% if already saved, load the ft_timelockanalysis files

savedFiles = dir(fullfile(saveDir,'data_*_tla_*.mat'));
for i = 1:length(savedFiles)
  fprintf('Loading %s...',savedFiles(i).name);
  load(fullfile(saveDir,savedFiles(i).name));
  fprintf('Done.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove the degrees of freedom? field. I'm not sure what it does.
data_SN = rmfield(data_SN,'dof');
data_SF = rmfield(data_SF,'dof');
data_SRO = rmfield(data_SRO,'dof');
data_SRS = rmfield(data_SRS,'dof');

data_TN = rmfield(data_TN,'dof');
data_TF = rmfield(data_TF,'dof');
data_TRO = rmfield(data_TRO,'dof');
data_TRS = rmfield(data_TRS,'dof');

%% decide who to kick out based on trial counts
[num_SN num_SF num_SRO num_SRS]
[num_TN num_TF num_TRO num_TRS]

lowNum = [];
thresh = 5;
% get the number of trials
for sub = 1:length(subjects)
  if num_SN(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_SF(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_SRO(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_SRS(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_TN(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_TF(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_TRO(sub) < thresh
    lowNum = [lowNum sub];
  end
  if num_TRS(sub) < thresh
    lowNum = [lowNum sub];
  end
end

% who we're rejecting
subjects(unique(lowNum))

% should be: 'SLORK003'    'SLORK013'    'SLORK018'    'SLORK029'

% noisy frontal EEG: 'SLORK012', 'SLORK022', 'SLORK025'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

badBehSub = {};

% item d' < 1.0
%badBehSub = {'SLORK010','SLORK011','SLORK030'};
%badBehSub = {'SLORK011','SLORK017'};

% SLORK012 and SLORK030's EEG is really noisy in LAS and RAS
%badBehSub = {'SLORK012','SLORK030'};

badSub = unique([badBehSub subjects(lowNum)]);

fprintf('Number of events included in EEG analyses (%d subjects; threshold: %d events):\n',sum(~ismember(subjects,badSub)),thresh);
sn_mean = mean(num_SN(~ismember(subjects,badSub)));
sn_sd = std(num_SN(~ismember(subjects,badSub)),0,1);
sn_sem = std(num_SN(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('SN: M=%.3f, SD=%.3f, SEM=%.3f\n',sn_mean,sn_sd,sn_sem);
sf_mean = mean(num_SF(~ismember(subjects,badSub)));
sf_sd = std(num_SF(~ismember(subjects,badSub)),0,1);
sf_sem = std(num_SF(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('SF: M=%.3f, SD=%.3f, SEM = %.3f\n',sf_mean,sf_sd,sf_sem);
sro_mean = mean(num_SRO(~ismember(subjects,badSub)));
sro_sd = std(num_SRO(~ismember(subjects,badSub)),0,1);
sro_sem = std(num_SRO(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('SRO: M=%.3f, SD=%.3f, SEM=%.3f\n',sro_mean,sro_sd,sro_sem);
srs_mean = mean(num_SRS(~ismember(subjects,badSub)));
srs_sd = std(num_SRS(~ismember(subjects,badSub)),0,1);
srs_sem = std(num_SRS(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('SRS: M=%.3f, SD=%.3f, SEM=%.3f\n',srs_mean,srs_sd,srs_sem);

tn_mean = mean(num_TN(~ismember(subjects,badSub)));
tn_sd = std(num_TN(~ismember(subjects,badSub)),0,1);
tn_sem = std(num_TN(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('TN: M=%.3f, SD=%.3f, SEM=%.3f\n',tn_mean,tn_sd,tn_sem);
tf_mean = mean(num_TF(~ismember(subjects,badSub)));
tf_sd = std(num_TF(~ismember(subjects,badSub)),0,1);
tf_sem = std(num_TF(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('TF: M=%.3f, SD=%.3f, SEM = %.3f\n',tf_mean,tf_sd,tf_sem);
tro_mean = mean(num_TRO(~ismember(subjects,badSub)));
tro_sd = std(num_TRO(~ismember(subjects,badSub)),0,1);
tro_sem = std(num_TRO(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('TRO: M=%.3f, SD=%.3f, SEM=%.3f\n',tro_mean,tro_sd,tro_sem);
trs_mean = mean(num_TRS(~ismember(subjects,badSub)));
trs_sd = std(num_TRS(~ismember(subjects,badSub)),0,1);
trs_sem = std(num_TRS(~ismember(subjects,badSub)),0,1)/sqrt(sum(~ismember(subjects,badSub)));
fprintf('TRS: M=%.3f, SD=%.3f, SEM=%.3f\n',trs_mean,trs_sd,trs_sem);

%% set up analysis strings

str_SN = [];
str_SF = [];
str_SRO = [];
str_SRS = [];
str_TN = [];
str_TF = [];
str_TRO = [];
str_TRS = [];

for sub = 1:length(subjects)
  if ismember(subjects{sub},badSub)
    continue
  else
    if sub == 1
      str_SN = sprintf('data_SN(%d)',sub);
      str_SF = sprintf('data_SF(%d)',sub);
      str_SRO = sprintf('data_SRO(%d)',sub);
      str_SRS = sprintf('data_SRS(%d)',sub);
      str_TN = sprintf('data_TN(%d)',sub);
      str_TF = sprintf('data_TF(%d)',sub);
      str_TRO = sprintf('data_TRO(%d)',sub);
      str_TRS = sprintf('data_TRS(%d)',sub);
    elseif sub > 1
      str_SN = sprintf('%s,data_SN(%d)',str_SN,sub);
      str_SF = sprintf('%s,data_SF(%d)',str_SF,sub);
      str_SRO = sprintf('%s,data_SRO(%d)',str_SRS,sub);
      str_SRS = sprintf('%s,data_SRS(%d)',str_SRS,sub);
      str_TN = sprintf('%s,data_TN(%d)',str_TN,sub);
      str_TF = sprintf('%s,data_TF(%d)',str_TF,sub);
      str_TRO = sprintf('%s,data_TRO(%d)',str_TRS,sub);
      str_TRS = sprintf('%s,data_TRS(%d)',str_TRS,sub);
    end
  end
end
if strcmp(str_SF(1),',')
  str_SN = str_SN(2:end);
  str_SF = str_SF(2:end);
  str_SRO = str_SRO(2:end);
  str_SRS = str_SRS(2:end);
  str_TN = str_TN(2:end);
  str_TF = str_TF(2:end);
  str_TRO = str_TRO(2:end);
  str_TRS = str_TRS(2:end);
end

%% get the grand average and add the electrode locations again

cfg = [];
cfg.keepindividual = 'yes';
eval(sprintf('ga_SN = ft_timelockgrandaverage(cfg,%s);',str_SN));
eval(sprintf('ga_SF = ft_timelockgrandaverage(cfg,%s);',str_SF));
eval(sprintf('ga_SRO = ft_timelockgrandaverage(cfg,%s);',str_SRO));
eval(sprintf('ga_SRS = ft_timelockgrandaverage(cfg,%s);',str_SRS));
eval(sprintf('ga_TN = ft_timelockgrandaverage(cfg,%s);',str_TN));
eval(sprintf('ga_TF = ft_timelockgrandaverage(cfg,%s);',str_TF));
eval(sprintf('ga_TRO = ft_timelockgrandaverage(cfg,%s);',str_TRO));
eval(sprintf('ga_TRS = ft_timelockgrandaverage(cfg,%s);',str_TRS));

ga_SN.elec = elec;
ga_SF.elec = elec;
ga_SRO.elec = elec;
ga_SRS.elec = elec;
ga_TN.elec = elec;
ga_TF.elec = elec;
ga_TRO.elec = elec;
ga_TRS.elec = elec;

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
% cfg.layout = ft_prepare_layout(cfg,ga_N);
% figure
% ft_multiplotER(cfg,ga_N,ga_F,ga_RS);
% figure
% ft_multiplotER(cfg,ga_TCR,ga_THSC,ga_THSI);

%cfg.interplimits = 'electrodes';
%ft_topoplotER(cfg,ga_N);

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
ft_singleplotER(cfg,ga_N,ga_RS,ga_F);
legend('N','RSrc','F','Location','SouthEast');
% P O/N
cfg.ylim = [-2 5];
cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'LPS','RPS'})});
figure
ft_singleplotER(cfg,ga_N,ga_RS,ga_F);
legend('N','RSrc','F','Location','NorthEast');

%% subplot with all subjects' ERPs

%roi = {'LAS','RAS'};
roi = {'LPS','RPS'};
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

% both
testType = 'both_rkn';
figure
for sub = 1:length(subjects)
  subplot(numRows,numCols,sub);
  plot(data_N(sub).time,mean(data_N(sub).avg(chan,:),1),'k');
  hold on
  if ~isempty(data_F(sub).avg)
    plot(data_F(sub).time,mean(data_F(sub).avg(chan,:),1),'b');
  end
  plot(data_RS(sub).time,mean(data_RS(sub).avg(chan,:),1),'r');
  if ismember(subjects{sub},badSub)
    title(sprintf('*BAD*: %s; N:%d, F:%d, RSrc:%d',subjects{sub},size(data_N(sub).trial,1),size(data_F(sub).trial,1),size(data_RS(sub).trial,1)));
  else
    title(sprintf('%s; N:%d, F:%d, RSrc:%d',subjects{sub},size(data_N(sub).trial,1),size(data_F(sub).trial,1),size(data_RS(sub).trial,1)));
  end
  %ylim([0 1.9e-13]);
  axis([-0.2 1.5 -10 10]);
  hold off
end
set(gcf,'Name',[testType,': ',sprintf(repmat('%s ',1,length(roi)),roi{:})])
subplot(numRows,numCols,length(subjects) + 1);
text(0.5,0.9,[testType,', ',sprintf(repmat('%s ',1,length(roi)),roi{:})],'color','k')
text(0.5,0.7,'N','color','k')
text(0.5,0.5,'F','color','r') ;
text(0.5,0.3,'RSrc','color','b')
axis off

%% plot the conditions

cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 9;
cfg.layout = ft_prepare_layout(cfg,ga_N);
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
testType = 'both_rkn';
for i = 1:length(plot_rois)
  roi = plot_rois{i};
  cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
  figure
  %hold on
  ft_singleplotER(cfg,ga_N,ga_RS,ga_F);
  hold on
  plot([latency(1) latency(2)],[0 0],'k--'); % horizontal
  if ismember('LAS',roi) || ismember('RAS',roi)
    plot([0 0],[minMaxVolt(1,1) minMaxVolt(1,2)],'k--'); % vertical
    if plotLegend
      legend('N','RSrc','F','Location','SouthEast');
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
      legend('N','RSrc','F','Location','NorthEast');
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

%% Plot average time contrasts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both: LAS + RAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.interactive = 'no';
minMaxVolt = [-1 1];
%cfg.colorbar = 'yes';
cfg.layout = ft_prepare_layout(cfg,ga_F);
%cfg.marker = 'on';
%cfg.marker = 'numbers';
cfg.highlight = 'on';
cfg.highlightsize = 10;
%cfg.comment = 'xlim';
%cfg.commentpos = 'title';
cfg.comment = 'no';
plotTitle = 0;
plotColorbar = 0;
if plotColorbar
  colorbar_str = '_cb';
else
  colorbar_str = '';
end

cfg.colormap = colormap('jet'); % default; blue to red
%cfg.colormap = colormap('hot'); % dark to light; better for b&w printers

testType = 'both_rkn';

% create contrast
ga_FvsCR = ga_F;
ga_FvsCR.avg = ga_F.avg - ga_N.avg;
ga_FvsCR.individual = ga_F.individual - ga_N.individual;
% make a plot
figure
roi = {'LAS','RAS'};
cond = {'F','N'};
cfg.xlim = [0.3 0.5];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_FvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('Know - Correct Rejections');
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_RSvsCR = ga_RS;
ga_RSvsCR.avg = ga_RS.avg - ga_N.avg;
ga_RSvsCR.individual = ga_RS.individual - ga_N.individual;
% make a plot
figure
roi = {'LAS','RAS'};
cond = {'RSrc','N'};
cfg.xlim = [0.3 0.5];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_RSvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('Remember Source - Correct Rejections');
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_RSvsKNOW = ga_RS;
ga_RSvsKNOW.avg = ga_RS.avg - ga_F.avg;
ga_RSvsKNOW.individual = ga_RS.individual - ga_F.individual;
% make a plot
figure
roi = {'LAS','RAS'};
cond = {'RSrc','F'};
cfg.xlim = [0.3 0.5];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_RSvsKNOW);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('Remember Source - Know');
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both: LPS + RPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create contrast
ga_FvsCR = ga_F;
ga_FvsCR.avg = ga_F.avg - ga_N.avg;
ga_FvsCR.individual = ga_F.individual - ga_N.individual;
% make a plot
figure
roi = {'LPS','RPS'};
cond = {'F','N'};
cfg.xlim = [0.5 0.8];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_FvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('Know - Correct Rejections');
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_RSvsCR = ga_RS;
ga_RSvsCR.avg = ga_RS.avg - ga_N.avg;
ga_RSvsCR.individual = ga_RS.individual - ga_N.individual;
% make a plot
figure
roi = {'LPS','RPS'};
cond = {'RSrc','N'};
cfg.xlim = [0.5 0.8];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_RSvsCR);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('Remember Source - Correct Rejections');
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

% create contrast
ga_RSvsKNOW = ga_RS;
ga_RSvsKNOW.avg = ga_RS.avg - ga_F.avg;
ga_RSvsKNOW.individual = ga_RS.individual - ga_F.individual;
% make a plot
figure
roi = {'LPS','RPS'};
cond = {'RSrc','F'};
cfg.xlim = [0.5 0.8];
cfg.zlim = [minMaxVolt(1) minMaxVolt(2)];
cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
ft_topoplotER(cfg,ga_RSvsKNOW);
set(gcf,'Name',sprintf('%s - %s',cond{1},cond{2}))
if plotColorbar
  h = colorbar;
  set(get(h,'YLabel'),'string','Voltage (\muV)');
end
if plotTitle
  titleStr = sprintf('Remember Source - Know');
  title(titleStr);
end
publishfig(gcf,0);
figfilename = sprintf('topo_ga_%s_%s%s%d_%d%s.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),cfg.xlim(1)*1000,cfg.xlim(2)*1000,colorbar_str,figFileExt);
if saveFigs
  print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
end

%% correlations dp

% get the electrode numbers
elecGroupsNum = cell(size(elecGroups));
for i = 1:length(elecGroups)
  elecGroupsNum{i} = str2double(strrep(elecGroups{i},'e',''));
end

% accuracy
d_item_side = abs([1.6145, 2.0823, 2.6991, 1.7986, 1.9783, 1.0971, 2.0902, 2.0895, 0.6578, 0.6442, 1.678, 1.7853, 1.6032, 2.4474, 2.3005, 1.7255, 1.566, 1.2529, 1.0774, 3.0318, 1.9626, 1.5029, 0.8762, 2.7775, 1.553, 1.3498, 0.8494, 0.844, 0.7109, 1.2076]);
d_source_side = abs([1.5528, 1.7537, 1.9537, 1.7535, 1.687, 0.0788, 1.2857, 2.0013, 0.4771, 0.6571, 1.1276, 1.8787, 1.1546, 3.0483, 1.4379, 0.7955, 1.7102, 0.8845, 1.1945, 1.7385, 1.7149, 2.3272, 0.0968, 3.2844, 1.2291, 0.9366, -0.6108, 0.6999, 0.2037, 0.622]);
d_item_q = abs([1.8692, 2.0083, 3.0126, 1.1594, 1.6155, 1.1885, 1.6652, 1.7839, 0.9898, 1.0621, 1.5352, 2.1099, 1.4867, 2.6831, 2.3211, 1.7398, 1.7557, 1.5432, 1.4436, 2.4062, 2.2162, 1.7094, 1.39, 2.8908, 1.301, 1.0106, 1.0331, 1.2649, 1.6038, 1.5393]);
d_source_q = abs([0.7439, 1.0587, 1.5439, 0.7418, 0.8039, -0.3197, 1.0668, 1.317, 0.0736, 0.5171, 1.3582, 1.2274, 0.738, 1.6176, 2.1666, 0.5236, 1.4651, 1.1363, 0.3559, 0.8078, 1.0654, 1.8929, -0.2507, 1.3646, 1.2062, 0.3971, -0.1053, 1.2216, 0.511, 0.8607]);

anaComponent = 'FN400';
%anaComponent = 'PON';

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
  roi = {'LPS','RPS'};
  %roi = {'LPS'};
  %roi = {'RPS'};
  chans = cat(2,elecGroupsNum{ismember(elecGroupsStr,roi)});
end

timeSamp = abs(prepost(1) * sampleRate) + [((timeMS(1) / (1000/sampleRate)) + 1) ((timeMS(2) / (1000/sampleRate)) + 1)];

% get the voltage for the correct region
volt_SN = [];
volt_SF = [];
volt_SRO = [];
volt_SRS = [];
volt_TN = [];
volt_TF = [];
volt_TRO = [];
volt_TRS = [];
for sub = 1:length(subjects)
  if ismember(subjects{sub},badSub)
    continue
  else
    volt_SN = [volt_SN mean(mean(mean(data_SN(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_SF = [volt_SF mean(mean(mean(data_SF(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_SRO = [volt_SRO mean(mean(mean(data_SRO(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_SRS = [volt_SRS mean(mean(mean(data_SRS(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TN = [volt_TN mean(mean(mean(data_TN(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TF = [volt_TF mean(mean(mean(data_TF(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TRO = [volt_TRO mean(mean(mean(data_TRO(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TRS = [volt_TRS mean(mean(mean(data_TRS(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
  end
end

% correlate accuracy with voltage differences

cfg_plot = [];
% set up how the lines will look
cfg_plot.linespec = 'ko';
cfg_plot.linewidth = 1.5;
cfg_plot.marksize = 8;
cfg_plot.markcolor = 'w';

% RSrc - N
% RSrc - F
% F - N

if strcmp(anaComponent,'FN400')
  % FN400: item accuracy with RSrc - N
  fprintf('FN400\n');
  
  % Item
  
  % side
  testType = 'side_itemdp';
  cond = {'RSrc','N'};
  x1 = d_item_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_itemdp';
  cond = {'RSrc','F'};
  x1 = d_item_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_itemdp';
  cond = {'F','N'};
  x1 = d_item_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_itemdp';
  cond = {'RSrc','N'};
  x1 = d_item_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_itemdp';
  cond = {'RSrc','F'};
  x1 = d_item_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_itemdp';
  cond = {'F','N'};
  x1 = d_item_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % Source
  
  % side
  testType = 'side_sourcedp';
  cond = {'RSrc','N'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_sourcedp';
  cond = {'RSrc','F'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_sourcedp';
  cond = {'F','N'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_sourcedp';
  cond = {'RSrc','N'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_sourcedp';
  cond = {'RSrc','F'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_sourcedp';
  cond = {'F','N'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  

elseif strcmp(anaComponent,'PON')
  % P O/N: source accuracy with HSC - HSI
  fprintf('P O/N\n');
  
  % Item
  
  % side
  testType = 'side_itemdp';
  cond = {'RSrc','N'};
  x1 = d_item_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_itemdp';
  cond = {'RSrc','F'};
  x1 = d_item_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_itemdp';
  cond = {'F','N'};
  x1 = d_item_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_itemdp';
  cond = {'RSrc','N'};
  x1 = d_item_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_itemdp';
  cond = {'RSrc','F'};
  x1 = d_item_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_itemdp';
  cond = {'F','N'};
  x1 = d_item_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % Source
  
  % side
  testType = 'side_sourcedp';
  cond = {'RSrc','N'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_sourcedp';
  cond = {'RSrc','F'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_sourcedp';
  cond = {'F','N'};
  x1 = d_source_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
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
  fprintf('Side:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_sourcedp';
  cond = {'RSrc','N'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_sourcedp';
  cond = {'RSrc','F'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_sourcedp';
  cond = {'F','N'};
  x1 = d_source_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
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
  fprintf('Question:     Item d'' with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
end


%% correlations IRK familiarity

% get the electrode numbers
elecGroupsNum = cell(size(elecGroups));
for i = 1:length(elecGroups)
  elecGroupsNum{i} = str2double(strrep(elecGroups{i},'e',''));
end

% accuracy
irk_fam_corr_side = abs([0.2154, 0.0566, 0.45, 0.1778, 0.2321, 0.14, 0.0833, 0.0964, 0.1705, 0.2845, 0.1233, 0.0132, 0.2811, 0.0526, 0.2143, 0.4202, 0.0169, 0.2091, 0.1358, 0.1346, 0.3281, 0.1299, 0.0877, 0.759, 0.1358, 0.0054, 0.0769, 0.1429, 0.1, 0.2545]);
irk_fam_inc_side = abs([0.0781, 0.024, 0.1014, 0.0565, 0.082, 0.1495, 0.0331, 0.0076, 0.0857, 0.1769, 0.082, 0.0043, 0.1859, 0.0037, 0.0703, 0.2837, 0.0041, 0.1077, 0.0417, 0.0709, 0.0778, 0.0292, 0.1043, 0.0498, 0.0488, 0.0083, 0.1091, 0.0917, 0.1075, 0.0696]);
irk_fam_corr_q = abs([0.0641, 0.0135, 0.2949, 0.1169, 0.1622, 0.1714, 0.0899, 0.033, 0.1, 0.2833, 0.1, 0.0089, 0.2673, 0.0132, 0.2667, 0.4981, 0.0492, 0.25, 0.0814, 0.0641, 0.1571, 0.1389, 0.0667, 0.4881, 0.1, 0.0047, 0.1429, 0.0833, 0.1127, 0.16]);
irk_fam_inc_q = abs([0.125, 0.0087, 0.125, 0.0973, 0.1239, 0.2577, 0.0556, 0.0234, 0.1053, 0.209, 0.0403, 0.0045, 0.1719, 0.0044, 0.0584, 0.3392, 0.041, 0.1343, 0.0777, 0.0776, 0.027, 0.0231, 0.11, 0.1138, 0.064, 0.0043, 0.1591, 0.0179, 0.0426, 0.1024]);

anaComponent = 'FN400';
%anaComponent = 'PON';

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
  roi = {'LPS','RPS'};
  %roi = {'LPS'};
  %roi = {'RPS'};
  chans = cat(2,elecGroupsNum{ismember(elecGroupsStr,roi)});
end

timeSamp = abs(prepost(1) * sampleRate) + [((timeMS(1) / (1000/sampleRate)) + 1) ((timeMS(2) / (1000/sampleRate)) + 1)];

% get the voltage for the correct region
volt_SN = [];
volt_SF = [];
volt_SRO = [];
volt_SRS = [];
volt_TN = [];
volt_TF = [];
volt_TRO = [];
volt_TRS = [];
for sub = 1:length(subjects)
  if ismember(subjects{sub},badSub)
    continue
  else
    volt_SN = [volt_SN mean(mean(mean(data_SN(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_SF = [volt_SF mean(mean(mean(data_SF(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_SRO = [volt_SRO mean(mean(mean(data_SRO(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_SRS = [volt_SRS mean(mean(mean(data_SRS(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TN = [volt_TN mean(mean(mean(data_TN(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TF = [volt_TF mean(mean(mean(data_TF(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TRO = [volt_TRO mean(mean(mean(data_TRO(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
    volt_TRS = [volt_TRS mean(mean(mean(data_TRS(sub).trial(:,chans,timeSamp(1):timeSamp(2)),3),2),1)];
  end
end

% correlate accuracy with voltage differences

cfg_plot = [];
% set up how the lines will look
cfg_plot.linespec = 'ko';
cfg_plot.linewidth = 1.5;
cfg_plot.marksize = 8;
cfg_plot.markcolor = 'w';

% RSrc - N
% RSrc - F
% F - N

if strcmp(anaComponent,'FN400')
  % FN400: item accuracy with RSrc - N
  fprintf('FN400\n');
  
  % Src Cor
  
  % side
  testType = 'side_irk_fam_corr';
  cond = {'RSrc','N'};
  x1 = irk_fam_corr_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_corr';
  cond = {'RSrc','F'};
  x1 = irk_fam_corr_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_corr';
  cond = {'F','N'};
  x1 = irk_fam_corr_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_corr';
  cond = {'RSrc','N'};
  x1 = irk_fam_corr_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_corr';
  cond = {'RSrc','F'};
  x1 = irk_fam_corr_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_corr';
  cond = {'F','N'};
  x1 = irk_fam_corr_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % Src Inc
  
  % side
  testType = 'side_irk_fam_inc';
  cond = {'RSrc','N'};
  x1 = irk_fam_inc_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_inc';
  cond = {'RSrc','F'};
  x1 = irk_fam_inc_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_inc';
  cond = {'F','N'};
  x1 = irk_fam_inc_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_inc';
  cond = {'RSrc','N'};
  x1 = irk_fam_inc_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_inc';
  cond = {'RSrc','F'};
  x1 = irk_fam_inc_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_inc';
  cond = {'F','N'};
  x1 = irk_fam_inc_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  

elseif strcmp(anaComponent,'PON')
  % P O/N: source accuracy with HSC - HSI
  fprintf('P O/N\n');
  
  % Src Cor
  
  % side
  testType = 'side_irk_fam_corr';
  cond = {'RSrc','N'};
  x1 = irk_fam_corr_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_corr';
  cond = {'RSrc','F'};
  x1 = irk_fam_corr_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_corr';
  cond = {'F','N'};
  x1 = irk_fam_corr_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_corr';
  cond = {'RSrc','N'};
  x1 = irk_fam_corr_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_corr';
  cond = {'RSrc','F'};
  x1 = irk_fam_corr_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_corr';
  cond = {'F','N'};
  x1 = irk_fam_corr_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SC: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SC');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SC with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % Src Inc
  
  % side
  testType = 'side_irk_fam_inc';
  cond = {'RSrc','N'};
  x1 = irk_fam_inc_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_inc';
  cond = {'RSrc','F'};
  x1 = irk_fam_inc_side(~ismember(subjects,badSub));
  y1 = (volt_SRS - volt_SF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % side
  testType = 'side_irk_fam_inc';
  cond = {'F','N'};
  x1 = irk_fam_inc_side(~ismember(subjects,badSub));
  y1 = (volt_SF - volt_SN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Side:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_inc';
  cond = {'RSrc','N'};
  x1 = irk_fam_inc_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_inc';
  cond = {'RSrc','F'};
  x1 = irk_fam_inc_q(~ismember(subjects,badSub));
  y1 = (volt_TRS - volt_TF);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
  % question
  testType = 'task_irk_fam_inc';
  cond = {'F','N'};
  x1 = irk_fam_inc_q(~ismember(subjects,badSub));
  y1 = (volt_TF - volt_TN);
  [rho p] = corr([x1' y1']);
  % plot it
  figure
  [m,b] = linefit(x1,y1,1,'r-');
  hold on
  plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
  hold off
  %title(sprintf('IRK Fam SI: %s - %s',cond{1},cond{2}));
  xlabel('IRK Fam SI');
  ylabel(sprintf('Voltage (%s): %s - %s','\muV',cond{1},cond{2}));
  %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
  axis([0 ceil(max(x1)) -4 4])
  axis square
  publishfig(gcf,0);
  figfilename = sprintf('corr_ga_%s_%s%s%d_%d.%s',testType,sprintf(repmat('%s_',1,length(cond)),cond{:}),sprintf(repmat('%s_',1,length(roi)),roi{:}),timeMS(1),timeMS(2),figFileExt);
  if saveFigs
    print(gcf,figPrintFormat,fullfile(saveDirFigs,figfilename));
  end
  fprintf('Question:     IRK Fam SI with %s - %s: Pearson''s r = %.4f (r^2 = %.4f), p = %.4f, m = %.3f, b = %.3f\n',cond{1},cond{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
  
end

%% get the trial counts for each subject to use in subsampling

% nTrial_N = num_N(~ismember(subjects,badSub));
% nTrial_F = num_F(~ismember(subjects,badSub));
% nTrial_RS = num_RS(~ismember(subjects,badSub));
% nTrial_SH = num_SH(~ismember(subjects,badSub));
% nTrial_TCR = num_TCR(~ismember(subjects,badSub));
% nTrial_THSC = num_THSC(~ismember(subjects,badSub));
% nTrial_THSI = num_THSI(~ismember(subjects,badSub));
% nTrial_TH = num_TH(~ismember(subjects,badSub));
% 
% %save('~/Desktop/trialCounts.mat','nTrial_*');
% 
% % nTrial_N = num_N(~ismember(subjects,badSub));
% % nTrial_HSC = num_HSC(~ismember(subjects,badSub));
% % nTrial_HSI = num_HSI(~ismember(subjects,badSub));
% % nTrial_H = num_H(~ismember(subjects,badSub));
% 
% % exp 1
% load('~/Desktop/trialCounts.mat');
% 
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% nExp3Sub = length(nTrial_N);
% nExp1Sub = sum(~ismember(subjects,badSub));
% subsRand = randperm(min([nExp3Sub,nExp1Sub]));
% 
% ind_N = nan(length(nTrial_N),1);
% ind_H = nan(length(nTrial_SH),1);
% ind_HSC = nan(length(nTrial_F),1);
% ind_HSI = nan(length(nTrial_RS),1);
% counter = 1;
% for sub = 1:length(subsRand)
%   ind_N(sub) = randperm(max([nTrial_N(subsRand(sub)),num_N(counter)]));
%   ind_N(sub) = ind_N(sub,1:min([nTrial_N(subsRand(sub)),num_N(counter)]));
%   
%   ind_H(sub) = randperm(max([nTrial_SH(subsRand(sub)),num_H(counter)]));
%   ind_H(sub) = ind_H(sub,1:min([nTrial_SH(subsRand(sub)),num_H(counter)]));
%   
%   ind_HSC(sub) = randperm(max([nTrial_F(subsRand(sub)),num_HSC(counter)]));
%   ind_HSC(sub) = ind_HSC(sub,1:min([nTrial_F(subsRand(sub)),num_HSC(counter)]));
%   
%   ind_HSI(sub) = randperm(max([nTrial_RS(subsRand(sub)),num_HSI(counter)]));
%   ind_HSI(sub) = ind_HSI(sub,1:min([nTrial_RS(subsRand(sub)),num_HSI(counter)]));
%   
%   counter = counter + 1;
% end
% 
% % exp 2
% load('~/Desktop/trialCounts.mat');
% 
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% nExp3Sub = length(nTrial_N);
% nExp2Sub = sum(~ismember(subjects,badSub));
% subs = randperm(min([nExp3Sub,nExp1Sub]));
% 
% ind_N = nan(length(nTrial_TCR),1);
% ind_H = nan(length(nTrial_TH),1);
% ind_HSC = nan(length(nTrial_THSC),1);
% ind_HSI = nan(length(nTrial_THSI),1);
% counter = 1;
% for sub = 1:length(subsRand)
%   ind_N(sub) = randperm(max([nTrial_TCR(subsRand(sub)),num_N(counter)]));
%   ind_N(sub) = ind_N(sub,1:min([nTrial_TCR(subsRand(sub)),num_N(counter)]));
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
%   cfg.trials = ind_N;
%   data_N_subsample(sub) = ft_preprocessing(cfg,data_N(sub));
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

% cfg = [];
% roi = {'LAS','RAS'};
% %roi = {'LAS'};
% %roi = {'RAS'};
% cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
% cfg.latency     = [0.3 0.5];

cfg = [];
roi = {'LPS','RPS'};
%roi = {'LPS'};
%roi = {'RPS'};
cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,roi)});
cfg.latency     = [0.5 0.8];

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

% both
timesel_N = find(ga_N.time >= cfg.latency(1) & ga_N.time <= cfg.latency(2));
timesel_F = find(ga_F.time >= cfg.latency(1) & ga_F.time <= cfg.latency(2));
timesel_RS = find(ga_RS.time >= cfg.latency(1) & ga_RS.time <= cfg.latency(2));
values_N = mean(mean(ga_N.individual(:,chan,timesel_N),2),3);
values_F = mean(mean(ga_F.individual(:,chan,timesel_F),2),3);
values_RS = mean(mean(ga_RS.individual(:,chan,timesel_RS),2),3);
sem_N = std(values_N)/sqrt(length(values_N));
sem_F = std(values_F)/sqrt(length(values_F));
sem_RS =  std(values_RS)/sqrt(length(values_RS));

fprintf('\t%s\t%s\t%s\t%s\n','N','F','RSrc');
% ga
[mean(values_N,1) mean(values_F,1) mean(values_RS,1)]
% sub avg
%[values_N values_F values_RS]

% both
stat_ttest_RSvsKNOW = ft_timelockstatistics(cfg,ga_RS,ga_F);
stat_ttest_FvsCR = ft_timelockstatistics(cfg,ga_F,ga_N);
stat_ttest_RSvsCR = ft_timelockstatistics(cfg,ga_RS,ga_N);

% side
fprintf('RS (M=%.3f; SEM=%.3f) vs KNOW (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(values_RS,1),sem_RS,mean(values_F,1),sem_F,stat_ttest_RSvsKNOW.df,stat_ttest_RSvsKNOW.stat,stat_ttest_RSvsKNOW.prob);
fprintf('KNOW (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(values_F,1),sem_F,mean(values_N,1),sem_N,stat_ttest_FvsCR.df,stat_ttest_FvsCR.stat,stat_ttest_FvsCR.prob);
fprintf('RS (M=%.3f; SEM=%.3f) vs CR  (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',mean(values_RS,1),sem_RS,mean(values_N,1),sem_N,stat_ttest_RSvsCR.df,stat_ttest_RSvsCR.stat,stat_ttest_RSvsCR.prob);
fprintf('\n');

% % plot to see the effect in each subject
% figure
% plot([values_F values_RS]','o-');
% xlim([0.5 2.5])
% title(sprintf('%.1fs--%.1fs, %s',cfg.latency(1),cfg.latency(2),sprintf(repmat('%s ',1,length(roi)),roi{:})));
% ylabel('Microvolts (\muV)');
% set(gca,'XTickLabel',{'','F','','RSrc',''})

% % matlab dependent samples ttest
% KNOWminRS = values_F - values_RS;
% [h,p,ci,stats] = ttest(KNOWminRS, 0, 0.05); % H0: mean = 0, alpha 0.05

%% mean amplitude line plots

cfg_plot = [];

if ismember('LAS',roi) || ismember('RAS',roi)
  conditions = {{'RSrc','F','N'}};
  cfg_plot.yminmax = [-5 -1];
elseif ismember('LPS',roi) || ismember('RPS',roi)
  conditions = {{'RSrc','F','N'}};
  cfg_plot.yminmax = [1 5];
end
%cfg_plot.yminmax = [floor(min([mean(values_HSC,1),mean(values_HSI,1),mean(values_N,1)])) ceil(max([mean(values_HSC,1),mean(values_HSI,1),mean(values_N,1)]))];

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
  
  testType = 'both_rkn';
  figure
  if length(cond) == 3
    % lines
    plot([mean(values_RS,1),mean(values_F,1),mean(values_N,1)],cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth);
    hold on
    % errorbars
    h = errorbar(1,mean(values_RS,1),sem_RS,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(2,mean(values_F,1),sem_F,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    
    h = errorbar(3,mean(values_N,1),sem_N,cfg_plot.side.linespec,'LineWidth',cfg_plot.errwidth);
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
    plot(1,mean(values_RS,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(2,mean(values_F,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
    plot(3,mean(values_N,1),cfg_plot.side.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.side.markcolor);
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
cfg.cond = {'CR','KNOW','RS'};
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
cfg.cond = {'CR','KNOW','RS'};
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
cfg.latency = [0.8 1.5];
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
stat_clus_RSvsKNOW = ft_timelockstatistics(cfg,ga_RS,ga_F);
stat_clus_RSvsCR = ft_timelockstatistics(cfg,ga_RS,ga_N);
stat_clus_FvsCR = ft_timelockstatistics(cfg,ga_F,ga_N);

%% plot the cluster statistics

cfg = [];
% p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
cfg.highlightsymbolseries = ['*','+','.','.','.'];
cfg.highlightcolorpos = [0.5 0 1];
cfg.highlightcolorneg = [0 0.5 0];
cfg.layout = ft_prepare_layout(cfg,ga_F);
cfg.contournum = 0;
cfg.emarker = '.';
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.zlim = [-5 5];
%cfg.interplimits = 'electrodes';
ft_clusterplot(cfg,stat_clus_RSvsKNOW); % none?
ft_clusterplot(cfg,stat_clus_RSvsCR); % pos front, neg back?
ft_clusterplot(cfg,stat_clus_FvsCR); % weak pos front, neg back?

%% older stuff

% create contrast
ga_FvsRS = ga_F;
ga_FvsRS.avg = ga_F.avg - ga_RS.avg;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
pos = stat_clus_FvsRS.posclusterslabelmat==1;
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
  ft_topoplotER(cfg,ga_FvsRS);
end

% create contrast
ga_HSCvsCR = ga_HSC;
ga_HSCvsCR.avg = ga_HSC.avg - ga_N.avg;
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
ga_HSIvsCR.avg = ga_HSI.avg - ga_N.avg;
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
