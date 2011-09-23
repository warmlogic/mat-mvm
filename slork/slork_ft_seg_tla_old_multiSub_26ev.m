% Do a nonparametric statistics clustering analysis for timelocked EEG
% (ERPs)

% See Maris and Oostenveld (2007) for a description of statistics

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
    'SLORK033'};
% %    'SLORK001'... % bad subject
% %    'SLORK021'
% %    'SLORK028'

%% set up parameters
homeDir = getenv('HOME');
%serverDir = '/Volumes/curranlab/Data/SLORK/eeg/eppp_nobackup/-1.0_2.0_26ev';
%serverLocalDir = '/Volumes/RAID/curranlab/Data/SLORK/eeg/eppp_nobackup/-1.0_2.0_26ev';

serverDir = '/Volumes/curranlab/Data/SLORK/eeg/nspp_nobackup/-1.0_2.0_26ev';
serverLocalDir = '/Volumes/RAID/curranlab/Data/SLORK/eeg/nspp_nobackup/-1.0_2.0_26ev';

if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  dataroot = fullfile(homeDir,'data/SLORK/eeg');
end

saveDir = fullfile(dataroot,'ft_data');
if ~exist(saveDir,'dir')
  mkdir(saveDir);
end

elecfile = fullfile(homeDir,'eeg/GSN_HydroCel_129_short.sfp');
locsFormat = 'besa_sfp';
elec = read_sens(elecfile,'fileformat',locsFormat);

% pre- and post-stimulus times to save, in seconds
prepost = [-1 2];

sampleRate = 250;

%eventValues = sort({'SCR','SHSC','SHSI','TCR','THSC','THSI'});

%eventValues = sort({'SCRM','SCRS','SFAF','SFAO','SFAS','SMM','SMS','SSCF','SSCO','SSCS','SSIF','SSIO','SSIS','TCRM','TCRS','TFAF','TFAO','TFAS','TMM','TMS','TSCF','TSCO','TSCS','TSIF','TSIO','TSIS'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in FieldTrip structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sub = 1:length(subjects)
  fprintf('Working on %s...\n',subjects{sub});
  
  filename = sprintf('%s_data_preproc_%0.1f_%.1f.mat',subjects{sub},prepost(1),prepost(2));
  dataSaveFile = fullfile(saveDir,filename);
  
  if ~exist(dataSaveFile,'file')
    data_preproc = seg2ft_nsmat(dataroot,subjects{sub},prepost,sampleRate,elecfile);
    
    fprintf('Saving %s...\n',subjects{sub});
    save(dataSaveFile,'data_preproc','prepost','elec');
  end
  
  eventValuesFile = fullfile(saveDir,'eventValues.mat');
  if ~exist(eventValuesFile,'file')
    if ~exist('eventValues','var')
      eventValues = {};
      for i = 1:length(data_preproc)
        eventValues = [eventValues data_preproc(i).cfg.trialdef.eventvalue];
      end
    end
    save(eventValuesFile,'eventValues','prepost','elec');
  end
  
  fprintf('Done.\n');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ft_timelockanalysis; make sure events go in alphabetical order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear data_SCRM data_SCRS data_SFAF data_SFAO data_SFAS data_SMM data_SMS data_SSCF data_SSCO data_SSCS data_SSIF data_SSIO data_SSIS data_TCRM data_TCRS data_TFAF data_TFAO data_TFAS data_TMM data_TMS data_TSCF data_TSCO data_TSCS data_TSIF data_TSIO data_TSIS

%clear data_SCR data_SHSC data_SHSI data_TCR data_THSC data_THSI

load(fullfile(saveDir,'eventValues.mat'));

clear data_* data

cfg_pp = [];
% baseline correct
cfg_pp.blc = 'yes';
cfg_pp.blcwindow  = [-0.2 0];
% apply lowpass filter at 35 Hz
cfg_pp.lpfilter = 'yes';
cfg_pp.lpfreq = 35;

cfg_proc = [];
cfg_proc.keeptrials = 'yes';

% cfg_tlb = [];
% cfg_tlb.baseline = [-0.2 0];
% cfg_tlb.channel = 'all';

% for keeping track of the number of trials that each subject has
for i = 1:length(eventValues)
  eventVal = eventValues{i};
  eval(sprintf('num_%s = zeros(%d,1);',eventVal,length(subjects)));
end

%matlabpool open

%parfor sub = 1:length(subjects)
for sub = 1:length(subjects)
  filename = sprintf('%s_data_preproc_%0.1f_%.1f.mat',subjects{sub},prepost(1),prepost(2));
  dataSaveFile = fullfile(saveDir,filename);
  fprintf('Loading %s...\n',subjects{sub});
  data = load(dataSaveFile);
  
  for ev = 1:length(data.data_preproc)
    eventVal_sub =  data.data_preproc(ev).cfg.trialdef.eventvalue;
    
    for i = 1:length(eventValues)
      eventVal = eventValues{i};
      
      if strcmp(eventVal_sub,eventVal)
        if data.data_preproc(ev).hdr.nTrials > 0
          
          tpp_str = sprintf('data_%s_pp = ft_preprocessing(cfg_pp,data.data_preproc(ev));',eventVal);
          eval(tpp_str);
          
          tla_str = sprintf('data_%s(sub) = ft_timelockanalysis(cfg_proc,data_%s_pp);',eventVal,eventVal);
          eval(tla_str);
          
          %%%%%%%%
          
          %tla_str = sprintf('data_%s_raw = ft_timelockanalysis(cfg_proc,data.data_preproc(ev));',eventVal);
          %eval(tla_str);
          
          %tlb_str = sprintf('data_%s(sub) = ft_timelockbaseline(cfg_tlb,data_%s_raw);',eventVal,eventVal);
          %eval(tlb_str);
          
          % get the var back
          %tla_str = sprintf('data_%s_var(sub) = ft_timelockanalysis(cfg_proc,data_%s(sub));',eventVal,eventVal);
          %eval(tla_str);
          
          %eval(sprintf('clear data_%s_raw',eventVal));
          clear data_*_pp
        else
          eval(sprintf('data_%s(sub).avg = [];',eventVal));
          eval(sprintf('data_%s(sub).var = [];',eventVal));
          eval(sprintf('data_%s(sub).fsample = [];',eventVal));
          eval(sprintf('data_%s(sub).time = [];',eventVal));
          eval(sprintf('data_%s(sub).dof = [];',eventVal));
          eval(sprintf('data_%s(sub).label = {};',eventVal));
          eval(sprintf('data_%s(sub).trial = [];',eventVal));
          eval(sprintf('data_%s(sub).dimord = [];',eventVal));
          eval(sprintf('data_%s(sub).elec = [];',eventVal));
          eval(sprintf('data_%s(sub).cfg = [];',eventVal));
          
        end
        
        % keep track of the number of trials for each subject
        eval(sprintf('num_%s(sub) = size(data_%s(sub).trial,1);',eventVal,eventVal));
        
        break
      end
    end
  end
  %% by hand
  %   data_SCRM(sub) = ft_timelockanalysis(cfg,data.data_preproc(1));
  %   data_SCRS(sub) = ft_timelockanalysis(cfg,data.data_preproc(2));
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
  
  %% ending
  fprintf('Done.\n');
end

% for i = 1:length(eventValues)
%   eventVal = eventValues{i};
%   eval(sprintf('num_%s(sub) = size(data_%s(sub).trial,1);',eventVal,eventVal));
% end

%matlabpool close

% % ft_appenddata must be done after preprocessing, before ft_timelockanalysis
% cfg = [];
% for sub = 1:length(subjects)
%   data_SHSC_pp(sub) = ft_appenddata(cfg,data_SSCS_pp(sub),data_SSCO_pp(sub),data_SSCF_pp(sub));
%   data_SHSI_pp(sub) = ft_appenddata(cfg,data_SSIS_pp(sub),data_SSIO_pp(sub),data_SSIF_pp(sub));
%   data_SCR_pp(sub) = ft_appenddata(cfg,data_SCRS_pp(sub),data_SCRM_pp(sub));
%   data_THSC_pp(sub) = ft_appenddata(cfg,data_TSCS_pp(sub),data_TSCO_pp(sub),data_TSCF_pp(sub));
%   data_THSI_pp(sub) = ft_appenddata(cfg,data_TSIS_pp(sub),data_TSIO_pp(sub),data_TSIF_pp(sub));
%   data_TCR_pp(sub) = ft_appenddata(cfg,data_TCRS_pp(sub),data_TCRM_pp(sub));
% end

% save so we don't have to preprocess later
for i = 1:length(eventValues)
  eventVal = eventValues{i};
  tlaFilename = sprintf('data_%s_tla_%0.1f_%.1f.mat',eventVal,prepost(1),prepost(2));
  tlaSaveFile = fullfile(saveDir,tlaFilename);
  eval(sprintf('save(tlaSaveFile,''data_%s'',''num_%s'');',eventVal,eventVal));
end

% savedFiles = dir(fullfile(saveDir,'data_*_tla_*.mat'));
% for i = 1:length(savedFiles)
%   fprintf('Loading %s...',savedFiles(i).name);
%   load(fullfile(saveDir,savedFiles(i).name));
%   fprintf('Done.\n');
% end

% load full conditions (i.e., 6 events)
savedFiles = dir(fullfile(strrep(saveDir,'_26ev',''),'data_*_tla_*.mat'));
for i = 1:length(savedFiles)
  fprintf('Loading %s...',savedFiles(i).name);
  load(fullfile(strrep(saveDir,'_26ev',''),savedFiles(i).name));
  fprintf('Done.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove the degrees of freedom? field. I'm not sure what it does, but the tutorial says to.
for i = 1:length(eventValues)
  eventVal = eventValues{i};
  eval(sprintf('data_%s = rmfield(data_%s,''dof'');',eventVal,eventVal));
end

% events actually of interest: S/TSCS vs S/TSI* (S, O, F) vs S/TCR* (S, M)
%[num_SCRS num_SCRM num_SSCS num_SSCO num_SSCF num_SSIS num_SSIO num_SSIF]
%[num_TCRS num_TCRM num_TSCS num_TSCO num_TSCF num_TSIS num_TSIO num_TSIF]

% compare to S/TCR S/THSC S/THSI
%[num_SCR num_SHSC num_SHSI]
%[num_TCR num_THSC num_THSI]

lowNum = [];
thresh = 20;
% % get the number of trials
% for sub = 1:length(subjects)
%   for i = 1:length(eventValues)
%     eventVal = eventValues{i};
%     if eval(sprintf('num_%s(%d) < %d',eventVal,sub,thresh))
%       lowNum = [lowNum sub];
%     end
%   end
% end
% 
% % who we're rejecting
% subjects(unique(lowNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d' < 1.0
%badBehSub = {'SLORK010','SLORK011','SLORK025','SLORK030','SLORK031','SLORK032'};

% 029's ERPs go REALLY negative
%badBehSub = {'SLORK013','SLORK029'};
badBehSub = {};

badSub = unique([badBehSub subjects(lowNum)]);

%% set up analysis strings
for i = 1:length(eventValues)
  eventVal = eventValues{i};
  eval(sprintf('str_%s = [];',eventVal));
end

for sub = 1:length(subjects)
  if ismember(subjects{sub},badSub)
    continue
  else
    if sub == 1
      for i = 1:length(eventValues)
        eventVal = eventValues{i};
        this_str = sprintf('data_%s(%d)',eventVal,sub);
        if ~isempty(eval(sprintf('%s.trial',this_str)))
          eval(sprintf('str_%s = ''%s'';',eventVal,this_str));
        end
      end
      % str_SSCS = 'data_SSCS(1),data_SSCS(2)';
    elseif sub > 1
      for i = 1:length(eventValues)
        eventVal = eventValues{i};
        this_str = sprintf('data_%s(%d)',eventVal,sub);
        if ~isempty(eval(sprintf('%s.trial',this_str)))
          eval(sprintf('str_%s = [%s,'',%s''];',eventVal,sprintf('str_%s',eventVal),this_str));
        end
      end
    end
  end
end
for i = 1:length(eventValues)
  eventVal = eventValues{i};
  if eval(sprintf('strcmp(str_%s(1),'','')',eventVal));
    eval(sprintf('str_%s = str_%s(2:end)',eventVal,eventVal));
  end
end

%% get the grand average
cfg = [];
cfg.keepindividual = 'yes';

for i = 1:length(eventValues)
  eventVal = eventValues{i};
  if ~isempty(eval(sprintf('str_%s',eventVal)))
    if isempty(strfind(eval(sprintf('str_%s',eventVal)),','))
      fprintf('Only 1 subject in %s grand average!\n',eventVal);
    end
    %eval(sprintf('data_%s = rmfield(data_%s,''elec'');',eventVal,eventVal));
    eval(sprintf('ga_%s = ft_timelockgrandaverage(cfg,%s);',eventVal,eval(sprintf('str_%s',eventVal))));
  end
end

%% add the electrode locations again
for i = 1:length(eventValues)
  eventVal = eventValues{i};
  eval(sprintf('ga_%s.elec = elec;',eventVal));
end

%% set up channel groups
elecGroups = {...
    {'e32','e33','e38','e39','e43','e44','e128'},... % LAI
    {'e1','e114','e115','e120','e121','e122','e125'},... % RAI
    {'e12','e13','e19','e20','e24','e28','e29'},... % LAS
    {'e4','e5','e111','e112','e117','e118','e124'},... % RAS
    {'e37','e42','e52','e53','e54','e60','e61'},... % LPS
    {'e78','e79','e85','e86','e87','e92','e93'},... % RPS
    {'e57','e58','e63','e64','e65','e68','e69'},... %LPI
    {'e89','e90','e94','e95','e96','e99','e100'} % RPI
             };
elecGroupsStr = {'LAI','RAI','LAS','RAS','LPS','RPS','LPI','RPI'};

%% plot the conditions
cfg = [];
cfg.showlabels = 'yes'; 
cfg.fontsize = 9; 
%cfg.ylim = [-3e-13 3e-13];
cfg.interactive = 'yes';
cfg.layout = ft_prepare_layout(cfg,ga_SCR);

figure
ft_multiplotER(cfg,ga_SCR,ga_SHSC,ga_SHSI);
figure
ft_multiplotER(cfg,ga_TCR,ga_THSC,ga_THSI);


%cfg.interplimits = 'electrodes';
%ft_topoplotER(cfg,ga_SCR);

%% plot the conditions

% events actually of interest: S/TSCS vs S/TSI* (S, O, F) vs S/TCR* (S, M);
% compare to S/TCR S/THSC S/THSI

cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 9;
cfg.linewidth = 1;
cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'RAS'})});
figure
ft_singleplotER(cfg,ga_SCR,ga_SHSC,ga_SHSI,ga_SSCS);
legend('SCR','SHSC','SHSI','SSCS','Location','NorthWest');

%[num_SCRS num_SCRM num_SSCS num_SSCO num_SSCF num_SSIS num_SSIO num_SSIF]

figure
ft_singleplotER(cfg,ga_SCRS,ga_SSCS,ga_SSIS);
legend('SCRS','SSCS','SSIS','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_SCRM,ga_SSCO,ga_SSIO);
legend('SCRM','SSCO','SSIO','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_SCRM,ga_SSCF,ga_SSIF);
legend('SCRM','SSCF','SSIF','Location','NorthWest');

figure
ft_singleplotER(cfg,ga_SCRS,ga_SCRM);
legend('SCRS','SCRM','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_SSCS,ga_SSCO,ga_SSCF);
legend('SSCS','SSCO','SSCF','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_SSIS,ga_SSIO,ga_SSIF);
legend('SSIS','SSIO','SSIF','Location','NorthWest');

%[num_TCRS num_TCRM num_TSCS num_TSCO num_TSCF num_TSIS num_TSIO num_TSIF]

figure
ft_singleplotER(cfg,ga_TCRS,ga_TSCS,ga_TSIS);
legend('TCRS','TSCS','TSIS','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_TCRM,ga_TSCO,ga_TSIO);
legend('TCRM','TSCO','TSIO','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_TCRM,ga_TSCF,ga_TSIF);
legend('TCRM','TSCF','TSIF','Location','NorthWest');

figure
ft_singleplotER(cfg,ga_TCRS,ga_TCRM);
legend('TCRS','TCRM','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_TSCS,ga_TSCO,ga_TSCF);
legend('TSCS','TSCO','TSCF','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_TSIS,ga_TSIO,ga_TSIF);
legend('TSIS','TSIO','TSIF','Location','NorthWest');




% create contrast
ga_SHSCvsSHSI = ga_SHSC;
%ga_SHSCvsSHSI = rmfield(ga_SHSCvsSHSI,'individual');
ga_SHSCvsSHSI.avg = ga_SHSC.avg - ga_SHSI.avg;
%ga_SHSCvsSHSI.individual = ga_SHSC.individual - ga_SHSI.individual;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
%pos = stat_clus_SHSCvsSHSI.posclusterslabelmat==1;
for k = 1:length(j)-1;
  subplot(4,5,k);
  cfg = [];
  cfg.interactive = 'yes';
  cfg.zlim = [-1 1];
  cfg.xlim=[j(k) j(k+1)];
  %pos_int = mean(pos(:,m(k):m(k+1))')'; % mean(pos(:,m(k):m(k+1)),2);
  cfg.highlight = 'on';
  %cfg.highlightchannel = find(pos_int==1);
  cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,{'RAS'})});
  %cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  ft_topoplotER(cfg,ga_SHSCvsSHSI);
end

% create contrast
ga_SHSCvsSCR = ga_SHSC;
%ga_SHSCvsSHSI = rmfield(ga_SHSCvsSHSI,'individual');
ga_SHSCvsSCR.avg = ga_SHSC.avg - ga_SCR.avg;
%ga_SHSCvsSHSI.individual = ga_SHSC.individual - ga_SHSI.individual;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
%pos = stat_clus_SHSCvsSHSI.posclusterslabelmat==1;
for k = 1:length(j)-1;
  subplot(4,5,k);
  cfg = [];
  cfg.interactive = 'yes';
  cfg.zlim = [-1 1];
  cfg.xlim=[j(k) j(k+1)];
  %pos_int = mean(pos(:,m(k):m(k+1))')'; % mean(pos(:,m(k):m(k+1)),2);
  cfg.highlight = 'on';
  %cfg.highlightchannel = find(pos_int==1);
  cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,{'RAS'})});
  %cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  ft_topoplotER(cfg,ga_SHSCvsSCR);
end

% create contrast
ga_SHSIvsSCR = ga_SHSI;
%ga_SHSCvsSHSI = rmfield(ga_SHSCvsSHSI,'individual');
ga_SHSIvsSCR.avg = ga_SHSI.avg - ga_SCR.avg;
%ga_SHSCvsSHSI.individual = ga_SHSC.individual - ga_SHSI.individual;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
%pos = stat_clus_SHSCvsSHSI.posclusterslabelmat==1;
for k = 1:length(j)-1;
  subplot(4,5,k);
  cfg = [];
  cfg.interactive = 'yes';
  cfg.zlim = [-1 1];
  cfg.xlim=[j(k) j(k+1)];
  %pos_int = mean(pos(:,m(k):m(k+1))')'; % mean(pos(:,m(k):m(k+1)),2);
  cfg.highlight = 'on';
  %cfg.highlightchannel = find(pos_int==1);
  cfg.highlightchannel = cat(2,elecGroups{ismember(elecGroupsStr,{'RAS'})});
  %cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  ft_topoplotER(cfg,ga_SHSIvsSCR);
end

cfg.channel = 'e75';
figure
ft_singleplotER(cfg,ga_SCR,ga_SSCS,ga_SHSI);
legend('SCR','SHSC','SHSI','Location','NorthWest');

figure
ft_singleplotER(cfg,ga_SCR,ga_SSCS,ga_SHSI);
legend('SCR','SHSC','SHSI','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_TCR,ga_TSCS,ga_THSI);
legend('TCR','THSC','THSI','Location','NorthWest');

cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'LAS'})});
figure
ft_singleplotER(cfg,ga_SCR,ga_SSCS,ga_SHSI);
legend('SCR','SHSC','SHSI','Location','NorthWest');
figure
ft_singleplotER(cfg,ga_TCR,ga_TSCS,ga_THSI);
legend('TCR','THSC','THSI','Location','NorthWest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% descriptive statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% for iSub = 1:length(subjects)
%   subplot(1,3,iSub);
%   plot(ga_HSC.time,squeeze(ga_HSC.individual(iSub,52,:)));
%   hold on
%   plot(ga_HSI.time,squeeze(ga_HSI.individual(iSub,52,:)),'r');
%   title(strcat('subject ',num2str(iSub)));
%   %ylim([0 1.9e-13]);
%   xlim([-0.5 1.5]);
% end
% subplot(1,3,length(subjects)+1); 
% text(0.5,0.5,'HSC','color','b') ;text(0.5,0.3,'HSI','color','r')
% axis off

% chan = 52;
% time = [0.5 0.8];
% timesel_HSC = find(ga_HSC.time >= time(1) & ga_HSC.time <= time(2));
% values_HSC = mean(ga_HSC.individual(:,chan,timesel_HSC),3);
% timesel_HSI = find(ga_HSI.time >= time(1) & ga_HSI.time <= time(2));
% values_HSI = mean(ga_HSI.individual(:,chan,timesel_HSC),3);
% % plot to see the effect in each subject
% M = [values_HSC,values_HSI];
% figure; plot(M','o-'); xlim([0.5 2.5])
% 
% % dependent samples ttest
% HSCminHSI = values_HSC - values_HSI;
% [h,p,ci,stats] = ttest(HSCminHSI, 0, 0.05); % H0: mean = 0, alpha 0.05

% fieldtrip ttest
cfg = [];
cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'RAS'})});
cfg.latency     = [0.3 0.5];
%cfg.channel = cat(2,elecGroups{ismember(elecGroupsStr,{'LPS'})});
%cfg.latency     = [0.5 0.8];
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
cfg.ivar               = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar               = 2; % the 2nd row in cfg.design contains the subject number

% side
stat_ttest_SHSCvsSHSI = timelockstatistics(cfg,ga_SHSC,ga_SHSI);
stat_ttest_SHSCvsSCR = timelockstatistics(cfg,ga_SHSC,ga_SCR);
stat_ttest_SHSIvsSCR = timelockstatistics(cfg,ga_SHSI,ga_SCR);

% task
stat_ttest_THSCvsTHSI = timelockstatistics(cfg,ga_THSC,ga_THSI);
stat_ttest_THSCvsTCR = timelockstatistics(cfg,ga_THSC,ga_TCR);
stat_ttest_THSIvsTCR = timelockstatistics(cfg,ga_THSI,ga_TCR);

% side vs task
stat_ttest_SHSCvsTHSC = timelockstatistics(cfg,ga_SHSC,ga_THSC);
stat_ttest_SHSIvsTHSI = timelockstatistics(cfg,ga_SHSI,ga_THSI);
stat_ttest_SCRvsTCR = timelockstatistics(cfg,ga_SCR,ga_TCR);

chan = str2double(strrep(cfg.channel,'e',''));

% side
timesel_SHSC = find(ga_SHSC.time >= cfg.latency(1) & ga_SHSC.time <= cfg.latency(2));
values_SHSC = mean(ga_SHSC.individual(:,chan,timesel_SHSC),3);
timesel_SHSI = find(ga_SHSI.time >= cfg.latency(1) & ga_SHSI.time <= cfg.latency(2));
values_SHSI = mean(ga_SHSI.individual(:,chan,timesel_SHSC),3);
timesel_SCR = find(ga_SCR.time >= cfg.latency(1) & ga_SCR.time <= cfg.latency(2));
values_SCR = mean(ga_SCR.individual(:,chan,timesel_SCR),3);

% task
timesel_THSC = find(ga_THSC.time >= cfg.latency(1) & ga_THSC.time <= cfg.latency(2));
values_THSC = mean(ga_THSC.individual(:,chan,timesel_THSC),3);
timesel_THSI = find(ga_THSI.time >= cfg.latency(1) & ga_THSI.time <= cfg.latency(2));
values_THSI = mean(ga_THSI.individual(:,chan,timesel_THSC),3);
timesel_TCR = find(ga_TCR.time >= cfg.latency(1) & ga_TCR.time <= cfg.latency(2));
values_TCR = mean(ga_TCR.individual(:,chan,timesel_TCR),3);

[mean(mean(values_SCR)) mean(mean(values_SHSC)) mean(mean(values_SHSI))]
%[mean(values_SCR,2) mean(values_SHSC,2) mean(values_SHSI,2)]

[mean(mean(values_TCR)) mean(mean(values_THSC)) mean(mean(values_THSI))]
%[mean(values_TCR,2) mean(values_THSC,2) mean(values_THSI,2)]

for i = 1:length(eventValues)
  eventVal = eventValues{i};
  eval(sprintf('timesel_%s = find(ga_%s.time >= %f & ga_%s.time <= %f);',eventVal,eventVal,cfg.latency(1),eventVal,cfg.latency(2)));
  eval(sprintf('values_%s = mean(ga_%s.individual(:,chan,timesel_%s),3);',eventVal,eventVal,eventVal));
end

[mean(mean(values_SCR)) mean(mean(values_SHSC)) mean(mean(values_SHSI)) mean(mean(values_SSCS))]
%[mean(values_SCR,2) mean(values_SHSC,2) mean(values_SHSI,2) mean(values_SSCS,2)]

[mean(mean(values_TCR)) mean(mean(values_THSC)) mean(mean(values_THSI)) mean(mean(values_TSCS))]
%[mean(values_TCR,2) mean(values_THSC,2) mean(values_THSI,2) mean(values_TSCS,2)]
%[num_SCRS num_SCRM num_SSCS num_SSCO num_SSCF num_SSIS num_SSIO num_SSIF]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cluster statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
stat_clus_SHSCvsSHSI = timelockstatistics(cfg,ga_SHSC,ga_SHSI);
stat_clus_SHSCvsSCR = timelockstatistics(cfg,ga_SHSC,ga_SCR);
stat_clus_SHSIvsSCR = timelockstatistics(cfg,ga_SHSI,ga_SCR);

% task
stat_clus_THSCvsTHSI = timelockstatistics(cfg,ga_THSC,ga_THSI);
stat_clus_THSCvsTCR = timelockstatistics(cfg,ga_THSC,ga_TCR);
stat_clus_THSIvsTCR = timelockstatistics(cfg,ga_THSI,ga_TCR);

% side vs task
stat_clus_SHSCvsTHSC = timelockstatistics(cfg,ga_SHSC,ga_THSC);
stat_clus_SHSIvsTHSI = timelockstatistics(cfg,ga_SHSI,ga_THSI);
stat_clus_SCRvsTCR = timelockstatistics(cfg,ga_SCR,ga_TCR);

%% plot the cluster statistics

cfg = [];
% p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
cfg.highlightsymbolseries = ['*','+','.','.','.'];
cfg.highlightcolorpos = [0.5 0 1];
cfg.highlightcolorneg = [0 0.5 0];
cfg.layout = ft_prepare_layout(cfg,ga_SHSC);
cfg.contournum = 0;
cfg.emarker = '.';
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.zlim = [-5 5];
%cfg.interplimits = 'electrodes';
ft_clusterplot(cfg,stat_clus_SHSCvsSHSI); % late FN400? P O/N
ft_clusterplot(cfg,stat_clus_SHSCvsSCR); % P O/N, late pos front
ft_clusterplot(cfg,stat_clus_SHSIvsSCR); % late P O/N dipole? LPS not sig

ft_clusterplot(cfg,stat_clus_THSCvsTHSI); % P O/N
ft_clusterplot(cfg,stat_clus_THSCvsTCR); % P O/N, late pos front
ft_clusterplot(cfg,stat_clus_THSIvsTCR); % late pos front dipole?

ft_clusterplot(cfg,stat_clus_SHSCvsTHSC); % none
ft_clusterplot(cfg,stat_clus_SHSIvsTHSI); % none
ft_clusterplot(cfg,stat_clus_SCRvsTCR); % none

%% older stuff

% create contrast
ga_SHSCvsSHSI = ga_SHSC;
%ga_SHSCvsSHSI = rmfield(ga_SHSCvsSHSI,'individual');
ga_SHSCvsSHSI.avg = ga_SHSC.avg - ga_SHSI.avg;
% make a plot
figure
j = [0:0.05:1.0];
m = round(linspace(1,250,length(j)));
%pos = stat_clus_SHSCvsSHSI.posclusterslabelmat==1;
for k = 1:length(j)-1;
  subplot(4,5,k);
  cfg = [];
  cfg.xlim=[j(k) j(k+1)];
  cfg.zlim = [-5 5];
  %pos_int = mean(pos(:,m(k):m(k+1))')'; % mean(pos(:,m(k):m(k+1)),2);
  cfg.highlight = 'on';
  %cfg.highlightchannel = find(pos_int==1);
  cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  ft_topoplotER(cfg,ga_SHSCvsSHSI.avg);
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
