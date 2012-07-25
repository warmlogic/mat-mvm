%dataroot = fullfile(getenv('HOME'),'Downloads','FRCE_data_classification');

dataroot = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data';
dataroot = fullfile(dataroot,'Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto');

% % dpss 4:128 log spaced
% dataroot = fullfile(dataroot,'pow_mtmconvol_dpss_pow_-500_1500_4_128');
% hanning 4:1:64
% dataroot = fullfile(dataroot,'pow_mtmconvol_hanning_pow_-500_1500_4_64');
% % erp
%dataroot = fullfile(dataroot,'tla_-1250_2250');

% wavelet, 4:1:100
dataroot = fullfile(dataroot,'pow_wavelet_w5_pow_-625_1615_4_100');

% load in the dataset
%adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_40/analysisDetails.mat';
%adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/tla_-1250_2250/analysisDetails.mat';
adFile = fullfile(dataroot,'analysisDetails.mat');
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
%ana.eventValues = {exper.eventValues};
%ana.eventValues = {{'AudForg','AudReca'},{'VisForg','VisReca'},{'Forg','Reca'},{'Aud','Vis'}};
%ana.eventValues = {{'AudForg','AudReca'}};
% ana.eventValues = {{'VisForg','VisReca'}};
%ana.eventValues = {{'VisForg','VisReca'},{'AudForg','AudReca'}};
%ana.eventValues = {{'Forg','Reca'}};
ana.eventValues = {{'Aud','Vis'}};

% pre-defined ROIs in this function
ana = mm_ft_elecGroups(ana);

% redefine which subjects to load
exper.subjects = {'FRCE 05'};

[data_freq] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

%[data_raw] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'raw');

%equatetrials = true;
equatetrials = false;

addtrials = false; %double number of trials in condition that starts out most frequent

subtract_trials = false; % subtract auditory trials to simulate skewed classes

%% beh: accuracy

ses = 1;

for typ = 1:length(ana.eventValues)
  % find the lowest number of trials within each subject
  evNum = nan(length(ana.eventValues{typ}),length(exper.subjects));
  for sub = 1:length(exper.subjects)
    for cnd = 1:length(ana.eventValues{typ});
    evNum(cnd,sub) = size(data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.powspctrm,1);
    end
  end
  fprintf('trial counts\n');
  disp(evNum);
  fprintf('trial count average\n');
  disp(mean(evNum,2));
  fprintf('accuracy\n');
  disp(evNum(2,:) ./ sum(evNum,1));
  fprintf('accuracy average\n');
  disp(mean(evNum(2,:) ./ sum(evNum,1)));
end

%% ERP: process the data

%subInd = find(ismember(exper.subjects,cfg_prep.subjs));
ses = 1;

if equatetrials
  % initialize random number generator
  rng('shuffle','twister');
  
  % find the lowest number of trials within each subject
  lowEvNum = nan(length(exper.subjects),1);
  for sub = 1:length(exper.subjects)
    lowEvNum(sub) = Inf;
    for i = 1:length(ana.eventValues{typ})
      if size(data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trial,2) < lowEvNum(sub)
        lowEvNum(sub) = size(data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trial,2);
      end
    end
  end
end

% initialize the struct
data_erp = struct;

cfg_ana = [];
%cfg_ana.chanStr = {'right'};
%cfg_ana.chanStr = {'LPS'};
%cfg_ana.chanStr = {'center74'};
cfg_ana.chanStr = {'center91'};
%cfg_ana.chanStr = {'all129'};

% set the analysis parameters
cfg = [];
cfg.parameter = 'trial';
cfg.keeptrials = 'yes'; % classifiers operate on individual trials
cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)}));

% get the ERP data
for sub = 1:length(exper.subjects)
  for i = 1:length(ana.eventValues{typ})
    if equatetrials
      % choose a random subset to equate trial numbers
      evNums = randperm(size(data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trial,2));
      cfg.trials = sort(evNums(1:lowEvNum(sub)));
    end
    data_erp.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data = ft_timelockanalysis(cfg,data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data);
  end
end

%% ERP: classification

ses = 1;
cfg_ana = [];
%cfg_ana.latency = [-0.5 0];
%cfg_ana.latency = [0.5 1.0];
%cfg_ana.latency = [0.3 1.0];
%cfg_ana.latency = [.15 0.4];
cfg_ana.latency = [0 1.0];

%cfg_ana.chanStr = 'right';
%cfg_ana.chanStr = {'center74'};
cfg_ana.chanStr = {'center91'};

accuracy = nan(length(exper.subjects),1);
binomial = nan(length(exper.subjects),1);
continTab = cell(length(exper.subjects),1);

cfg = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
%cfg.parameter = 'trial';
cfg.method = 'crossvalidate';
cfg.nfolds = 5;
cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)}));


cfg.latency = cfg_ana.latency;
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg_ana.nChan = length(cfg.channel);

% z-transform, feature selection (map data to a different space), SVM
%cfg_ana.method = 'cspFs10SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_filterer('maxfeatures',10) ft_mv_svm};

%cfg_ana.method = 'cspSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_svm};
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan,'numpatterns',1,'outputdatatype','logpowcsp','filttype','CSP0') ft_mv_svm};

%cfg_ana.method = 'fs10SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_filterer('maxfeatures',10) ft_mv_svm};

% z-transform, feature selection in the original space
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0.1)};

% basic SVM
%cfg_ana.method = 'SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_svm};

% unregularized linear regression
%cfg_ana.method = 'unregLinR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0,'family','gaussian')};
% Didn't work

% unregularized logistic regression
%cfg_ana.method = 'unregLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0,'family','binomial')};
% crashes

% L1 regularized logistic regression
%cfg_ana.method = 'L1regLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',1,'lambda',1,'family','binomial')};
% doesn't perform very well

% L2 regularized logistic regression
cfg_ana.method = 'L2regLogR';
cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0,'lambda',1,'family','binomial')};
% seems to have the best performance

% elastic net logistic regression
%cfg_ana.method = 'elasNetLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0.5,'lambda',1,'family','binomial')};
% terrible performance

% optimize lambda using a 5-fold inner cross validation
%cfg_ana.method = 'optLamLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0.5, 'validator',ft_mv_crossvalidator('nfolds',5,'metric','accuracy'),'family','binomial')};
%cfg_ana.method = 'L2optLamLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0, 'validator',ft_mv_crossvalidator('nfolds',5,'metric','accuracy'),'family','binomial')};

for sub = 1:length(exper.subjects)
  fprintf('%s\n',exper.subjects{sub});
  
  data1 = data_erp.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data;
  data2 = data_erp.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data;
  
  cfg.design = [ones(size(data1.trial,1),1); 2*ones(size(data2.trial,1),1)]';
  
  stat = ft_timelockstatistics(cfg,data1,data2);
  
  accuracy(sub) = stat.statistic.accuracy;
  binomial(sub) = stat.statistic.binomial;
  %continTab{sub} = stat.cv.performance('contingency');
  continTab{sub} = stat.statistic.contingency;
  
  % print out some performance info
  fprintf('\n%s\naccuracy = %.2f%%, p = %.4f',exper.subjects{sub},accuracy(sub)*100,binomial(sub));
  if binomial(sub) < .05
    fprintf(' ***\n');
  else
    fprintf('\n');
  end
  
  fprintf('\t\tPredicted\n');
  fprintf('\t\t%s\t%s\n',ana.eventValues{typ}{1},ana.eventValues{typ}{2});
  fprintf('True\t%s\t%d\t%d\n',ana.eventValues{typ}{1},continTab{sub}(1,1),continTab{sub}(1,2));
  fprintf('\t%s\t%d\t%d\n',ana.eventValues{typ}{2},continTab{sub}(2,1),continTab{sub}(2,2));
  
end

%% ERP: save the results to a file

%vsStr = sprintf('%s%s',ana.eventValues{typ}{1},sprintf(repmat('vs%s',1,1),ana.eventValues{typ}{2}));
vsStr = sprintf(repmat('%s',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:});
chanStr = sprintf(repmat('%s_',1,length(cfg_ana.chanStr)),cfg_ana.chanStr{:});
if strcmp(cfg.avgovertime,'yes')
  timeAvgStr = 'Avg';
else
  timeAvgStr = '';
end
if strcmp(cfg.avgoverchan,'yes')
  chanAvgStr = 'Avg';
else
  chanAvgStr = '';
end
outfile = sprintf('ft_%s_%s_%dfold_%dsub_%s%s%d_%d%s.txt',...
  vsStr,cfg_ana.method,cfg.nfolds,length(exper.subjects),...
  chanStr,chanAvgStr,cfg.latency(1)*1000,cfg.latency(2)*1000,timeAvgStr);
fprintf('Saving results to %s...\n',outfile);
fid = fopen(fullfile(dataroot,outfile),'w+');

% print details and performance
fprintf(fid,'%s\n',cfg_ana.method);
fprintf(fid,'%s\t%d folds\n',cfg.method,cfg.nfolds);
fprintf(fid,'start: %d ms\n',cfg.latency(1)*1000);
fprintf(fid,'end: %d ms\n',cfg.latency(2)*1000);

fprintf(fid,'avgOverChan\t%s\n',cfg.avgoverchan);
fprintf(fid,'avgOverTime\t%s\n',cfg.avgovertime);
fprintf(fid,'\n');

fprintf(fid,'Probabilities\n');
%fprintf(fid,'\tAverage\n');
for sub = 1:length(exper.subjects)
  fprintf(fid,'%s%s',exper.subjects{sub},sprintf(repmat('\t%.2f',1,size(accuracy,2)),accuracy(sub,:)));
  %fprintf(fid,'\t%.2f\n',mean(accuracy(sub,:),2));
  fprintf(fid,'\n');
end
fprintf(fid,'Average%s',sprintf(repmat('\t%.2f',1,size(accuracy,2)),mean(accuracy,1)));
%fprintf(fid,'\t%.2f\n',mean(mean(accuracy,1),2));
fprintf(fid,'\n');

% print p-values
fprintf(fid,'p-values\n');
for sub = 1:length(exper.subjects)
  fprintf(fid,'%s%s\n',exper.subjects{sub},sprintf(repmat('\t%f',1,size(binomial,2)),binomial(sub,:)));
end

fclose(fid);

% %% plot it
%
% cfg = [];
% cfg.interplimits = 'electrodes';
% cfg.interpolation = 'linear';
% cfg.layout = 'GSN-HydroCel-129.sfp';
% cfg.xlim = cfg_ana.latency;
% cfg.parameter = 'model1';
% cfg.comments = '';
% cfg.colorbar = 'yes';
% ft_topoplotER(cfg,stat);






%% Time-Frequency: process the data

%subInd = find(ismember(exper.subjects,cfg_prep.subjs));
ses = 1;

if equatetrials
  ana.equateTrials = true;
  
  % initialize random number generator
  rng('shuffle','twister');
  
  for typ = 1:length(ana.eventValues)
    % find the lowest number of trials within each subject
    lowEvNum = Inf(length(exper.subjects),1);
    for sub = 1:length(exper.subjects)
      for i = 1:length(ana.eventValues{typ})
        if size(data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trial,2) < lowEvNum(sub)
          lowEvNum(sub) = size(data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trial,2);
        end
      end
    end
  end
else
  ana.equateTrials = false;
end

% initialize the struct
data_freq = struct;

cfg_ana = [];
%cfg_ana.chanStr = {'right'};
%cfg_ana.chanStr = {'LPS'};
%cfg_ana.chanStr = {'center74'};
cfg_ana.chanStr = {'center91'};
%cfg_ana.chanStr = {'all129'};

% set the analysis parameters
cfg = [];
%cfg.parameter = 'trial';
cfg.keeptrials = 'yes'; % classifiers operate on individual trials
cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)}));

% % multitaper
% cfg.output       = 'pow';
% cfg.method       = 'mtmconvol';
% %cfg.foi          = (4:1:64);
% % logarythmically spaced (4 to 128)
% cfg.foi = (2^(1/8)).^(16:56);
% % % logarythmically spaced (2 to 64)
% %cfg.foi = (2^(1/8)).^(16:48);
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
% %cfg.toi          = (0.3:0.02:1.0);
% %cfg.toi          = (0.1:0.02:0.6);
% cfg.toi          = (-0.3:0.02:1.0);
% %cfg.taper        = 'hanning';
% cfg.taper = 'dpss';
% % tapsmofrq is not used for hanning taper; it is used for dpss
% cfg.tapsmofrq = 0.4*cfg.foi;

% wavelet
cfg.method = 'wavelet';
cfg.width = 6;
cfg.toi = (-0.3:0.02:1.625);
%cfg.foi = (4:1:64);
cfg.foi = (2^(1/8)).^(16:53);

% logarythmically spaced (4 to 100)
% cfg.foi = (2^(1/8)).^(16:53);
% logarythmically spaced (2 to 64)
% cfg.foi = (2^(1/8)).^(8:48);

% get the TF data
for typ = 1:length(ana.eventValues)
  for sub = 1:length(exper.subjects)
    for i = 1:length(ana.eventValues{typ})
      if equatetrials
        ana.equateTrials = 'yes';
        % choose a random subset to equate trial numbers
        evNums = randperm(size(data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trial,2));
        cfg.trials = sort(evNums(1:lowEvNum(sub)));
      end
      data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data = ft_freqanalysis(cfg,data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data);
    end
  end
end


%% Add trials to simulate skewed classes.

if addtrials
  
  ana.addtrials = true;
  
  ses=1;
  
  % initialize random number generator
  rng('shuffle','twister');
  
  for typ = 1:length(ana.eventValues)
    % find the highest number of trials within each subject
    highEvNum = zeros(length(exper.subjects),1);
    
    for sub = 1:length(exper.subjects)
      for i = 1:length(ana.eventValues{typ})
        if size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1) > highEvNum(sub)
          highEvNum(sub) = size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1);
        end
      end
      for i = 1:length(ana.eventValues{typ})
        % if we're working on the smaller event, get a random sub-selection
        if size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1) == highEvNum(sub)
            
            
            data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm((highEvNum(sub)+1):(highEvNum(sub)*2),:,:,:) = data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm(1:(highEvNum(sub)),:,:,:);
            data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.cumtapcnt((highEvNum(sub)+1):(highEvNum(sub)*2),:) = data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.cumtapcnt(1:(highEvNum(sub)),:);
            data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trialinfo((highEvNum(sub)+1):(highEvNum(sub)*2),:) = data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trialinfo(1:(highEvNum(sub)),:);
    
%           fprintf('%s: Subselecting %d (of %d) trials for %s...\n',exper.subjects{sub},lowEvNum(sub),size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1),ana.eventValues{typ}{i});
%           evNums = randperm(size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1));
%           % old version
%           data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data = ft_selectdata(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data,'rpt',sort(evNums(1:lowEvNum(sub))));
%           % new version, rpt selection doesn't work as of ft-20120116
%           %cfg = [];
%           %cfg.rpt = sort(evNums(1:lowEvNum(sub)));
%           %data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data = ft_selectdata(cfg,data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data);
        end
      end % i
    end % sub
    
  end % typ
else
  ana.addtrials = false;
end

%% subtract auditory trials to simulate skewed classes

if subtract_trials 
  
  ana.subtract_trials = true;
  
  ses=1;
  
  % initialize random number generator
  rng('shuffle','twister');
  
  typ = 1;  %auditory 
  
  lowEvNum = zeros(length(exper.subjects));
    
    for sub = 1:length(exper.subjects)
          lowEvNum(sub) = size(data_freq.(ana.eventValues{1}{2}).sub(sub).ses(ses).data.trialinfo,1)/2;%set low number equal to half the visual trials
           
%         lowEvNum(sub) = size(data_freq.Vis.sub.ses.data.trialinfo,1)/2;%set low number equal to half the visual trials
       
          fprintf('%s: Subselecting %d (of %d) trials for %s...\n',exper.subjects{sub},lowEvNum(sub),size(data_freq.(ana.eventValues{1}{1}).sub(sub).ses(ses).data.powspctrm,1),ana.eventValues{1}{1});
          evNums = randperm(size(data_freq.(ana.eventValues{1}{1}).sub(sub).ses(ses).data.powspctrm,1));
          data_freq.(ana.eventValues{1}{1}).sub(sub).ses(ses).data = ft_selectdata(data_freq.(ana.eventValues{1}{1}).sub(sub).ses(ses).data,'rpt',sort(evNums(1:lowEvNum(sub))));
          
    end % sub
    
else
  ana.equateTrials = false;
end


%% Time-Frequency: if data is already processed but not equated

% TODO: can happen automatically in ft_mv_crossvalidator ('balance',true)

if equatetrials
  
  ana.equateTrials = true;
  
  ses=1;
  
  % initialize random number generator
  rng('shuffle','twister');
  
  for typ = 1:length(ana.eventValues)
    % find the lowest number of trials within each subject
    lowEvNum = Inf(length(exper.subjects),1);
    
    for sub = 1:length(exper.subjects)
      for i = 1:length(ana.eventValues{typ})
        if size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1) < lowEvNum(sub)
          lowEvNum(sub) = size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1);
        end
      end
      for i = 1:length(ana.eventValues{typ})
        % if we're working on the smaller event, get a random sub-selection
        if size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1) > lowEvNum(sub)
          fprintf('%s: Subselecting %d (of %d) trials for %s...\n',exper.subjects{sub},lowEvNum(sub),size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1),ana.eventValues{typ}{i});
          evNums = randperm(size(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.powspctrm,1));
          % old version
          data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data = ft_selectdata(data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data,'rpt',sort(evNums(1:lowEvNum(sub))));
          % new version, rpt selection doesn't work as of ft-20120116
          %cfg = [];
          %cfg.rpt = sort(evNums(1:lowEvNum(sub)));
          %data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data = ft_selectdata(cfg,data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data);
        end
      end % i
    end % sub
    
  end % typ
else
  ana.equateTrials = false;
end

%% Time-Frequency: Change in freq relative to baseline using absolute power

% save the original data
%data_freq_orig = data_freq;

cfg_fb = [];
cfg_fb.baseline = [-0.3 -0.1];
cfg_fb.baselinetype = 'absolute';
%cfg_fb.baselinetype = 'relative';

ana.blc = true;
ana.blc_method = cfg_fb.baselinetype;

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        fprintf('%s, %s, %s, ',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
        data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
      end
    end
  end
end
fprintf('Done.\n');

%% Time-Frequency: classification

% TODO:
%
% use ft_selectdata (old or new?) to pull the data out of the ft struct
%
% use DMLT's functions to do more nuanced classification than FT can do
%
% understand the binomial test results
%
% low and high alpha? what were the differences?
%
% normalized power? (use my mm_ft_loadData)

if ~isfield(ana,'blc')
  ana.blc = false;
end
if ~isfield(ana,'equateTrials')
  ana.equateTrials = false;
end

ses = 1;

cfg_ana = [];

%cfg_ana.freqs = [4 8; 6 10; 8 12; 10 15; 12 19; 15 25; 19 30; 25 35; 30 40; 35 64; 64 100];
%cfg_ana.freqs = [4 8; 6 10; 8 12; 10 15; 12 19; 15 25; 19 30; 25 35; 30 40; 35 64];
%cfg_ana.freqs = [4 8; 8.1 12; 12.1 19; 19.1 30; 30.1 42; 42.1 64];
%cfg_ana.freqs = [4 8; 8.1 14; 14.1 28; 28.1 42; 42.1 64];
%cfg_ana.freqs = [4 8; 6 12; 8.1 14; 12 18; 14.1 22; 18.1 28; 24.1 36; 28.1 42; 36.1 50; 42.1 64];
cfg_ana.freqs = [4 8; 8.1 14; 14.1 21; 21.1 28; 28.1 42; 42.1 64];
%cfg_ana.freqs = [4 64];
%cfg_ana.freqs = [4 8; 8 14; 14 28; 28 64];
%cfg_ana.freqs = [4 8; 8.1 14; 14.1 32; 32.1 64; 64.1 100];

%cfg_ana.latencies = [0 1.0];
cfg_ana.latencies = [0 0.5; 0.5 1.0];
%cfg_ana.latencies = [1.0 1.45];
%cfg_ana.latencies = [0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
%cfg_ana.latencies = [-0.2 0.0; 0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];

%cfg_ana.chanStr = {{'right'}};
%cfg_ana.chanStr = {{'center74'}};
cfg_ana.chanStr = {{'center91'}};
% cfg_ana.chanStr = {{'center74'},{'center91'}};
%cfg_ana.chanStr = {{'center74'},{'center91'},{'midline'},{'left'},{'right'},{'anterior'},{'posterior'}};
%cfg_ana.chanStr = {{'LAS'},{'RAS'},{'LPS'},{'RPS'},{'LAI'},{'RAI'},{'LPI'},{'RPI'},{'FI'},{'FS'},{'PS'},{'PI'}};
%cfg_ana.chanStr = {{'center74'},{'center91'},{'midline'},{'left'},{'right'},{'anterior'},{'posterior'},{'LAS'},{'RAS'},{'LPS'},{'RPS'},{'LAI'},{'RAI'},{'LPI'},{'RPI'},{'FI'},{'FS'},{'PS'},{'PI'}};
%cfg_ana.chanStr = {{'all129'}};

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.method = 'crossvalidate';
% cfg.nfolds = 10;
cfg.nfolds = 5;
%cfg.nfolds = Inf;
%cfg.nfolds = 0.8;

cfg.statistic = {'accuracy' 'binomial' 'contingency'};
%cfg.statistic = {'accuracy' 'binomial'};

%cfg.latency = [-0.2 0];
%cfg.latency = [0 1.0];
%cfg.latency = [0 0.6];
%cfg.latency = [0 0.3];
%cfg.latency = [0.3 0.6];
%cfg.latency = [0.5 1.0];

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
%cfg.avgoverfreq = 'no';
cfg.avgoverfreq = 'yes';

if length(cfg_ana.chanStr) == 1
  cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr{1})}));
  cfg_ana.nChan = length(cfg.channel);
end

% z-transform, feature selection (map data to a different space), SVM
%cfg_ana.method = 'cspFs10SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_filterer('maxfeatures',10) ft_mv_svm};

%cfg_ana.method = 'cspSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_svm};
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan,'numpatterns',3,'outputdatatype','logpowcsp','filttype','CSP0') ft_mv_svm};

%cfg_ana.method = 'filt100SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_filterer('maxfeatures',100) ft_mv_svm('verbose',true)};

%cfg_ana.method = 'CSPgsCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_gridsearch('verbose',true,'mva',ft_mv_svm,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars','C','vals',logspace(-3,3,7))};

%cfg_ana.method = 'gsFiltCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_gridsearch('verbose',true,'mva',{ft_mv_filterer ft_mv_svm},'validator',ft_mv_crossvalidator('nfolds',4,'metric','accuracy'),'mvidx',[1 2],'vars',{'maxfeatures' 'C'},'vals',{1:2:10 [1 10]})};
%cfg.mva = {ft_mv_standardizer ft_mv_gridsearch('verbose',true,'mva',{ft_mv_filterer ft_mv_svm},'validator',ft_mv_crossvalidator('nfolds',4,'metric','accuracy'),'mvidx',[1 2],'vars',{'maxfeatures' 'C'},'vals',{1:2:10 logspace(-3,3,7)})};

%cfg_ana.method = 'filt100gsCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_filterer('maxfeatures',100) ft_mv_gridsearch('verbose',true,'mva',ft_mv_svm,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars','C','vals',logspace(-3,3,7))};
%cfg_ana.method = 'gsCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_gridsearch('verbose',true,'mva',ft_mv_svm,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars','C','vals',logspace(-3,2,6))};

% z-transform, feature selection in the original space
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0.1)};

% basic SVM
%cfg_ana.method = 'SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_svm};

% unregularized linear regression
%cfg_ana.method = 'unregLinR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0,'family','gaussian')};
% Didn't work

% unregularized logistic regression
%cfg_ana.method = 'unregLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0,'family','binomial')};
% crashes

%cfg_ana.method = 'L1gsLam';
%cfg.mva = {ft_mv_standardizer ft_mv_gridsearch('verbose',true,'mva',ft_mv_glmnet,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars',{'alpha','lambda'},'vals',{(0:0.1:1) (0.1:0.01:0.2)})};
% not working

% L1 regularized logistic regression
% cfg_ana.method = 'L1regLogRLam01';
% cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',1,'lambda',0.1,'family','binomial')};
%cfg_ana.method = 'L1regLogRLam01_ns';
%cfg.mva = {ft_mv_glmnet('alpha',1,'lambda',0.1,'family','binomial')};
% lam=1 doesn't perform very well
% lam=0.1 and 0.05 performs well
% lambda=0.2 doesn't perform well

% cfg.mva = {dml.standardizer dml.graphnet('family','binomial','L1',0.1)};
cfg.mva = {dml.standardizer dml.enet('family','binomial','alpha',1)};
% cfg.mva = {dml.standardizer dml.glmnet('family','binomial','alpha',1)};
cfg_ana.method = 'L1regLog-av_sub-skew_resample';

%cfg_ana.method = 'L1regLogRLamN100';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',1,'lambda',[],'nlambda',100,'lambda_min',0.05,'family','binomial')};
% not working

% L2 regularized logistic regression
%cfg_ana.method = 'L2regLogRLam09';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0,'lambda',0.9,'family','binomial')};
%cfg_ana.method = 'L2regLogRLam1';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0,'lambda',1,'family','binomial')};
% alpha=0, lambda=1 has good performance
% alpha=0, lambda=0.1 ? (Slooow)
% lam 0.9 doesn't perform much better than 1.0

% elastic net logistic regression
%cfg_ana.method = 'elasNetLogRAlp99Lam01';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0.99,'lambda',0.1,'family','binomial')};
%cfg_ana.method = 'elasNetLogRAlp01Lam09';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0.1,'lambda',0.9,'family','binomial')};
% meh performance at alpha=0, lambda=2
% decent performance at alpha=0.1, lambda=0.9
% bad performance at alpha=0.2, lambda=0.9
% bad performance at alpha=0.1, lambda=0.5
% terrible performance at alpha=0.5, lambda=1
% great performance at alpha=0.99, lambda=0.1

% L1 optimize lambda using a 5-fold inner cross validation
%cfg_ana.method = 'L1regLogRLamOpt';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',1,'validator',ft_mv_crossvalidator('verbose',true,'nfolds',4,'metric','accuracy'),'family','binomial')};

%cfg_ana.method = 'elasNetLogRAlp99LamOpt';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0.99,'validator',ft_mv_crossvalidator('verbose',true,'nfolds',4,'metric','accuracy'),'family','binomial')};

% L2 optimize lambda using a 5-fold inner cross validation
%cfg_ana.method = 'L2regLogRLamOpt';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0,'validator',ft_mv_crossvalidator('verbose',true'nfolds',5,'metric','accuracy'),'family','binomial')};
% lumping all trials in one category, really bad performance


% %%%%%%%%% testing
% cfg_ana.method = 'L1regLogRLam01';
% obj = ft_mv_analysis({ft_mv_standardizer ft_mv_glmnet('alpha',1,'lambda',0.1,'family','binomial')});
% 
% obj = obj.train(X,Y);
% 
% %%%%%%%%%%%%%%%%%
% 
% gn = ft_mv_glmnet('validator',ft_mv_crossvalidator('nfolds',5),'family','gaussian');
% 
% cv = ft_mv_crossvalidator('mva',{ft_mv_standardizer ft_mv_svm},'nfolds',10,'metric','accuracy');
% cv = cv.train(X,Y);
% % display classification accuracy
% cv.metric = 'accuracy';
% disp(cv.performance);
% % display outcome of mcnemar test
% cv.sigtest = 'mcnemar';
% disp(cv.significance);
% 
% % investigate 'trainfolds' (or 'testfolds') property
% 
% %%%%%%%%%%%%%%%%
% 
% ft_mv_hmm;
% %%%%%%%%%%%%%%%%%%%%%


for typ = 1:length(ana.eventValues)
  if length(ana.eventValues{typ}) ~= 2
    error('Currently we can only classify 2 events, you included %d (%s)',length(ana.eventValues{typ}),vsStr);
  end
  vsStr = sprintf(repmat('%s',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:});
  vsStrTab = sprintf(repmat('\t%s',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:});
  
  % initialize some values
  accuracy = nan(length(exper.subjects),length(cfg_ana.chanStr),size(cfg_ana.latencies,1),size(cfg_ana.freqs,1));
  binomial = nan(length(exper.subjects),length(cfg_ana.chanStr),size(cfg_ana.latencies,1),size(cfg_ana.freqs,1));
  continTab = cell(length(exper.subjects),length(cfg_ana.chanStr),size(cfg_ana.latencies,1),size(cfg_ana.freqs,1));
  stat_all = struct([]);
  
  for chn = 1:length(cfg_ana.chanStr)
    cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr{chn})}));
    cfg_ana.nChan = length(cfg.channel);
    chanStr = sprintf(repmat('%s_',1,length(cfg_ana.chanStr{chn})),cfg_ana.chanStr{chn}{:});
    chanStrTab = sprintf(repmat('\t%s',1,length(cfg_ana.chanStr{chn})),cfg_ana.chanStr{chn}{:});
    if isempty(cfg.channel)
      error('No channels selected for%s!',chanStrTab);
    end
    
    for sub = 1:length(exper.subjects)
      fprintf('%s\n',exper.subjects{sub});
      data1 = data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data;
      data2 = data_freq.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data;
      
      for lat = 1:size(cfg_ana.latencies,1)
        cfg.latency = cfg_ana.latencies(lat,:);
        
        for frq = 1:size(cfg_ana.freqs,1)
          cfg.frequency = cfg_ana.freqs(frq,:);
          
          data = ft_selectdata_old(data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data,data_freq.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data,...
            'param',cfg.parameter,...
            'foilim',cfg_ana.freqs(frq,:),...
            'toilim',cfg_ana.latencies(lat,:),...
            'channel',cfg.channel,...
            'avgoverchan',cfg.avgoverchan,...
            'avgoverfreq',cfg.avgoverfreq,...
            'avgovertime',cfg.avgovertime);
          
          %data = ft_selectdata_new(cfg,data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data,data_freq.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data);
          
          if any(isinf(data.(cfg.parameter)(:)))
            warning('Inf encountered; replacing by zeros');
            data.(cfg.parameter)(isinf(data.(cfg.parameter)(:))) = 0;
          end
          
          if any(isnan(data.(cfg.parameter)(:)))
            warning('Nan encountered; replacing by zeros');
            data.(cfg.parameter)(isnan(data.(cfg.parameter)(:))) = 0;
          end
          
          % make the data 2D
          cfg.design = [ones(size(data1.(cfg.parameter),1),1); 2*ones(size(data2.(cfg.parameter),1),1)]';
          reshapevec = [size(data.(cfg.parameter),1), (size(data.(cfg.parameter),2) * size(data.(cfg.parameter),3) * size(data.(cfg.parameter),4))];
          
          % % do we need to make the data 3D?
          % cfg.design = repmat(cfg.design',[1,1,size(data.time,2)]);
          % reshapevec = [size(data.(cfg.parameter),1), (size(data.(cfg.parameter),2) * size(data.(cfg.parameter),3)), size(data.(cfg.parameter),4)];
          
          % debug
          keyboard
          
          X = reshape(data.(cfg.parameter),reshapevec);
          Y = cfg.design';
          
          m = dml.standardizer;
          m = m.train(X);
          Z = m.test(X);
          
          % search through L1 parameters
          v = dml.graphnet.lambdapath(Z,Y,'binomial',50,1e-2);
          m = dml.gridsearch('validator',dml.crossvalidator('type','split','stat','accuracy','mva',{dml.standardizer dml.graphnet('family','binomial','restart',false)}),'vars','L1','vals',v,'verbose',true);
          tic;
          m = m.train(X,Y);
          toc
          %m.statistic
          figure;
          plot(m.configs,m.outcome); xlabel('L1'); ylabel('accuracy');
          Z2 = m.test(X);
          
          m = dml.crossvalidator('mva',{dml.standardizer dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'verbose',true);
          m = m.train(X,Y);
          Z2 = m.test(X);
          
          % gridsearch SVM for best C
          m = dml.gridsearch('validator',dml.crossvalidator('type','split','stat','accuracy','mva',{dml.standardizer dml.svm('restart',false)}),'vars','C','vals',fliplr(logspace(-4,1,6)),'verbose',true);
          m = m.train(X,Y);
          Z2 = m.test(X);
          %m.statistic % doesn't work for gridsearch, not sure how to evaluate
          
          % Naive Bayes
          m = dml.crossvalidator('mva',{dml.standardizer dml.naive},'stat',{'accuracy','binomial','contingency'},'verbose',true);
          m = m.train(X,Y);
          m.statistic
          
          % simple elastic net
          m = dml.graphnet('family','binomial','L1',0.1);
          m = m.train(X,Y);
          m.statistic
          
          % faster elastic net; no params: finds best lambda
          m = dml.enet;
          tic;
          m = m.train(X,Y);
          toc;
          figure;
          bar(m.model.weights);
          figure;
          imagesc(data.time,1:size(data.label,1),squeeze(reshape(m.model.weights,[size(data.(cfg.parameter),2), size(data.(cfg.parameter),3), size(data.(cfg.parameter),4)])));
          axis xy
          
          % resample Naive Bayes
          m = dml.crossvalidator('mva',{dml.standardizer dml.naive},'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true);
          m = m.train(X,Y);
          m.statistic
          
          % same as fieldtrip defaults
          m = dml.crossvalidator('mva',dml.analysis({dml.standardizer('verbose',true) dml.svm('verbose',true)}),...
            'stat',{'accuracy','binomial','contingency'},...
            'type','nfold','folds',5);
          m = m.train(X,Y);
          m.statistic
          
          % run it
          stat_all(sub,chn,lat,frq).stat = ft_freqstatistics(cfg,data1,data2);
          
          % store some values
          accuracy(sub,chn,lat,frq) = stat_all(sub,chn,lat,frq).stat.statistic.accuracy;
          binomial(sub,chn,lat,frq) = stat_all(sub,chn,lat,frq).stat.statistic.binomial;
          continTab{sub,chn,lat,frq} = stat_all(sub,chn,lat,frq).stat.statistic.contingency;
          
          % print out some performance info
          fprintf('\n%s, %s, %s (%d chans), %.1f-%.1f s, %.1f-%.1f Hz\n',exper.subjects{sub},vsStr,chanStr(1:end-1),cfg_ana.nChan,cfg.latency(1),cfg.latency(2),cfg.frequency(1),cfg.frequency(2));
          fprintf('accuracy = %.2f%%',stat_all(sub,chn,lat,frq).stat.statistic.accuracy*100);
          fprintf(', binomial = %.4f',stat_all(sub,chn,lat,frq).stat.statistic.binomial);
          if stat_all(sub,chn,lat,frq).stat.statistic.binomial < .05
            fprintf(' ***\n');
          else
            fprintf('\n');
          end
          
          % print out the contingency table
          fprintf('\t\tPredicted\n');
          fprintf('\t\t%s\t%s\n',ana.eventValues{typ}{1},ana.eventValues{typ}{2});
          fprintf('True\t%s\t%d\t%d\n',ana.eventValues{typ}{1},continTab{sub,chn,lat,frq}(1,1),continTab{sub,chn,lat,frq}(1,2));
          fprintf('\t%s\t%d\t%d\n',ana.eventValues{typ}{2},continTab{sub,chn,lat,frq}(2,1),continTab{sub,chn,lat,frq}(2,2));
          fprintf('\n');
          
        end % frq
      end % lat
    end % sub
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time-Frequency: save the results to a file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(cfg.avgoverfreq,'yes')
      freqAvgStr = 'Avg';
    else
      freqAvgStr = '';
    end
    if strcmp(cfg.avgovertime,'yes')
      timeAvgStr = 'Avg';
    else
      timeAvgStr = '';
    end
    if size(cfg_ana.latencies,1) > 1
      timeAvgStr = sprintf('%sSplit%d',timeAvgStr,size(cfg_ana.latencies,1));
    end
    if strcmp(cfg.avgoverchan,'yes')
      chanAvgStr = 'Avg';
    else
      chanAvgStr = '';
    end
    
    outfile = sprintf('ft_%s_%s_%dfreq%s_%s%s%d_%d%s',...
      vsStr,cfg_ana.method,size(cfg_ana.freqs,1),freqAvgStr,...
      chanStr,chanAvgStr,cfg_ana.latencies(1,1)*1000,cfg_ana.latencies(end,end)*1000,timeAvgStr);
    
    fprintf('Saving results to %s...\n',outfile);
    fid = fopen(fullfile(dataroot,[outfile,'.txt']),'w+');
    
    fprintf(fid,'ROI%s\t(%d channels)\n',chanStrTab,cfg_ana.nChan);
    fprintf(fid,'Events%s\n',vsStrTab);
    
    % print details and performance
    fprintf(fid,'%s\n',cfg_ana.method);
    
    fprintf(fid,'tf method\t%s\n',ft_findcfg(data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.cfg,'method'));
    fprintf(fid,'taper\t%s\n',ft_findcfg(data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.cfg,'taper'));
    
    fprintf(fid,'all freqs\t%.2f-%.2f Hz\n',data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.freq(1),data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.freq(end));
    fprintf(fid,'all times\t%.2f-%.2f s\n',data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.time(1),data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.time(end));
    fprintf(fid,'this freq range\t%.2f-%.2f Hz\n',cfg_ana.freqs(1),cfg_ana.freqs(end));
    fprintf(fid,'this freq bins%s\n',sprintf(repmat('\t%.2f-%.2f Hz',1,size(cfg_ana.freqs,1)),cfg_ana.freqs'));
    fprintf(fid,'this time range\t%.2f-%.2f s\t%s\n',cfg_ana.latencies(1,1),cfg_ana.latencies(end,end),timeAvgStr);
    fprintf(fid,'this time bins%s\n',sprintf(repmat('\t%.2f-%.2f s',1,size(cfg_ana.latencies,1)),cfg_ana.latencies'));
    
    if isfield(ana,'blc')
      if ana.blc
        if ~isfield(ana,'blc_method')
          ana.blc_method = 'unknown';
        end
        fprintf(fid,'tf baseline corrected\t%s\n',ana.blc_method);
      else
        fprintf(fid,'not tf baseline corrected\n');
      end
    else
      fprintf(fid,'probably not tf baseline corrected\n');
    end
    
    if isfield(ana,'equateTrials')
      if ana.equateTrials
        fprintf(fid,'equated trial counts\n');
      else
        fprintf(fid,'non-equated trial counts\n');
      end
    else
      fprintf(fid,'probably non-equated trial counts\n');
    end
    
    fprintf(fid,'mva\t%s\n',cfg_ana.method);
    fprintf(fid,'method\t%s\n',cfg.method);
    fprintf(fid,'folds\t%.2f\n',cfg.nfolds);
    
    fprintf(fid,'avgOverChan\t%s\n',cfg.avgoverchan);
    fprintf(fid,'avgOverTime\t%s\n',cfg.avgovertime);
    fprintf(fid,'avgOverFreq\t%s\n',cfg.avgoverfreq);
    fprintf(fid,'\n');
    
    fprintf(fid,'Probabilities\n');
    fprintf(fid,'\t%s',sprintf(repmat('\t%.1f-%.1f Hz',1,size(cfg_ana.freqs,1)),cfg_ana.freqs'));
    fprintf(fid,'\tAverage\n');
    for sub = 1:length(exper.subjects)
      fprintf(fid,'%s',exper.subjects{sub});
      for lat = 1:size(cfg_ana.latencies,1)
        % print the accuracy for eqach frequency band for this latency
        fprintf(fid,'\t%.1f-%.1fs%s',cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2),sprintf(repmat('\t%.3f',1,size(cfg_ana.freqs,1)),accuracy(sub,chn,lat,:)));
        % average across frequency bands
        fprintf(fid,'\t%.3f\n',mean(accuracy(sub,chn,lat,:),4));
      end
    end
    fprintf(fid,'Average');
    for lat = 1:size(cfg_ana.latencies,1)
      % average acorss subjects within this latency
      fprintf(fid,'\t%.1f-%.1fs%s',cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2),sprintf(repmat('\t%.3f',1,size(cfg_ana.freqs,1)),mean(accuracy(:,chn,lat,:),1)));
      % average across frequency bands and subjects
      fprintf(fid,'\t%.3f\n',mean(mean(accuracy(:,chn,lat,:),4),1));
    end
    % average across all latencies and subjects
    fprintf(fid,'\tAverage%s',sprintf(repmat('\t%.3f',1,size(cfg_ana.freqs,1)),mean(mean(accuracy(:,chn,:,:),3),1)));
    % average across everything
    fprintf(fid,'\t%.2f\n',mean(mean(mean(accuracy(:,chn,:,:),1),3),4));
    fprintf(fid,'\n');
    
    % print p-values
    fprintf(fid,'binomial test values\n');
    fprintf(fid,'\t%s\n',sprintf(repmat('\t%.1f-%.1f Hz',1,size(cfg_ana.freqs,1)),cfg_ana.freqs'));
    for sub = 1:length(exper.subjects)
      fprintf(fid,'%s',exper.subjects{sub});
      for lat = 1:size(cfg_ana.latencies,1)
        % p-values for each subject, latency, and frequnecy band
        fprintf(fid,'\t%.1f-%.1fs%s\n',cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2),sprintf(repmat('\t%f',1,size(cfg_ana.freqs,1)),binomial(sub,chn,lat,:)));
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % statistics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % t-test to see if accuracies across subjects are above 50%. within a
    % frequency band, within a time region, test accuracy vs 0.5.
    
    if length(exper.subjects) > 1
      chance50 = 0.5*ones(size(exper.subjects));
      alpha = 0.05;
      
      fprintf(fid,'\nVersus chance\n');
      fprintf(fid,'%s\n',sprintf(repmat('\t%.1f-%.1f Hz\t',1,size(cfg_ana.freqs,1)),cfg_ana.freqs'));
      fprintf(fid,'%s\n',sprintf(repmat('\tp\tt(%d)',1,size(cfg_ana.freqs,1)),(size(exper.subjects,1)-1)*ones(size(cfg_ana.freqs,1),1)));
      
      for lat = 1:size(cfg_ana.latencies,1)
        fprintf(fid,'%.1f-%.1f s',cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2));
        for frq = 1:size(cfg_ana.freqs,1)
          [h,p,ci,stats] = ttest(squeeze(accuracy(:,chn,lat,frq)),chance50,alpha,'both');
          fprintf('%s, %.1f-%.1f s,\t%.1f-%.1f Hz,\tavg=%.2f%%:\t',chanStr(1:end-1),cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2),cfg_ana.freqs(frq,1),cfg_ana.freqs(frq,2),mean(accuracy(:,chn,lat,frq),1)*100);
          fprintf('p=%.4f, t(%d)=%.4f',p,stats.df,stats.tstat);
          if h == 1
            fprintf(' <---');
          end
          fprintf('\n');
          
          fprintf(fid,'\t%.4f\t%.4f',p,stats.tstat);
          
        end
        fprintf('\n');
        fprintf(fid,'\n');
      end
    else
      fprintf(fid,'Only 1 subject, cannot compare to chance-level classification performance.\n');
      fprintf('Only 1 subject, cannot compare to chance-level classification performance.\n');
    end
    
    % close the file
    fclose(fid);
    
    % save the details to a mat file
    save(fullfile(dataroot,outfile),'stat_all');
    
    fprintf('Done\n\n');
    
  end % channel ROI
end % typ

%% Time-Frequency: make feature subplots of each subject's frequency bands

clim = [-0.05, 0.05];
figsize = [900 900];

cfg_plot = [];
cfg_plot.layout = 'GSN-HydroCel-129.sfp';
cfg_plot.parameter = 'model1';
cfg_plot.colorbar = 'no';
cfg_plot.zlim = [-0.0005 0.0005];
cfg_plot.marker = 'on';
%cfg_plot.marker = 'labels';
cfg_plot.fontsize = 9;

cfg_plot.plotstyle = 'topo';
%cfg_plot.plotstyle = 'multi';
%cfg_plot.plotstyle = 'F-TxC';
%cfg_plot.plotstyle = 'C-TxF';

for chn = 1:length(cfg_ana.chanStr)
  for sub = 1:length(exper.subjects)
    for lat = 1:size(cfg_ana.latencies,1)
      for frqband = 1:size(cfg_ana.freqs,1)
        
        % only plot if the p-value was reasonable
        if stat_all(sub,chn,lat,frqband).stat.statistic.binomial < 0.05
          if strcmp(cfg_plot.plotstyle,'topo') || strcmp(cfg_plot.plotstyle,'F-TxC')
            figure;
            count = 0;
            nFreq = length(stat_all(sub,chn,lat,frqband).stat.freq);
            nChan = length(stat_all(sub,chn,lat,frqband).stat.label);
            
            for frq = 1:nFreq
              % Model 1
              count = count + 1;
              
              if strcmp(cfg_plot.plotstyle,'topo')
                if nFreq <= 5
                  nRow = 1;
                  nCol = nFreq;
                elseif nFreq > 5 && nFreq <= 10
                  nRow = 2;
                  nCol = ceil(nFreq / 2);
                elseif nFreq > 10 && nFreq <= 15
                  nRow = 3;
                  nCol = ceil(nFreq / 3);
                elseif nFreq > 15
                  nRow = 4;
                  nCol = ceil(nFreq / 4);
                end
                
                cfg_plot.comment = sprintf('%s, %.1f Hz',exper.subjects{sub},stat_all(sub,chn,lat,frqband).stat.freq(frq));
                subplot(nRow,nCol,count);
                ft_topoplotTFR(cfg_plot,stat_all(sub,chn,lat,frqband).stat);
              elseif strcmp(cfg_plot.plotstyle,'F-TxC')
                if nFreq <= 5
                  nRow = 1;
                  nCol = nFreq;
                elseif nFreq > 8 && nFreq <= 16
                  nRow = 2;
                  nCol = ceil(nFreq / 2);
                elseif nFreq > 16
                  nRow = 3;
                  nCol = ceil(nFreq / 3);
                end
                
                %subplot(nRow,nCol,count);
                subplot(nRow,nCol,count);
                imagesc(stat_all(sub,chn,lat,frqband).stat.time,1:nChan,squeeze(stat_all(sub,chn,lat,frqband).stat.model1(:,frq,:)),clim);
                axis xy;
                set(gca,'YTick',1:nChan);
                set(gca,'YTickLabel',stat_all(sub,chn,lat,frqband).stat.label);
                title(sprintf('%s, %.1f Hz',exper.subjects{sub},stat_all(sub,chn,lat,frqband).stat.freq(frq)));
                xlabel('Time (s)');
                ylabel('Electrode number');
                %colorbar;
                
                %         % Model 2
                %         %count = count + 1;
                %         %subplot(nFreq,2,count);
                %         subplot(2,nFreq,count + nFreq);
                %         imagesc(stat_all(sub,chn,lat,frqband).stat.time,1:nChan,squeeze(stat_all(sub,chn,lat,frqband).stat.model2(:,frq,:)),clim);
                %         axis xy;
                %         set(gca,'YTick',1:nChan);
                %         set(gca,'YTickLabel',stat_all(sub,chn,lat,frqband).stat.label);
                %         title(sprintf('%s, Model 2, %.1f Hz',exper.subjects{sub},stat_all(sub,chn,lat,frqband).stat.freq(frq)));
                %         xlabel('Time (s)');
                %         ylabel('Electrode number');
                %         %colorbar;
                
              end
            end
          elseif strcmp(cfg_plot.plotstyle,'multi')
            figure;
            cfg_plot.showlabels = 'yes';
            cfg_plot.interactive = 'yes';
            ft_multiplotTFR(cfg_plot,stat_all(sub,chn,lat,frqband).stat);
          elseif strcmp(cfg_plot.plotstyle,'C-TxF')
            figure;
            count = 0;
            [chanNum,j] = find(stat_all(sub,chn,lat,frqband).stat.model1 ~= 0);
            chanNum = unique(chanNum);
            
            if length(chanNum) <= 8
              nRow = 1;
              nCol = length(chanNum);
            elseif length(chanNum) > 8 && length(chanNum) <= 16
              nRow = 2;
              nCol = ceil(length(chanNum) / 2);
            elseif length(chanNum) > 16 && length(chanNum) < 30
              nRow = 3;
              nCol = ceil(length(chanNum) / 3);
            elseif length(chanNum) > 30
              nRow = 4;
              nCol = ceil(length(chanNum) / 4);
            end
            
            for c = 1:length(chanNum)
              count = count + 1;
              
              subplot(nRow,nCol,count);
              imagesc(stat_all(sub,chn,lat,frqband).stat.time,stat_all(sub,chn,lat,frqband).stat.freq,squeeze(stat_all(sub,chn,lat,frqband).stat.model1(c,:,:)),clim);
              axis xy;
              %set(gca,'YTickLabel',stat_all(sub,chn,lat,frqband).stat.freq);
              title(sprintf('%s, %s',exper.subjects{sub},stat_all(sub,chn,lat,frqband).stat.label{c}));
              xlabel('Time (s)');
              ylabel('Frequency (Hz)');
              
            end
          end
          set(gcf,'Name',sprintf('%s %.1f-%.1f s, %.1f-%.1f Hz: acc=%.1f%%, p=%.3f',exper.subjects{sub},cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2),stat_all(sub,chn,lat,frqband).stat.freq(1),stat_all(sub,chn,lat,frqband).stat.freq(end),stat_all(sub,chn,lat,frqband).stat.performance*100,stat_all(sub,chn,lat,frqband).stat.statistic.binomial));
          pos = get(gcf, 'Position');
          set(gcf, 'Units', 'pixels', 'Position', [pos(1), pos(2), figsize(2), figsize(1)]);
        end
      end
    end
  end
end

% %% Time-Frequency: plot it
%
% cfg_plot             = [];
% cfg_plot.layout = 'GSN-HydroCel-129.sfp';
% cfg_plot.zparam      = 'model1';
% cfg_plot.comment     = '';
% cfg_plot.colorbar    = 'yes';
% ft_topoplotTFR(cfg_plot,stat_all(sub,chn,lat,frqband).stat);
