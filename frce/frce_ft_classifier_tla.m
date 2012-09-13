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

