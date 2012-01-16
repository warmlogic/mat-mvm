%dataroot = fullfile(getenv('HOME'),'Downloads','FRCE_data_classification');

%dataroot = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_64/';

dataroot = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/tla_-1250_2250/';

%center73 = [3 4 5 6 7 9 10 11 12 13 15 16 18 19 20 22 23 24 27 28 29 30 31 34 35 36 37 40 41 42 46 47 51 52 53 54 55 59 60 61 62 66 67 71 72 76 77 78 79 80 84 85 86 87 91 92 93 97 98 102 103 104 105 106 109 110 111 112 116 117 118 123 124];

% load in the dataset
%adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_40/analysisDetails.mat';
%adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/tla_-1250_2250/analysisDetails.mat';
adFile = fullfile(dataroot,'analysisDetails.mat');
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,1);
%ana.eventValues = {exper.eventValues};
%ana.eventValues = {{'AudForg','AudReca'},{'VisForg','VisReca'},{'Forg','Reca'},{'Aud','Vis'}};
%ana.eventValues = {{'AudForg','AudReca'}};
%ana.eventValues = {{'VisForg','VisReca'}};
ana.eventValues = {{'Forg','Reca'}};
%ana.eventValues = {{'Aud','Vis'}};

% pre-defined in this function
ana = mm_ft_channelgroups(ana);

% redefine which subjects to load
%exper.subjects = {'FRCE 05'};

%[data_freq] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

[data_raw] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'raw');

equateTrials = true;

%% ERP: process the data

%subInd = find(ismember(exper.subjects,cfg_prep.subjs));
ses = 1;

if equateTrials
  % initialize random number generator
  %rng('shuffle');
  RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(clock*100)));
  
  % find the lowest number of trials within each subject
  lowEvNum = nan(length(exper.subjects),1);
  for sub = 1:length(exper.subjects)
    lowEvNum(sub) = Inf;
    for i = 1:length(ana.eventValues{1})
      if size(data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.trial,2) < lowEvNum(sub)
        lowEvNum(sub) = size(data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.trial,2);
      end
    end
  end
end

% initialize the struct
data_erp = struct;

% set the analysis parameters
cfg = [];
cfg.parameter = 'trial';
cfg.keeptrials = 'yes'; % classifiers operate on individual trials
%cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'FS'})});
cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'center73'})});

% get the ERP data
for sub = 1:length(exper.subjects)
  for i = 1:length(ana.eventValues{1})
    if equateTrials
      % choose a random subset to equate trial numbers
      evNums = randperm(size(data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.trial,2));
      cfg.trials = sort(evNums(1:lowEvNum(sub)));
    end
    data_erp.(ana.eventValues{1}{i}).sub(sub).ses(ses).data = ft_timelockanalysis(cfg,data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  end
end

%% ERP: classification

ses = 1;
cfg_ana = [];
cfg_ana.latency = [-0.5 0];
%cfg_ana.latency = [0.5 1.0];
%cfg_ana.latency = [0.3 1.0];
%cfg_ana.latency = [.15 0.4];
%cfg_ana.latency = [0 1.0];

%cfg_ana.chanStr = 'right';
cfg_ana.chanStr = {'center73'};

accuracy = nan(length(exper.subjects),1);
pval = nan(length(exper.subjects),1);
continTab = cell(length(exper.subjects),1);

cfg = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
%cfg.parameter = 'trial';
cfg.method = 'crossvalidate';
cfg.nfolds = 5;
cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)});


cfg.latency = cfg_ana.latency;
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg_ana.nChan = length(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)}));

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

for sub = 1:length(exper.subjects)
  fprintf('%s\n',exper.subjects{sub});
  
  data1 = data_erp.(ana.eventValues{1}{1}).sub(sub).ses(ses).data;
  data2 = data_erp.(ana.eventValues{1}{2}).sub(sub).ses(ses).data;
  
  cfg.design = [ones(size(data1.trial,1),1); 2*ones(size(data2.trial,1),1)]';
  
  stat = ft_timelockstatistics(cfg,data1,data2);
  
  accuracy(sub) = stat.performance;
  pval(sub) = stat.pvalue;
  continTab{sub} = stat.cv.performance('contingency');
  
  % print out some performance info
  fprintf('\n%s\naccuracy = %.2f%%, p = %.4f',exper.subjects{sub},accuracy(sub)*100,pval(sub));
  if pval(sub) < .05
    fprintf(' ***\n');
  else
    fprintf('\n');
  end
  
  fprintf('\t\tPredicted\n');
  fprintf('\t\t%s\t%s\n',ana.eventValues{1}{1},ana.eventValues{1}{2});
  fprintf('True\t%s\t%d\t%d\n',ana.eventValues{1}{1},continTab{sub}(1,1),continTab{sub}(1,2));
  fprintf('\t%s\t%d\t%d\n',ana.eventValues{1}{2},continTab{sub}(2,1),continTab{sub}(2,2));
  
end

%% ERP: save the results to a file

%vsStr = sprintf('%s%s',ana.eventValues{1}{1},sprintf(repmat('vs%s',1,1),ana.eventValues{1}{2}));
vsStr = sprintf(repmat('%s',1,length(ana.eventValues{1})),ana.eventValues{1}{:});
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
  fprintf(fid,'%s%s\n',exper.subjects{sub},sprintf(repmat('\t%f',1,size(pval,2)),pval(sub,:)));
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

if equateTrials
  % initialize random number generator
  %rng('shuffle');
  RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(clock*100)));
  
  % find the lowest number of trials within each subject
  lowEvNum = Inf(length(exper.subjects),1);
  for sub = 1:length(exper.subjects)
    for i = 1:length(ana.eventValues{1})
      if size(data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.trial,2) < lowEvNum(sub)
        lowEvNum(sub) = size(data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.trial,2);
      end
    end
  end
end

% initialize the struct
data_freq = struct;

cfg_ana = [];
%cfg_ana.chanStr = {'right'};
%cfg_ana.chanStr = {'LPS'};
cfg_ana.chanStr = {'center73'};

% set the analysis parameters
cfg = [];
%cfg.parameter = 'trial';
cfg.keeptrials = 'yes'; % classifiers operate on individual trials
cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)});

% % multitaper
% cfg.output       = 'pow';
% cfg.method       = 'mtmconvol';
% cfg.taper        = 'hanning';
% cfg.foi          = (4:1:64);
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
% %cfg.toi          = (0.3:0.02:1.0);
% %cfg.toi          = (0.1:0.02:0.6);
% cfg.toi          = (-0.3:0.02:1.0);

% wavelet
cfg.method = 'wavelet';
cfg.width = 6;
cfg.toi = (-0.3:0.02:1.0);
cfg.foi = (4:1:64);


% get the TF data
for sub = 1:length(exper.subjects)
  for i = 1:length(ana.eventValues{1})
    if equateTrials
      % choose a random subset to equate trial numbers
      evNums = randperm(size(data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.trial,2));
      cfg.trials = sort(evNums(1:lowEvNum(sub)));
    end
    data_freq.(ana.eventValues{1}{i}).sub(sub).ses(ses).data = ft_freqanalysis(cfg,data_raw.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  end
end

%% Time-Frequency: Change in freq relative to baseline using absolute power

cfg_fb = [];
cfg_fb.baseline = [-0.3 -0.1];
cfg_fb.baselinetype = 'absolute';

data_freq_orig = data_freq;

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        fprintf('%s, %s, %s, ',exper.subjects{sub},exper.sessions{ses},ana.eventValues{typ}{evVal});
        data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqbaseline(cfg_fb,data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
      end
    end
  end
end

%% Time-Frequency: classification

% % TF data is already processed
% if equateTrials
%   % initialize random number generator
%   %rng('shuffle');
%   RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(clock*100)));
%   
%   % find the lowest number of trials within each subject
%   lowEvNum = Inf(length(exper.subjects),1);
%   for sub = 1:length(exper.subjects)
%     for i = 1:length(ana.eventValues{1})
%       if size(data_freq.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.powspctrm,1) < lowEvNum(sub)
%         lowEvNum(sub) = size(data_freq.(ana.eventValues{1}{i}).sub(sub).ses(ses).data.powspctrm,1);
%       end
%     end
%   end
% end

ses = 1;

cfg_ana = [];
cfg_ana.freqs = [4 7;6 10;7 12;10 15;12 19;15 25; 19 30;25 35;30 40;35 64];
%cfg_ana.freqs = [4 8; 8 12; 12 28; 28 64];

accuracy = nan(length(exper.subjects),size(cfg_ana.freqs,1));
pval = nan(length(exper.subjects),size(cfg_ana.freqs,1));
continTab = cell(length(exper.subjects),size(cfg_ana.freqs,1));

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.method = 'crossvalidate';
cfg.nfolds = 5;

%cfg_ana.chanStr = 'right';
cfg_ana.chanStr = {'center73'};
cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)});

%cfg.latency = [-0.3 0];
%cfg.latency = [0 1.0];
cfg.latency = [0 0.6];
%cfg.latency = [0 0.3];
%cfg.latency = [0.3 0.6];
%cfg.latency = [0.5 1.0];

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
%cfg.avgoverfreq = 'yes';

cfg_ana.nChan = length(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr)}));

% z-transform, feature selection (map data to a different space), SVM
%cfg_ana.method = 'cspFs10SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_filterer('maxfeatures',10) ft_mv_svm};

cfg_ana.method = 'cspSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_svm};
cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan,'numpatterns',1,'outputdatatype','logpowcsp','filttype','CSP0') ft_mv_svm};

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
%cfg_ana.method = 'L2regLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0,'lambda',1,'family','binomial')};
% seems to have the best performance

% elastic net logistic regression
%cfg_ana.method = 'elasNetLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0.5,'lambda',1,'family','binomial')};
% terrible performance

% optimize lambda using a 5-fold inner cross validation
%cfg_ana.method = 'optLamLogR';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0.5, 'validator',ft_mv_crossvalidator('nfolds',5,'metric','accuracy'),'family','binomial')};

for sub = 1:length(exper.subjects)
  fprintf('%s\n',exper.subjects{sub});
  
  data1 = data_freq.(ana.eventValues{1}{1}).sub(sub).ses(ses).data;
  data2 = data_freq.(ana.eventValues{1}{2}).sub(sub).ses(ses).data;
  
  for frq = 1:size(cfg_ana.freqs,1)
    
    cfg.frequency = cfg_ana.freqs(frq,:);
    
    cfg.design = [ones(size(data1.powspctrm,1),1); 2*ones(size(data2.powspctrm,1),1)]';
    
    stat = ft_freqstatistics(cfg,data1,data2);
    
    accuracy(sub,frq) = stat.performance;
    pval(sub,frq) = stat.pvalue;
    continTab{sub,frq} = stat.cv.performance('contingency');
    
    % print out some performance info
    fprintf('\n%s, %d-%d Hz\naccuracy = %.2f%%, p = %.4f',exper.subjects{sub},cfg.frequency(1),cfg.frequency(2),accuracy(sub,frq)*100,pval(sub,frq));
    if pval(sub,frq) < .05
     fprintf(' ***\n');
    else
     fprintf('\n');
    end
    
    % print out the contingency table
    fprintf('\t\tPredicted\n');
    fprintf('\t\t%s\t%s\n',ana.eventValues{1}{1},ana.eventValues{1}{2});
    fprintf('True\t%s\t%d\t%d\n',ana.eventValues{1}{1},continTab{sub,frq}(1,1),continTab{sub,frq}(1,2));
    fprintf('\t%s\t%d\t%d\n',ana.eventValues{1}{2},continTab{sub,frq}(2,1),continTab{sub,frq}(2,2));
    
  end
end

%% Time-Frequency: save the results to a file

%vsStr = sprintf('%s%s',ana.eventValues{1}{1},sprintf(repmat('vs%s',1,1),ana.eventValues{1}{2}));
vsStr = sprintf(repmat('%s',1,length(ana.eventValues{1})),ana.eventValues{1}{:});
chanStr = sprintf(repmat('%s_',1,length(cfg_ana.chanStr)),cfg_ana.chanStr{:});
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
if strcmp(cfg.avgoverchan,'yes')
  chanAvgStr = 'Avg';
else
  chanAvgStr = '';
end
outfile = sprintf('ft_%s_%s_%dfold_%dsub_%dfreq%s_%s%s%d_%d%s.txt',...
  vsStr,cfg_ana.method,cfg.nfolds,length(exper.subjects),size(cfg_ana.freqs,1),freqAvgStr,...
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
fprintf(fid,'avgOverFreq\t%s\n',cfg.avgoverfreq);
fprintf(fid,'\n');

fprintf(fid,'Probabilities\n');
fprintf(fid,'%s',sprintf(repmat('\t%d-%d Hz',1,size(cfg_ana.freqs,1)),cfg_ana.freqs'));
fprintf(fid,'\tAverage\n');
for sub = 1:length(exper.subjects)
  fprintf(fid,'%s%s',exper.subjects{sub},sprintf(repmat('\t%.2f',1,size(accuracy,2)),accuracy(sub,:)));
  fprintf(fid,'\t%.2f\n',mean(accuracy(sub,:),2));
end
fprintf(fid,'Average%s',sprintf(repmat('\t%.2f',1,size(accuracy,2)),mean(accuracy,1)));
fprintf(fid,'\t%.2f\n',mean(mean(accuracy,1),2));
fprintf(fid,'\n');

% print p-values
fprintf(fid,'p-values\n');
for sub = 1:length(exper.subjects)
  fprintf(fid,'%s%s\n',exper.subjects{sub},sprintf(repmat('\t%f',1,size(pval,2)),pval(sub,:)));
end

fclose(fid);

%% Time-Frequency: plot it

cfg             = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.zparam      = 'model1';
cfg.comment     = '';
cfg.colorbar    = 'yes';
ft_topoplotTFR(cfg,stat);
