%dataroot = fullfile(getenv('HOME'),'Downloads','FRCE_data_classification');

% % dpss 4:128 log spaced
% dataroot = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_dpss_pow_-500_1500_4_128/';
% hanning 4:1:64
%dataroot = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_64/';
dataroot = '/Users/matt/data/FRCE/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_64';
% % erp
%dataroot = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/tla_-1250_2250/';

%center73 = [3 4 5 6 7 9 10 11 12 13 15 16 18 19 20 22 23 24 27 28 29 30 31 34 35 36 37 40 41 42 46 47 51 52 53 54 55 59 60 61 62 66 67 71 72 76 77 78 79 80 84 85 86 87 91 92 93 97 98 102 103 104 105 106 109 110 111 112 116 117 118 123 124];

% load in the dataset
%adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_40/analysisDetails.mat';
%adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/Aud_AudForg_AudReca_Forg_Reca_Vis_VisForg_VisReca_eq0_art_nsAuto/tla_-1250_2250/analysisDetails.mat';
adFile = fullfile(dataroot,'analysisDetails.mat');
%[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,{'/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250','/Users/matt/data/FRCE'});
%ana.eventValues = {exper.eventValues};
%ana.eventValues = {{'AudForg','AudReca'},{'VisForg','VisReca'},{'Forg','Reca'},{'Aud','Vis'}};
ana.eventValues = {{'AudForg','AudReca'}};
%ana.eventValues = {{'VisForg','VisReca'}};
%ana.eventValues = {{'Forg','Reca'}};
%ana.eventValues = {{'Aud','Vis'}};

% pre-defined ROIs in this function
ana = mm_ft_elecGroups(ana);

% redefine which subjects to load
exper.subjects = {'FRCE 05'};

[data_freq] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

%[data_raw] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'raw');

equateTrials = true;

%% get the data so MVPA can use it

cfg_prep = [];

%   {'FRCE 03';
%   'FRCE 05';
%   'FRCE 06';
%   'FRCE 07';
%   'FRCE 08';
%   'FRCE 09'; % Grit said: good classifier performance
%   'FRCE 12';
%   'FRCE 13';
%   'FRCE 14'; % Grit said: good classifier performance
%   'FRCE 15';
%   };

cfg_prep.parameter = 'powspctrm';
cfg_prep.freqs = [4 8; 8 14; 14 28; 28 42; 42 64];
cfg_prep.times = [0 1.0];
cfg_prep.chanStr = {'center91'};
cfg_prep.chans = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_prep.chanStr)}));
ses = 1;

data = struct();

for frqband = 1:size(cfg_prep.freqs,1)
  fprintf('%.2f-%.2f Hz...\n',cfg_prep.freqs(frqband,1),cfg_prep.freqs(frqband,2));
  
  for sub = 1:length(exper.subjects)
    data.subjects(sub).pattern = [];
    data.subjects(sub).regressors = [];
    data.subjects(sub).selectors = [];
    data.subjects(sub).mask = [];
    
    
    for typ = 1:length(ana.eventValues)
      
      % initialize for this contrast type
      pattern_all = [];
      runs_all = [];
      regs_all = [];
      
      for cnd = 1:length(ana.eventValues{typ})
        fprintf('%s: %s\n',exper.subjects{sub},ana.eventValues{typ}{cnd});
        
        % collect the data for each frequency
        data_tf = [];
        runs = [];
        chansel = ismember(data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.label,cfg_prep.chans);
        freqsel = find(data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.freq >= cfg_prep.freqs(frqband,1) & data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.freq <= cfg_prep.freqs(frqband,2));
        timesel = find(data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.time >= cfg_prep.times(1) & data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.time <= cfg_prep.times(2));
        
        runCounter = 1;
        trialCounter = 0;
        for trl = 1:size(data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.(cfg_prep.parameter),1)
          % create a fake runs selector
          trialCounter = trialCounter + 1;
          if mod(trialCounter,5) == 0
            runCounter = runCounter + 1;
          end
          
          data_tf = cat(4,data_tf,data_freq.(ana.eventValues{typ}{cnd}).sub(sub).ses(ses).data.(cfg_prep.parameter)(trl,chansel,freqsel,timesel));
          runs = cat(2,runs,runCounter*ones(1,length(timesel)));
        end
        
        % create the pattern that MVPA needs
        pattern = permute(data_tf,[2 3 1 4]);
        %vDims = size(pattern);
        %pattern = reshape(pattern,prod(vDims(1:3)), vDims(4));
        pattern_all = cat(4,pattern_all,pattern);
        
        % regressors
        regs = zeros(length(ana.eventValues{typ}),length(runs));
        regs(cnd,:) = 1;
        regs_all = cat(2,regs_all,regs);
        
        % runs selector (TODO: NOT trial numbers, fix this)
        runs_all = cat(2,runs_all,runs);
        
      end % cnd
      
      % unjumble the fake runs
      [runs_all,i] = sort(runs_all);
      pattern_all = pattern_all(:,:,:,i);
      regs_all = regs_all(:,i);
      
      % initialize for this contrast type
      subj = init_subj('FRCE',exper.subjects{sub});
      
      % store the mask (channels x frequencies)
      vDims = size(pattern);
      nFeatures = prod(vDims(1:3));
      maskname = 'allFeatures';
      subj = init_object(subj,'mask',maskname);
      subj = set_mat(subj,'mask',maskname,ones([1 1 nFeatures]));
      
      % store the pattern
      patname = sprintf(repmat('%s_',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:});
      patname = patname(1:end-1);
      subj = load_matrix_pattern(subj,patname,maskname,pattern_all);
      
      % store the regressors (conditions x (timepoints x trials))
      subj = init_object(subj,'regressors','conds');
      subj = set_mat(subj,'regressors','conds',regs_all);
      subj = set_objfield(subj,'regressors','conds','condnames',ana.eventValues{typ});
      
      % store the runs selector
      subj = init_object(subj,'selector','runs');
      subj = set_mat(subj,'selector','runs',runs_all);
      
      % z-score the data
      subj = zscore_runs(subj,patname,'runs');
      
      % create the cross-validation indices
      subj = create_xvalid_indices(subj,'runs');
      
      % select some features
      subj = feature_select(subj,[patname,'_z'],'conds','runs_xval');
      
      %summarize(subj,'display_groups',false)
      %summarize(subj,'objtype','mask')
      
      class_args.train_funct_name = 'train_bp';
      class_args.test_funct_name = 'test_bp';
      class_args.nHidden = 0;
      
      [subj results] = cross_validation(subj,[patname,'_z'],'conds','runs_xval',[patname,'_z_thresh0.05'],class_args);
      
    end % typ
  end % sub
end % frqband

%% ERP: process the data

%subInd = find(ismember(exper.subjects,cfg_prep.subjs));
ses = 1;

if equateTrials
  % initialize random number generator
  rng('shuffle','twister');
  
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
%cfg_ana.latency = [-0.5 0];
%cfg_ana.latency = [0.5 1.0];
%cfg_ana.latency = [0.3 1.0];
%cfg_ana.latency = [.15 0.4];
cfg_ana.latency = [0 1.0];

%cfg_ana.chanStr = 'right';
%cfg_ana.chanStr = {'center74'};
cfg_ana.chanStr = {'center91'};

accuracy = nan(length(exper.subjects),1);
pval = nan(length(exper.subjects),1);
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
cfg.toi = (-0.3:0.02:1.0);
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
      if equateTrials
        ana.equateTrials = 'yes';
        % choose a random subset to equate trial numbers
        evNums = randperm(size(data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data.trial,2));
        cfg.trials = sort(evNums(1:lowEvNum(sub)));
      end
      data_freq.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data = ft_freqanalysis(cfg,data_raw.(ana.eventValues{typ}{i}).sub(sub).ses(ses).data);
    end
  end
end

%% Time-Frequency: if data is already processed but not equated
if equateTrials
  
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
end

%% Time-Frequency: Change in freq relative to baseline using absolute power

% save the original data
%data_freq_orig = data_freq;

cfg_fb = [];
cfg_fb.baseline = [-0.3 -0.1];
%cfg_fb.baselinetype = 'absolute';
cfg_fb.baselinetype = 'relative';

ana.blc = true;
ana.blc_method = cfg_fb.baselinetype;

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

ses = 1;

cfg_ana = [];

%cfg_ana.freqs = [4 8; 6 10; 8 12; 10 15; 12 19; 15 25; 19 30; 25 35; 30 40; 35 64; 64 100];
%cfg_ana.freqs = [4 8; 6 10; 8 12; 10 15; 12 19; 15 25; 19 30; 25 35; 30 40; 35 64];
%cfg_ana.freqs = [4 8; 8.1 12; 12.1 19; 19.1 30; 30.1 42; 42.1 64];
cfg_ana.freqs = [4 8; 8.1 14; 14.1 28; 28.1 42; 42.1 64];
%cfg_ana.freqs = [4 64];
%cfg_ana.freqs = [4 8; 8 14; 14 28; 28 64];
%cfg_ana.freqs = [4 8; 8.1 14; 14.1 32; 32.1 64; 64.1 100];

%cfg_ana.latencies = [0 1.0];
%cfg_ana.latencies = [0 0.5; 0.5 1.0; 1.0 1.45];
cfg_ana.latencies = [1.0 1.45];
%cfg_ana.latencies = [0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
%cfg_ana.latencies = [-0.2 0.0; 0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];

%cfg_ana.chanStr = {{'right'}};
cfg_ana.chanStr = {{'center74'}};
%cfg_ana.chanStr = {{'center91'}};
%cfg_ana.chanStr = {{'center74'},{'center91'},{'midline'},{'left'},{'right'},{'anterior'},{'posterior'}};
%cfg_ana.chanStr = {{'LAS'},{'RAS'},{'LPS'},{'RPS'},{'LAI'},{'RAI'},{'LPI'},{'RPI'},{'FI'},{'FS'},{'PS'},{'PI'}};
%cfg_ana.chanStr = {{'center74'},{'center91'},{'midline'},{'left'},{'right'},{'anterior'},{'posterior'},{'LAS'},{'RAS'},{'LPS'},{'RPS'},{'LAI'},{'RAI'},{'LPI'},{'RPI'},{'FI'},{'FS'},{'PS'},{'PI'}};
%cfg_ana.chanStr = {{'all129'}};

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.method = 'crossvalidate';
cfg.nfolds = 10;

%cfg.latency = [-0.2 0];
%cfg.latency = [0 1.0];
%cfg.latency = [0 0.6];
%cfg.latency = [0 0.3];
%cfg.latency = [0.3 0.6];
%cfg.latency = [0.5 1.0];

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
%cfg.avgoverfreq = 'yes';

% z-transform, feature selection (map data to a different space), SVM
%cfg_ana.method = 'cspFs10SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_filterer('maxfeatures',10) ft_mv_svm};

%cfg_ana.method = 'cspSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_svm};
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan,'numpatterns',3,'outputdatatype','logpowcsp','filttype','CSP0') ft_mv_svm};

%cfg_ana.method = 'filt10SVM';
%cfg.mva = {ft_mv_standardizer ft_mv_filterer('maxfeatures',10) ft_mv_svm('verbose',true)};

%cfg_ana.method = 'CSPgsCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',cfg_ana.nChan) ft_mv_gridsearch('verbose',true,'mva',ft_mv_svm,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars','C','vals',logspace(-3,3,7))};

%cfg_ana.method = 'gsFiltCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_gridsearch('verbose',true,'mva',{ft_mv_filterer ft_mv_svm},'validator',ft_mv_crossvalidator('nfolds',4,'metric','accuracy'),'mvidx',[1 2],'vars',{'maxfeatures' 'C'},'vals',{1:2:10 [1 10]})};
%cfg.mva = {ft_mv_standardizer ft_mv_gridsearch('verbose',true,'mva',{ft_mv_filterer ft_mv_svm},'validator',ft_mv_crossvalidator('nfolds',4,'metric','accuracy'),'mvidx',[1 2],'vars',{'maxfeatures' 'C'},'vals',{1:2:10 logspace(-3,3,7)})};

%cfg_ana.method = 'filt100gsCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_filterer('maxfeatures',100) ft_mv_gridsearch('verbose',true,'mva',ft_mv_svm,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars','C','vals',logspace(-3,3,7))};
%cfg_ana.method = 'gsCSVM';
%cfg.mva = {ft_mv_standardizer ft_mv_gridsearch('verbose',true,'mva',ft_mv_svm,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars','C','vals',logspace(-3,3,7))};

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
cfg_ana.method = 'L1regLogRLam01';
cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',1,'lambda',0.1,'family','binomial')};
%cfg_ana.method = 'L1regLogRLam01_ns';
%cfg.mva = {ft_mv_glmnet('alpha',1,'lambda',0.1,'family','binomial')};
% lam=1 doesn't perform very well
% lam=0.1 and 0.05 performs well
% lambda=0.2 doesn't perform well

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
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',1,'validator',ft_mv_crossvalidator('verbose',true,'nfolds',5,'metric','accuracy'),'family','binomial')};

% L2 optimize lambda using a 5-fold inner cross validation
%cfg_ana.method = 'L2regLogRLamOpt';
%cfg.mva = {ft_mv_standardizer ft_mv_glmnet('alpha',0,'validator',ft_mv_crossvalidator('verbose',true'nfolds',5,'metric','accuracy'),'family','binomial')};
% lumping all trials in one category, really bad performance


for typ = 1:length(ana.eventValues)
  vsStr = sprintf(repmat('%s',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:});
  vsStrTab = sprintf(repmat('\t%s',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:});
  if length(ana.eventValues{typ}) ~= 2
    error('Currently we can only classify 2 events, you included %d (%s)',length(ana.eventValues{typ}),vsStr);
  end
  
  % initialize some values
  accuracy = nan(length(exper.subjects),length(cfg_ana.chanStr),size(cfg_ana.latencies,1),size(cfg_ana.freqs,1));
  pval = nan(length(exper.subjects),length(cfg_ana.chanStr),size(cfg_ana.latencies,1),size(cfg_ana.freqs,1));
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
        
        for frqband = 1:size(cfg_ana.freqs,1)
          cfg.frequency = cfg_ana.freqs(frqband,:);
          
          cfg.design = [ones(size(data1.powspctrm,1),1); 2*ones(size(data2.powspctrm,1),1)]';
          
          % run it
          stat_all(sub,chn,lat,frqband).stat = ft_freqstatistics(cfg,data1,data2);
          
          % store some values
          accuracy(sub,chn,lat,frqband) = stat_all(sub,chn,lat,frqband).stat.performance;
          pval(sub,chn,lat,frqband) = stat_all(sub,chn,lat,frqband).stat.pvalue;
          continTab{sub,chn,lat,frqband} = stat_all(sub,chn,lat,frqband).stat.cv.performance('contingency');
          
          % print out some performance info
          fprintf('\n%s, %s, %s (%d chans), %.1f-%.1f s, %.1f-%.1f Hz\n',exper.subjects{sub},vsStr,chanStr(1:end-1),cfg_ana.nChan,cfg.latency(1),cfg.latency(2),cfg.frequency(1),cfg.frequency(2));
          fprintf('accuracy = %.2f%%, p = %.4f',stat_all(sub,chn,lat,frqband).stat.performance*100,stat_all(sub,chn,lat,frqband).stat.pvalue);
          if stat_all(sub,chn,lat,frqband).stat.pvalue < .05
            fprintf(' ***\n');
          else
            fprintf('\n');
          end
          
          % print out the contingency table
          fprintf('\t\tPredicted\n');
          fprintf('\t\t%s\t%s\n',ana.eventValues{typ}{1},ana.eventValues{typ}{2});
          fprintf('True\t%s\t%d\t%d\n',ana.eventValues{typ}{1},continTab{sub,chn,lat,frqband}(1,1),continTab{sub,chn,lat,frqband}(1,2));
          fprintf('\t%s\t%d\t%d\n',ana.eventValues{typ}{2},continTab{sub,chn,lat,frqband}(2,1),continTab{sub,chn,lat,frqband}(2,2));
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
    fprintf(fid,'p-values\n');
    fprintf(fid,'\t%s\n',sprintf(repmat('\t%.1f-%.1f Hz',1,size(cfg_ana.freqs,1)),cfg_ana.freqs'));
    for sub = 1:length(exper.subjects)
      fprintf(fid,'%s',exper.subjects{sub});
      for lat = 1:size(cfg_ana.latencies,1)
        % p-values for each subject, latency, and frequnecy band
        fprintf(fid,'\t%.1f-%.1fs%s\n',cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2),sprintf(repmat('\t%f',1,size(cfg_ana.freqs,1)),pval(sub,chn,lat,:)));
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % statistics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % t-test to see if accuracies are above 50%. within a frequency band,
    % within a time region, test accuracy vs 0.5.
    
    chance50 = 0.5*ones(size(exper.subjects));
    alpha = 0.05;
    
    fprintf(fid,'\nVersus chance\n');
    fprintf(fid,'%s\n',sprintf(repmat('\t%.1f-%.1f Hz\t',1,size(cfg_ana.freqs,1)),cfg_ana.freqs'));
    fprintf(fid,'%s\n',sprintf(repmat('\tp\tt(%d)',1,size(cfg_ana.freqs,1)),(size(exper.subjects,1)-1)*ones(size(cfg_ana.freqs,1),1)));
    
    for lat = 1:size(cfg_ana.latencies,1)
      fprintf(fid,'%.1f-%.1fs',cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2));
      for frqband = 1:size(cfg_ana.freqs,1)
        [h,p,ci,stats] = ttest(squeeze(accuracy(:,chn,lat,frqband)),chance50,alpha,'both');
        fprintf('%s, %.1f-%.1f s,\t%.1f-%.1f Hz,\tavg=%.2f%%:\t',chanStr(1:end-1),stat_all(1,chn,lat,frqband).stat.time(1),stat_all(1,chn,lat,frqband).stat.time(end),stat_all(1,chn,lat,frqband).stat.freq(1),stat_all(1,chn,lat,frqband).stat.freq(end),mean(accuracy(:,chn,lat,frqband),1)*100);
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
        if stat_all(sub,chn,lat,frqband).stat.pvalue < 0.05
          if strcmp(cfg_plot.plotstyle,'topo') || strcmp(cfg_plot.plotstyle,'F-TxC')
            figure;
            count = 0;
            nFreq = length(stat_all(sub,chn,lat,frqband).stat.freq);
            nChan = length(stat_all(sub,chn,lat,frqband).stat.label);
            
            for frqband = 1:nFreq
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
                
                cfg_plot.comment = sprintf('%s, %.1f Hz',exper.subjects{sub},stat_all(sub,chn,lat,frqband).stat.freq(frqband));
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
                imagesc(stat_all(sub,chn,lat,frqband).stat.time,1:nChan,squeeze(stat_all(sub,chn,lat,frqband).stat.model1(:,frqband,:)),clim);
                axis xy;
                set(gca,'YTick',1:nChan);
                set(gca,'YTickLabel',stat_all(sub,chn,lat,frqband).stat.label);
                title(sprintf('%s, %.1f Hz',exper.subjects{sub},stat_all(sub,chn,lat,frqband).stat.freq(frqband)));
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
          set(gcf,'Name',sprintf('%s %.1f-%.1f s, %.1f-%.1f Hz: acc=%.1f%%, p=%.3f',exper.subjects{sub},cfg_ana.latencies(lat,1),cfg_ana.latencies(lat,2),stat_all(sub,chn,lat,frqband).stat.freq(1),stat_all(sub,chn,lat,frqband).stat.freq(end),stat_all(sub,chn,lat,frqband).stat.performance*100,stat_all(sub,chn,lat,frqband).stat.pvalue));
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
