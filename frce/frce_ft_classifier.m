dataroot = fullfile(getenv('HOME'),'Downloads','FRCE_data_classification');

chans73 = [3 4 5 6 7 9 10 11 12 13 15 16 18 19 20 22 23 24 27 28 29 30 31 34 35 36 37 40 41 42 46 47 51 52 53 54 55 59 60 61 62 66 67 71 72 76 77 78 79 80 84 85 86 87 91 92 93 97 98 102 103 104 105 106 109 110 111 112 116 117 118 123 124];

% load in the dataset
adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_40/analysisDetails.mat';
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,1);
ana.eventValues = {exper.eventValues};
%ana.eventValues = {{'AudForg','AudReca'},{'VisForg','VisReca'},{'Forg','Reca'}};

% pre-defined in this function
ana = mm_ft_channelgroups(ana);

% redefine which subjects to load
%exper.subjects = {'FRCE 05'};

%[data_freq] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

[data_raw] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'raw');

equateTrials = 1;

%% process the data - ERP

%subInd = find(ismember(exper.subjects,cfg_prep.subjs));
ses = 1;

cfg = [];
cfg.parameter = 'trial';
cfg.keeptrials = 'yes'; % classifiers operate on individual trials
%cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS'})});
cfg.channel = eval(sprintf('{%s}',sprintf(repmat('''E%d'' ',size(chans73)),chans73)));

if equateTrials
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

%% ERP classification

ses = 1;
cfg_ana = [];
cfg_ana.latency = [0.3 1.0];

accuracy = nan(length(exper.subjects),1);
pval = nan(length(exper.subjects),1);
continTab = cell(length(exper.subjects),1);

cfg = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
%cfg.parameter = 'trial';
cfg.method = 'crossvalidate';
cfg.latency = cfg_ana.latency;
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',length(chans73)) ft_mv_svm};

for sub = 1:length(exper.subjects)
  
  VisForg = data_erp.(ana.eventValues{1}{1}).sub(sub).ses(ses).data;
  VisReca = data_erp.(ana.eventValues{1}{2}).sub(sub).ses(ses).data;
  
  cfg.design = [ones(size(VisForg.trial,1),1); 2*ones(size(VisReca.trial,1),1)];
  
  stat = ft_timelockstatistics(cfg,VisForg,VisReca);
  
  accuracy(sub) = stat.performance;
  pval(sub) = stat.pvalue;
  
  fprintf('accuracy = %.2f%%, p = %.4f\n',stat.performance*100,stat.pvalue);
  
  continTab{sub} = stat.cv.performance('contingency');
  fprintf('\t\tPredicted\n');
  fprintf('\t\t%s\t%s\n',ana.eventValues{1}{1},ana.eventValues{1}{2});
  fprintf('True\t%s\t%d\t%d\n',ana.eventValues{1}{1},continTab{sub}(1,1),continTab{sub}(1,2));
  fprintf('\t%s\t%d\t%d\n',ana.eventValues{1}{2},continTab{sub}(2,1),continTab{sub}(2,2));
  
end

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








%% Change in freq relative to baseline using absolute power

cfg_fb = [];
cfg_fb.baseline = [-0.3 -0.1];
cfg_fb.baselinetype = 'absolute';

%data_freq_orig = data_freq;

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

%% set config parameters for this analysis

cfg_prep = [];
cfg_prep.subjs = {
  'FRCE 03';
  'FRCE 05';
  'FRCE 06';
  'FRCE 07';
  'FRCE 08';
  'FRCE 09'; % Grit said: good classifier performance
  'FRCE 12';
  'FRCE 13';
  'FRCE 14'; % Grit said: good classifier performance
  'FRCE 15';
  };

cfg_prep.conds = {'VisForg','VisReca'}; % exper.eventValues;
cfg_prep.freqs = [4 7;6 10;7 12;10 15;12 19;15 25; 19 30;25 35;30 40];
cfg_prep.parameter = 'powspctrm';
cfg_prep.chans = chans73;
%cfg_prep.chans = (1:129);
cfg_prep.times = [0 1.0];
ses = 1;

%% ft classification

% try doing a cross validation

accuracy = nan(length(exper.subjects),length(cfg_prep.freqs));
pval = nan(length(exper.subjects),length(cfg_prep.freqs));
continTab = cell(length(exper.subjects),length(cfg_prep.freqs));

for sub = 1:length(cfg_prep.subjs)
  VisForg = data_freq.VisForg.sub(sub).ses(ses).data;
  VisReca = data_freq.VisReca.sub(sub).ses(ses).data;
  
  for frq = 1:size(cfg_prep.freqs,1)
    
    cfg = [];
    %cfg.layout = ft_prepare_layout([],ana);
    cfg.layout = 'GSN-HydroCel-129.sfp';
    cfg.method  = 'crossvalidate';
    cfg.nfolds = 5;
    cfg.design  = [ones(size(VisForg.powspctrm,1),1); 2*ones(size(VisReca.powspctrm,1),1)]';
    %cfg.channel = eval(sprintf('{%s}',sprintf(repmat('''E%d'' ',size(chans73)),(chans73))));
    cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS'})});
    %cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS','RPS'})});
    %cfg.latency = [0 1];
    %cfg.latency = [0 0.5];
    cfg.latency = [0.5 1];
    %cfg.frequency = [4 8];
    %cfg.frequency = [8 12];
    
    cfg.frequency = cfg_prep.freqs(frq,:);
    
    cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0.1)};
    
    stat = ft_freqstatistics(cfg,VisForg,VisReca);
    
    accuracy(sub,frq) = stat.performance;
    pval(sub,frq) = stat.pvalue;
    
    % print out some performance info
    fprintf('\n%s, %d-%d Hz, accuracy = %.2f%%, p = %.4f\n',cfg_prep.subjs{sub},cfg_prep.freqs(frq,1),cfg_prep.freqs(frq,2),accuracy(sub,frq)*100,pval(sub,frq));
    if pval(sub,frq) < .05
     fprintf(' ***\n\n\n');
    else
     fprintf('\n\n\n');
    end
    
    % print out the contingency table
    continTab{sub} = stat.cv.performance('contingency');
    fprintf('\t\tPredicted\n');
    fprintf('\t\t%s\t%s\n',ana.eventValues{1}{1},ana.eventValues{1}{2});
    fprintf('True\t%s\t%d\t%d\n',ana.eventValues{1}{1},continTab{sub}(1,1),continTab{sub}(1,2));
    fprintf('\t%s\t%d\t%d\n',ana.eventValues{1}{2},continTab{sub}(2,1),continTab{sub}(2,2));
    
  end
end

%% print out the results

fid = fopen(fullfile(dataroot,sprintf('fieldtrip_results_%s_%dfolds_%dsubjs_%dfreqs_%dchans.txt',cfg.method,cfg.nfolds,length(cfg_prep.subjs),length(cfg_prep.freqs),length(cfg.channels))),'w+');

% print performance
fprintf(fid,'%s',sprintf(repmat('\t%d-%d Hz',1,size(cfg_prep.freqs,1)),cfg_prep.freqs'));
fprintf(fid,'\tAverage\n');
for sub = 1:length(cfg_prep.subjs)
  fprintf(fid,'%s%s',cfg_prep.subjs{sub},sprintf(repmat('\t%.2f',1,size(accuracy,2)),accuracy(sub,:)));
  fprintf(fid,'\t%.2f\n',mean(accuracy(sub,:),2));
end
fprintf(fid,'Average%s',sprintf(repmat('\t%.2f',1,size(accuracy,2)),mean(accuracy,1)));
fprintf(fid,'\t%.2f\n',mean(mean(accuracy,1),2));

% print p-values
fprintf(fid,'p-values\n');
for sub = 1:length(cfg_prep.subjs)
  fprintf(fid,'%s%s\n',cfg_prep.subjs{sub},sprintf(repmat('\t%f',1,size(pval,2)),pval(sub,:)));
end

fclose(fid);

%% plot it
cfg             = [];
cfg.interplimits = 'electrodes'; 
%cfg.interpolation = 'linear';
%cfg.layout = ft_prepare_layout([],ana);
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.parameter      = 'model1';
cfg.comment     = '';
cfg.colorbar    = 'yes';
ft_topoplotTFR(cfg,stat);
