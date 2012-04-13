% sosi classification test

adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

% % redefine which subjects to load
% exper.subjects = {'SOSI005'};

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%ana.eventValues = {exper.eventValues};
ana.eventValues = {{'RHSC','RHSI','RCR'}};

[data_raw] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'raw');

%% prep

cfg_prep = [];
cfg_prep.subjs = {'SOSI005'};
cfg_prep.freqs = 1;
cfg_prep.time = [0.5 0.8];

chans73 = [3 4 5 6 7 9 10 11 12 13 15 16 18 19 20 22 23 24 27 28 29 30 31 34 35 36 37 40 41 42 46 47 51 52 53 54 55 59 60 61 62 66 67 71 72 76 77 78 79 80 84 85 86 87 91 92 93 97 98 102 103 104 105 106 109 110 111 112 116 117 118 123 124];

%% process the data - ERP

sub = find(ismember(exper.subjects,cfg_prep.subjs));
ses = 1;

cfg = [];
cfg.parameter = 'trial';
cfg.keeptrials = 'yes'; % classifiers operate on individual trials
%cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS'})});
cfg.channel = eval(sprintf('{%s}',sprintf(repmat('''E%d'' ',size(chans73)),(chans73))));
tRHSC = ft_timelockanalysis(cfg,data_raw.RHSC.sub(sub).ses(ses).data);
tRHSI = ft_timelockanalysis(cfg,data_raw.RHSI.sub(sub).ses(ses).data);

%% ERP classification

cfg = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.method = 'crossvalidate';
cfg.design = [ones(size(tRHSC.trial,1),1); 2*ones(size(tRHSI.trial,1),1)]';
cfg.latency = cfg_prep.time;
cfg.mva = {ft_mv_standardizer ft_mv_csp('numchan',length(chans73)) ft_mv_svm};
stat = ft_timelockstatistics(cfg,tRHSC,tRHSI);

fprintf('accuracy = %.2f%%, p = %.4f\n',stat.performance*100,stat.pvalue);

contin = stat.cv.performance('contingency');
fprintf('\t\tPredicted\n');
fprintf('\t\t%s\t%s\n',ana.eventValues{1}{1},ana.eventValues{1}{2});
fprintf('True\t%s\t%d\t%d\n',ana.eventValues{1}{1},contin(1,1),contin(1,2));
fprintf('\t%s\t%d\t%d\n',ana.eventValues{1}{2},contin(2,1),contin(2,2));

%% plot it

cfg             = [];
cfg.interplimits = 'electrodes'; 
cfg.interpolation = 'linear';
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.xlim        = cfg_prep.time;
cfg.parameter      = 'model1';
cfg.comments    = '';
cfg.colorbar    = 'yes';
ft_topoplotER(cfg,stat);

%% process the data - TF

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 4:1:8;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
cfg.channel = eval(sprintf('{%s}',sprintf(repmat('''E%d'' ',size(chans73)),(chans73))));
cfg.toi          = 0.5:0.1:0.8;
cfg.keeptrials   = 'yes'; % classifiers operate on individual trials
tfrRHSC           = ft_freqanalysis(cfg,data_raw.RHSC.sub(sub).ses(ses).data);
tfrRHSI          = ft_freqanalysis(cfg,data_raw.RHSI.sub(sub).ses(ses).data);

%% simple TF classification

cfg         = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.method  = 'crossvalidate';
cfg.design  = [ones(size(tfrRHSC.powspctrm,1),1); 2*ones(size(tfrRHSI.powspctrm,1),1)]';
cfg.mva = {ft_mv_standardizer ft_mv_glmnet('lambda',0.1)};
stat = ft_freqstatistics(cfg,tfrRHSC,tfrRHSI);

fprintf('accuracy = %.2f%%, p = %.4f\n',stat.performance*100,stat.pvalue);

contin = stat.cv.performance('contingency');
fprintf('\t\tPredicted\n');
fprintf('\t\t%s\t%s\n',ana.eventValues{1}{1},ana.eventValues{1}{2});
fprintf('True\t%s\t%d\t%d\n',ana.eventValues{1}{1},contin(1,1),contin(1,2));
fprintf('\t%s\t%d\t%d\n',ana.eventValues{1}{2},contin(2,1),contin(2,2));

%% plot it

cfg             = [];
cfg.interplimits = 'electrodes'; 
cfg.interpolation = 'linear';
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.zparam      = 'model1';
cfg.comment     = '';
cfg.colorbar    = 'yes';
ft_topoplotTFR(cfg,stat);


%% longer TF classification

ses=1;

accuracy = nan(length(cfg_prep.subjs),length(cfg_prep.freqs));
pval = nan(length(cfg_prep.subjs),length(cfg_prep.freqs));

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
    %cfg.channels = eval(sprintf('{%s}',sprintf(repmat('''E%d'' ',size(chans73)),(chans73))));
    %cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS'})});
    %cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS','RPS'})});
    cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS'})});
    cfg.latency = [0.5 0.8];
    
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
    contin = stat.cv.performance('contingency');
    fprintf('\t\tPredicted\n');
    fprintf('\t\t%s\t%s\n',ana.eventValues{1}{1},ana.eventValues{1}{2});
    fprintf('True\t%s\t%d\t%d\n',ana.eventValues{1}{1},contin(1,1),contin(1,2));
    fprintf('\t%s\t%d\t%d\n',ana.eventValues{1}{2},contin(2,1),contin(2,2));
    
  end
end

