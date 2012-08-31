% FRCE500 classification

dataroot = '/Volumes/curranlab/Data/FRCE500/2 Session Recall/EEG/nspp/-1000_1250/ft_data';
dataroot = fullfile(dataroot,'no_recall_recall_eq0_art_nsAuto');

% wavelet, 4:8
%dataroot = fullfile(dataroot,'pow_wavelet_w5_pow_-500_980_4_8');

% fourier, 4:100
dataroot = fullfile(dataroot,'pow_wavelet_w5_fourier_-500_980_4_100');

% load in the dataset
adFile = fullfile(dataroot,'analysisDetails.mat');
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

ana.eventValues = {{'no_recall','recall'}};

% pre-defined ROIs in this function
ana = mm_ft_elecGroups(ana);

% % redefine which subjects to load
% exper.subjects = {'FRCE500 12'};

%[data_raw] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'raw');
%[data_freq] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

%equatetrials = true;
equatetrials = false;

addtrials = false; %double number of trials in condition that starts out most frequent

subtract_trials = false; % subtract auditory trials to simulate skewed classes

cfg_ld = [];
%cfg_ld.keeptrials = 'no';
cfg_ld.keeptrials = 'yes';
cfg_ld.equatetrials = 'no';
%cfg_ld.equatetrials = 'yes';

% type of input (used in the filename to load)
cfg_ld.ftype = 'fourier';
% cfg_ld.ftype = 'pow';

% type of output: 'pow', 'coh', 'phase'
cfg_ld.output = 'pow';

% % transformation: 'log10', 'log', 'vec', 'dB'
cfg_ld.transform = 'log10';
% %cfg_ld.transform = 'vec';

% normalization of single or average trials
cfg_ld.norm_trials = 'single'; % Grandchamp & Delorme (2011)
%cfg_ld.norm_trials = 'average';

% baseline type
% % 'zscore', 'absolute', 'relchange', 'relative', 'condition' (use ft_freqcomparison)
cfg_ld.baseline_type = 'zscore';
%cfg_ld.baseline_type = 'dB'; % relative baseline
% cfg_ld.baseline_type = 'absolute'; % subtract baseline mean
% cfg_ld.baseline_type = 'relchange'; % subtract and divide by baseline mean
% cfg_ld.baseline_type = 'relative'; % divide by baseline mean

% baseline period
cfg_ld.baseline_time = [-0.4 -0.1];
%cfg_ld.baseline_time = [-0.2 0];

% at what data stage should it be baseline corrected?
cfg_ld.baseline_data = 'pow';

%cfg_ld.saveFile = true;
cfg_ld.saveFile = false;

% only keep induced data by removing evoked?
cfg_ld.rmevoked = 'no';

if isfield(cfg_ld,'equatetrials') && strcmp(cfg_ld.equatetrials,'yes')
  eq_str = '_eq';
else
  eq_str = '';
end
if isfield(cfg_ld,'keeptrials') && strcmp(cfg_ld.keeptrials,'yes')
  kt_str = '_trials';
else
  kt_str = '_avg';
end
if isfield(cfg_ld,'rmevoked') && strcmp(cfg_ld.rmevoked,'yes')
  indu_str = '_induced';
else
  indu_str = '_whole';
end
saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s%s%s%s.mat',cfg_ld.output,eq_str,kt_str,indu_str));

% if exist(saveFile,'file')
%   fprintf('Loading saved file: %s\n',saveFile);
%   load(saveFile);
% else
%   fprintf('Running mm_ft_loadData\n');
%   if exist('data_evoked','var')
%     [data_freq,exper] = mm_ft_loadData(cfg_ld,exper,dirs,ana,data_evoked);
%   else
%     [data_freq,exper] = mm_ft_loadData(cfg_ld,exper,dirs,ana);
%   end
%   if cfg_ld.saveFile
%     fprintf('Saving %s...\n',saveFile);
%     save(saveFile,sprintf('data_%s',cfg_ld.output),'exper','cfg_ld');
%   end
% end
fprintf('Done.\n');


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
cfg_ana.freqs = [4 8; 8.1 12; 12.1 19; 19.1 28; 28.1 42; 42.1 64; 64.1 100];
%cfg_ana.freqs = [4 8; 6 12; 8.1 14; 12 18; 14.1 22; 18.1 28; 24.1 36; 28.1 42; 36.1 50; 42.1 64];
%cfg_ana.freqs = [4 8];

%cfg_ana.latencies = [0 1.0];
% cfg_ana.latencies = [0 0.62];
%cfg_ana.latencies = [1.0 1.45];
%cfg_ana.latencies = [0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
%cfg_ana.latencies = [-0.2 0.0; 0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
cfg_ana.latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6];

%cfg_ana.chanStr = {{'center74'}};
cfg_ana.chanStr = {{'center91'}};
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

cfg.avgoverchan = 'no';
cfg.avgovertime = 'yes';
%cfg.avgoverfreq = 'no';
cfg.avgoverfreq = 'yes';

if length(cfg_ana.chanStr) == 1
  cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr{1})}));
  cfg_ana.nChan = length(cfg.channel);
end

% cfg.mva = {dml.standardizer dml.graphnet('family','binomial','L1',0.1)};
% cfg.mva = {dml.standardizer dml.enet('family','binomial','alpha',1)};
% cfg.mva = {dml.standardizer dml.glmnet('family','binomial','alpha',1)};
% cfg_ana.method = 'L1regLog-av_sub-skew_resample';
cfg.mva = {dml.crossvalidator('mva',{dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true)};
cfg_ana.method = 'enet-alpha02-resample';

exper.allsubjects = exper.subjects;

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
  
  for sub = 1:length(exper.allsubjects)
    
    % so we only load one at a time
    exper.subjects = exper.allsubjects(sub);
    
    fprintf('%s\n',exper.allsubjects{sub});
    
    [data_freq,exper] = mm_ft_loadData(cfg_ld,exper,dirs,ana);
    
    %data1 = data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data;
    %data2 = data_freq.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data;
    data1 = data_freq.(ana.eventValues{typ}{1}).sub(1).ses(ses).data;
    data2 = data_freq.(ana.eventValues{typ}{2}).sub(1).ses(ses).data;
    
    for chn = 1:length(cfg_ana.chanStr)
      cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr{chn})}));
      cfg_ana.nChan = length(cfg.channel);
      chanStr = sprintf(repmat('%s_',1,length(cfg_ana.chanStr{chn})),cfg_ana.chanStr{chn}{:});
      chanStrTab = sprintf(repmat('\t%s',1,length(cfg_ana.chanStr{chn})),cfg_ana.chanStr{chn}{:});
      if isempty(cfg.channel)
        error('No channels selected for%s!',chanStrTab);
      end
      
      for lat = 1:size(cfg_ana.latencies,1)
        cfg.latency = cfg_ana.latencies(lat,:);
        
        for frq = 1:size(cfg_ana.freqs,1)
          cfg.frequency = cfg_ana.freqs(frq,:);
          
          %data = ft_selectdata_old(data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data,data_freq.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data,...
          data = ft_selectdata_old(data1,data2,...
            'param',cfg.parameter,...
            'foilim',cfg_ana.freqs(frq,:),...
            'toilim',cfg_ana.latencies(lat,:),...
            'channel',cfg.channel,...
            'avgoverchan',cfg.avgoverchan,...
            'avgoverfreq',cfg.avgoverfreq,...
            'avgovertime',cfg.avgovertime);
          
          %data = ft_selectdata_new(cfg,data_freq.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data,data_freq.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data);
          %data = ft_selectdata_new(cfg,data1,data2);
          
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
          %keyboard
          
          X = reshape(data.(cfg.parameter),reshapevec);
          Y = cfg.design';
          
%           m = dml.standardizer;
%           m = m.train(X);
%           Z = m.test(X);
%           
%           % search through L1 parameters
%           v = dml.graphnet.lambdapath(Z,Y,'binomial',50,1e-2);
%           m = dml.gridsearch('validator',dml.crossvalidator('type','split','stat','accuracy','mva',{dml.standardizer dml.graphnet('family','binomial','restart',false)}),'vars','L1','vals',v,'verbose',true);
%           tic;
%           m = m.train(X,Y);
%           toc
%           %m.statistic
%           figure;
%           plot(m.configs,m.outcome); xlabel('L1'); ylabel('accuracy');
%           Z2 = m.test(X);
          
          %m = dml.crossvalidator('mva',{dml.standardizer dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'verbose',true);
          m = dml.crossvalidator('mva',{dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true);
          m = m.train(X,Y);
          %m.statistic
          %Z2 = m.test(X);
          
%           % gridsearch SVM for best C
%           m = dml.gridsearch('validator',dml.crossvalidator('type','split','stat','accuracy','mva',{dml.standardizer dml.svm('restart',false)}),'vars','C','vals',fliplr(logspace(-4,1,6)),'verbose',true);
%           m = m.train(X,Y);
%           Z2 = m.test(X);
%           %m.statistic % doesn't work for gridsearch, not sure how to evaluate
%           
%           % Naive Bayes
%           m = dml.crossvalidator('mva',{dml.standardizer dml.naive},'stat',{'accuracy','binomial','contingency'},'verbose',true);
%           m = m.train(X,Y);
%           m.statistic
%           
%           % simple elastic net
%           m = dml.graphnet('family','binomial','L1',0.1);
%           m = m.train(X,Y);
%           %m.statistic
%           Z2 = m.test(X);
%           
%           % faster elastic net; no params: finds best lambda
%           m = dml.enet;
%           tic;
%           m = m.train(X,Y);
%           toc;
%           figure;
%           bar(m.model.weights);
%           figure;
%           imagesc(data.time,1:size(data.label,1),squeeze(reshape(m.model.weights,[size(data.(cfg.parameter),2), size(data.(cfg.parameter),3), size(data.(cfg.parameter),4)])));
%           axis xy
%           
%           % resample Naive Bayes
%           m = dml.crossvalidator('mva',{dml.standardizer dml.naive},'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true);
%           m = m.train(X,Y);
%           m.statistic
%           
%           % same as fieldtrip defaults
%           m = dml.crossvalidator('mva',dml.analysis({dml.standardizer('verbose',true) dml.svm('verbose',true)}),...
%             'stat',{'accuracy','binomial','contingency'},...
%             'type','nfold','folds',5);
%           m = m.train(X,Y);
%           m.statistic
%           
%           % run it
%           stat_all(sub,chn,lat,frq).stat = ft_freqstatistics(cfg,data1,data2);
          
          % store some values
          %accuracy(sub,chn,lat,frq) = stat_all(sub,chn,lat,frq).stat.statistic.accuracy;
          %binomial(sub,chn,lat,frq) = stat_all(sub,chn,lat,frq).stat.statistic.binomial;
          %continTab{sub,chn,lat,frq} = stat_all(sub,chn,lat,frq).stat.statistic.contingency;
          accuracy(sub,chn,lat,frq) = m.statistic('accuracy');
          binomial(sub,chn,lat,frq) = m.statistic('binomial');
          continTab{sub,chn,lat,frq} = m.statistic('contingency');
          
          % print out some performance info
          fprintf('\n%s, %s, %s (%d chans), %.1f-%.1f s, %.1f-%.1f Hz\n',exper.allsubjects{sub},vsStr,chanStr(1:end-1),cfg_ana.nChan,cfg.latency(1),cfg.latency(2),cfg.frequency(1),cfg.frequency(2));
          fprintf('accuracy = %.2f%%',accuracy(sub,chn,lat,frq)*100);
          fprintf(', binomial = %.4f',binomial(sub,chn,lat,frq));
          if binomial(sub,chn,lat,frq) < .05
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
    end % chn
  end % sub
  
  % put the subjects back
  exper.subjects = exper.allsubjects;
  sub=1;
  
  for chn = 1:length(cfg_ana.chanStr)
    cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr{chn})}));
    cfg_ana.nChan = length(cfg.channel);
    chanStr = sprintf(repmat('%s_',1,length(cfg_ana.chanStr{chn})),cfg_ana.chanStr{chn}{:});
    chanStrTab = sprintf(repmat('\t%s',1,length(cfg_ana.chanStr{chn})),cfg_ana.chanStr{chn}{:});
    if isempty(cfg.channel)
      error('No channels selected for%s!',chanStrTab);
    end
    
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
    %save(fullfile(dataroot,outfile),'stat_all');
    
    fprintf('Done\n\n');
    
  end % channel ROI
end % typ
