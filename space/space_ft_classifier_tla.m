% FRCE500 classification

dataroot = '~/data/SPACE/EEG/Sessions/face_house_ratings/ftpp/-1000_1000/ft_data';
%dataroot = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/face_house_ratings/ftpp/-1000_1000/ft_data';
dataroot = fullfile(dataroot,'Face_Face_SA_Face_SU_Face_VA_Face_VU_House_House_SA_House_SU_House_VA_House_VU_eq0_art_ftManual_ftICA');

dataroot = fullfile(dataroot,'tla_-1000_1000');

% load in the dataset
adFile = fullfile(dataroot,'analysisDetails.mat');
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,false);

ana.eventValues = {{'Face','House'}};

% pre-defined ROIs in this function
ana = mm_ft_elecGroups(ana);

% redefine which subjects to load
exper.subjects = {'SPACE001'};

[data_tla] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'tla');

% equatetrials = true;
equatetrials = false;

addtrials = false; %double number of trials in condition that starts out most frequent

subtract_trials = false; % subtract auditory trials to simulate skewed classes

%% classification from the tutorial

cfg_ana = [];

cfg_ana.latencies = [0.0 0.2; 0.1 0.3; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
lat = 2;
typ = 1;

sub = 1;
ses = 1;

data1 = data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data;
data2 = data_tla.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data;

cfg = [];
cfg.parameter = 'trial';
cfg.layout = 'GSN-HydroCel-129.sfp';
% cfg.method = 'crossvalidate_mm';
% cfg.resample = true;
cfg.method = 'crossvalidate';

cfg.statistic = {'accuracy' 'binomial' 'contingency'};

cfg.latency = cfg_ana.latencies(lat,:);
cfg.channel = 'all';
cfg.avgoverchan = 'no';
cfg.avgovertime = 'yes';

% data1 = ft_selectdata_new(cfg,data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data);
% data2 = ft_selectdata_new(cfg,data_tla.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data);

% data1 = ft_selectdata_old(data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data,...
%   'param',cfg.parameter,...
%   'toilim',cfg.latency,...
%   'channel',cfg.channel,...
%   'avgoverchan',cfg.avgoverchan,...
%   'avgovertime',cfg.avgovertime);
% 
% data2 = ft_selectdata_old(data_tla.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data,...
%   'param',cfg.parameter,...
%   'toilim',cfg.latency,...
%   'channel',cfg.channel,...
%   'avgoverchan',cfg.avgoverchan,...
%   'avgovertime',cfg.avgovertime);


cfg.design = [ones(size(data1.(cfg.parameter),1),1); 2*ones(size(data2.(cfg.parameter),1),1)]';

%cfg.mva = {dml.standardizer dml.enet('family','binomial','alpha',0.2)};
%cfg.mva = {dml.standardizer dml.enet('family','binomial','alpha',1)};
%cfg.mva = {dml.crossvalidator('mva',{dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true)};
cfg.mva = dml.analysis({dml.enet('family','binomial','alpha',0.2)});
% ,'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true;

stat = ft_timelockstatistics(cfg,data1,data2);

% stat.mymodel = stat.model{1}.primal;
stat.mymodel = stat.model{1}.weights;
cfg = [];
cfg.parameter = 'mymodel';
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.xlim = cfg_ana.latencies(lat,:);
cfg.comments = '';
cfg.colorbar = 'yes';
cfg.interplimits= 'electrodes';
ft_topoplotER(cfg,stat);

%%

if equatetrials
  
  ana.equateTrials = true;
  
  ses=1;
  
  % initialize random number generator
  rng('shuffle','twister');
  
else
  ana.equateTrials = false;
end

%% ERP: classification

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

cfg_ana.freqs = [0 0];

%cfg_ana.latencies = [0 1.0];
% cfg_ana.latencies = [0 0.62];
%cfg_ana.latencies = [1.0 1.45];
%cfg_ana.latencies = [0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
%cfg_ana.latencies = [-0.2 0.0; 0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
cfg_ana.latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
% cfg_ana.latencies = [0.0 0.5; 0.5 1.0];

%cfg_ana.chanStr = {{'center74'}};
cfg_ana.chanStr = {{'center91'}};
%cfg_ana.chanStr = {{'all129'}};

cfg = [];
cfg.parameter = 'trial';
cfg.layout = 'GSN-HydroCel-129.sfp';
% cfg.method = 'crossvalidate_mm';
cfg.method = 'crossvalidate';
% cfg.nfolds = 10;
cfg.nfolds = 5;
%cfg.nfolds = Inf;
%cfg.nfolds = 0.8;

cfg.statistic = {'accuracy' 'binomial' 'contingency'};
cfg.resample = true;

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
%cfg.avgovertime = 'yes';

if length(cfg_ana.chanStr) == 1
  cfg.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.chanStr{1})}));
  cfg_ana.nChan = length(cfg.channel);
end

% cfg.mva = {dml.standardizer dml.graphnet('family','binomial','L1',0.1)};
% cfg.mva = {dml.standardizer dml.enet('family','binomial','alpha',1)};
% cfg.mva = {dml.standardizer dml.glmnet('family','binomial','alpha',1)};
% cfg_ana.method = 'L1regLog-av_sub-skew_resample';

%cfg.mva = dml.crossvalidator('mva',dml.analysis({dml.enet('family','binomial','alpha',0.2)}),'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true);
cfg.mva = dml.analysis({dml.enet('family','binomial','alpha',0.2)});
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
    
    [data_tla] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'tla');
    
    if equatetrials
      % find the lowest number of trials within each subject
      lowEvNum = Inf(length(exper.subjects),1);
      
      for i = 1:length(ana.eventValues{typ})
        if size(data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data.trial,1) < lowEvNum(1)
          lowEvNum(1) = size(data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data.trial,1);
        end
      end
      for i = 1:length(ana.eventValues{typ})
        % if we're working on the smaller event, get a random sub-selection
        if size(data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data.trial,1) > lowEvNum(1)
          fprintf('%s: Subselecting %d (of %d) trials for %s...\n',exper.subjects{1},lowEvNum(1),size(data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data.trial,1),ana.eventValues{typ}{i});
          evNums = randperm(size(data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data.trial,1));
          % old version
          data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data = ft_selectdata(data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data,'rpt',sort(evNums(1:lowEvNum(1))));
          % new version, rpt selection doesn't work as of ft-20120116
          %cfg = [];
          %cfg.rpt = sort(evNums(1:lowEvNum(1)));
          %data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data = ft_selectdata(cfg,data_tla.(ana.eventValues{typ}{i}).sub(1).ses(ses).data);
        end
      end % i
    end % if

    %data1 = data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data;
    %data2 = data_tla.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data;
    data1 = data_tla.(ana.eventValues{typ}{1}).sub(1).ses(ses).data;
    data2 = data_tla.(ana.eventValues{typ}{2}).sub(1).ses(ses).data;
    
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
          
          %data = ft_selectdata_old(data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data,data_tla.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data,...
          data = ft_selectdata_old(data1,data2,...
            'param',cfg.parameter,...
            'toilim',cfg.latency,...
            'channel',cfg.channel,...
            'avgoverchan',cfg.avgoverchan,...
            'avgovertime',cfg.avgovertime);
          
          %data = ft_selectdata_new(cfg,data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data,data_tla.(ana.eventValues{typ}{2}).sub(sub).ses(ses).data);
          %data = ft_selectdata_new(cfg,data1,data2);
          
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
          
          stat_all(sub,chn,lat,frq).stat = ft_timelockstatistics(cfg,data1,data2);
          
          % debug
          keyboard
          
          %m = dml.crossvalidator('mva',{dml.standardizer dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'verbose',true);
          m = dml.crossvalidator('mva',{dml.standardizer dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'resample',cfg.resample,'verbose',true);
          cfg_ana.method = 'enet-alpha02-resample';
          
          % % L1
          % m = dml.crossvalidator('mva',{dml.standardizer dml.enet('family','binomial','alpha',1)},'stat',{'accuracy','binomial','contingency'},'resample',cfg.resample,'verbose',true);
          % cfg_ana.method = 'enet-alpha1-resample';
          
          % L2
          %m = dml.crossvalidator('mva',{dml.standardizer dml.enet('family','binomial','alpha',0)},'stat',{'accuracy','binomial','contingency'},'resample',cfg.resample,'verbose',true);
          %m = dml.crossvalidator('mva',{dml.enet('family','binomial','alpha',0)},'stat',{'accuracy','binomial','contingency'},'resample',cfg.resample,'verbose',true);
          %cfg_ana.method = 'enet-alpha0-resample';
          
          m = m.train(X,Y);
          %m.statistic
          %stat = m.test(X);
          
          % debug
          %keyboard
          
%           % gridsearch SVM for best C
%           m = dml.gridsearch('validator',dml.crossvalidator('type','split','stat','accuracy','mva',{dml.standardizer dml.svm('restart',false)}),'vars','C','vals',fliplr(logspace(-4,1,6)),'verbose',true);
%           m = m.train(X,Y);
%           stat = m.test(X);
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
%           stat = m.test(X);
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
    
    freqAvgStr = '';
    outfile = sprintf('ft_%s_%s_%dfreq%s_%s%s%d_%d%s',...
      vsStr,cfg_ana.method,size(cfg_ana.freqs,1),freqAvgStr,...
      chanStr,chanAvgStr,cfg_ana.latencies(1,1)*1000,cfg_ana.latencies(end,end)*1000,timeAvgStr);
    
    fprintf('Saving results to %s...\n',outfile);
    fid = fopen(fullfile(dataroot,[outfile,'.txt']),'w+');
    
    fprintf(fid,'ROI%s\t(%d channels)\n',chanStrTab,cfg_ana.nChan);
    fprintf(fid,'Events%s\n',vsStrTab);
    
    % print details and performance
    fprintf(fid,'%s\n',cfg_ana.method);
    
    fprintf(fid,'all times\t%.2f-%.2f s\n',data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.time(1),data_tla.(ana.eventValues{typ}{1}).sub(sub).ses(ses).data.time(end));
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
