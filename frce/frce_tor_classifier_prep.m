% load in the dataset
adFile = '/Volumes/curranlab/Data/FRCE/EEG/Sessions/cueing paradigm/relabeled/eppp/-1250_2250/ft_data/VisForg_VisReca_eq0_art_nsAuto/pow_mtmconvol_hanning_pow_-500_1500_4_40/analysisDetails.mat';
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
ana.eventValues = {exper.eventValues};
ana = mm_ft_elecGroups(ana);
[data_freq] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'pow');

%[data_raw] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'raw');

chans73 = [3 4 5 6 7 9 10 11 12 13 15 16 18 19 20 22 23 24 27 28 29 30 31 34 35 36 37 40 41 42 46 47 51 52 53 54 55 59 60 61 62 66 67 71 72 76 77 78 79 80 84 85 86 87 91 92 93 97 98 102 103 104 105 106 109 110 111 112 116 117 118 123 124];

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

%% make a data struct (4D): subjects: trials (VisForg+VisReca) x channel x frequency x time

% this grabs all 129 channels

data = struct();

for sub = 1:length(cfg_prep.subjs)
  % initialize
  data.subjects(sub).data = [];
  data.subjects(sub).outcome = [];
  
  % collect the data concatenating acorss conditions
  for cnd = 1:length(cfg_prep.conds)
    fprintf('%s: %s: ',cfg_prep.subjs{sub},cfg_prep.conds{cnd});
    
    % collect the data for each frequency
    data_tf = [];
    for frq = 1:size(cfg_prep.freqs,1)
      fprintf('%d-%d Hz...',cfg_prep.freqs(frq,1),cfg_prep.freqs(frq,2));
      timesel = find(data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.time >= cfg_prep.times(1) & data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.time <= cfg_prep.times(2));
      freqsel = find(data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.freq >= cfg_prep.freqs(frq,1) & data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.freq <= cfg_prep.freqs(frq,2));
      data_tf = cat(3,data_tf,mean(data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.(cfg_prep.parameter)(:,:,freqsel,timesel),3));
    end
    
    % concatenate the data across conditions
    data.subjects(sub).data = cat(1,data.subjects(sub).data,data_tf);
    
    % concatenate the outcome vector across conditions
    if strcmp(cfg_prep.conds{cnd},'VisReca')
      data.subjects(sub).outcome = cat(1,data.subjects(sub).outcome,ones(size(data_tf,1),1));
    elseif strcmp(cfg_prep.conds{cnd},'VisForg')
      data.subjects(sub).outcome = cat(1,data.subjects(sub).outcome,-1*ones(size(data_tf,1),1));
    end
    fprintf('\n');
  end
end
fprintf('Done organizing %dD data struct.\n',ndims(data.subjects(1).data));

% % test
% sub=1;ses=1;
% trial = 5;chan = 12; frq = 2;
% figure
% plot(squeeze(data.subjects(sub).data(trial,chan,frq,:)));
% figure
% cnd = 1;
% timesel = find(data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.time >= cfg_prep.times(1) & data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.time <= cfg_prep.times(2));
% freqsel = find(data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.freq >= cfg_prep.freqs(frq,1) & data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.freq <= cfg_prep.freqs(frq,2));
% plot(squeeze(mean(data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.(cfg_prep.parameter)(trial,chan,freqsel,timesel),3)));

%% flatten the data struct (2D): subjects: trials (VisForg+VisReca) x (channel x frequency x time)

% this grabs only the 73 channels specified at the top

data_flat = struct();

ses = 1;
cnd = 1;

for sub = 1:length(cfg_prep.subjs)
  timesel = find(data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.time >= cfg_prep.times(1) & data_freq.(cfg_prep.conds{cnd}).sub(sub).ses(ses).data.time <= cfg_prep.times(2));
  
  % initialize
  data_flat.subjects(sub).data = [];
  data_flat.subjects(sub).colIndex = [];
  
  % get the data for each channel at each frequency and concatenate
  for chn = 1:length(cfg_prep.chans)
    for frq = 1:size(cfg_prep.freqs,1)
      % collect the data
      data_flat.subjects(sub).data = cat(2,data_flat.subjects(sub).data,squeeze(data.subjects(sub).data(:,cfg_prep.chans(chn),frq,:)));
      % construct an index for the columns: [channel; freq; timepoint]
      data_flat.subjects(sub).colIndex = cat(2,data_flat.subjects(sub).colIndex,[repmat(chn,1,length(timesel)); repmat(frq,1,length(timesel)); 1:length(timesel)]);
    end
  end
  
  data_flat.subjects(sub).outcome = data.subjects(sub).outcome;
end
fprintf('Done creating flat data struct.\n');

