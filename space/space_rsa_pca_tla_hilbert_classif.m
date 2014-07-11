% space RSA

%% load the analysis details

subDir = '';
dataDir = fullfile('SPACE','EEG','Sessions','ftpp',subDir);
% Possible locations of the data files (dataroot)
serverDir = fullfile(filesep,'Volumes','curranlab','Data');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dreamDir = fullfile(filesep,'data','projects','curranlab');
localDir = fullfile(getenv('HOME'),'data');

% pick the right dataroot
if exist('serverDir','var') && exist(serverDir,'dir')
  dataroot = serverDir;
  %runLocally = 1;
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
  %runLocally = 1;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
  %runLocally = 0;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');

subjects = {
  %'SPACE001'; % low trial counts
  'SPACE002';
  'SPACE003';
  'SPACE004';
  'SPACE005';
  'SPACE006';
  'SPACE007';
  %'SPACE008'; % didn't perform task correctly, didn't perform well
  'SPACE009';
  'SPACE010';
  'SPACE011';
  'SPACE012';
  'SPACE013';
  'SPACE014';
  'SPACE015';
  'SPACE016';
  %'SPACE017'; % old assessment: really noisy EEG, half of ICA components rejected
  'SPACE018';
  %'SPACE019';
  'SPACE020';
  'SPACE021';
  'SPACE022';
  'SPACE027';
  'SPACE029';
  'SPACE037';
  %'SPACE039'; % noisy EEG; original EEG analyses stopped here
  'SPACE023';
  'SPACE024';
  'SPACE025';
  'SPACE026';
  'SPACE028';
  %'SPACE030'; % low trial counts
  'SPACE032';
  'SPACE034';
  'SPACE047';
  'SPACE049';
  'SPACE036';
  };

% only one cell, with all session names
sesNames = {'session_1'};

allowRecallSynonyms = true;

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% set up similarity analysis

% first set up classifier, if needed

% accurateClassifSelect = true;
accurateClassifSelect = false;

if accurateClassifSelect
  dataTypes_train = {'Face', 'House'};
  equateTrainTrials = true;
  standardizeTrain = true;
  alpha = 0.2;
  
  % do both P1 and P2 need to be classified correctly to use this trial?
  classifRequireP1 = true;
  classifRequireP2 = true;
  
  classif_str = 'classif';
else
  classif_str = 'noClassif';
end

% then set up similarity comparisons

% dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass', ...
%   'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'};

dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass'};

latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
  0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
  0 0.3; 0.3 0.6; 0.6 0.9; ...
  0 0.5; 0.5 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0];

% latencies = [0.0 0.5; 0.5 1.0];

% column numbers in trialinfo
% trialNumCol = 5;
phaseCountCol = 4;
stimNumCol = 6;
categNumCol = 7;
% pairNumCol = 13;

thisROI = {'LPI2','LPS','LT','RPI2','RPS','RT'};
% thisROI = {'center109'};
% thisROI = {'all129'};
% thisROI = {'LPI', 'PI', 'RPI'};
% thisROI = {'LPS'};
% thisROI = {'LPS', 'RPS'};
% thisROI = {'LAS', 'RAS'};
% thisROI = {'Fz'};
% thisROI = {'Cz'};
% thisROI = {'Pz'};
% thisROI = {'PI'};
% thisROI = {'posterior'};
% thisROI = {'LPS', 'RPS', 'LPI', 'PI', 'RPI'};
% thisROI = {'E70', 'E83'};
% thisROI = {'E83'};

if iscell(thisROI)
  roi_str = sprintf(repmat('%s',1,length(thisROI)),thisROI{:});
elseif ischar(thisROI)
  roi_str = thisROI;
end

cfg_sel = [];
cfg_sel.avgoverchan = 'no';
% cfg_sel.avgovertime = 'no';
cfg_sel.avgovertime = 'yes';

cfg_sel.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,thisROI)});

sim_method = 'cosine';
% sim_method = 'correlation';
% sim_method = 'spearman';

% eigenvalue options

% % keep components with eigenvalue >= 1
% eig_criterion = 'kaiser';

% % compute the percent explained variance expected from each component if
% % all events are uncorrelated with each other; keep it if above this level.
% % So, each component would explain 100/n, where n is the number of
% % events/components.
% eig_criterion = 'analytic';

% keep components that cumulatively explain at least 85% of the variance
eig_criterion = 'CV85';

%% train on expo images, test on multistudy images

sesNum = 1;

% ana.trl_order.multistudy_image = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'catNum', 'targ', 'spaced', 'lag', 'presNum', 'pairOrd', 'pairNum', 'cr_recog_acc', 'cr_recall_resp', 'cr_recall_spellCorr'};

if accurateClassifSelect
  ana.eventValues = {{'expo_stim','multistudy_image'}};
  ana.eventValuesSplit = { ...
    { ...
    {'Face','House'} ...
    { ...
    %'img_onePres' ...
    'img_RgH_rc_spac_p1','img_RgH_rc_spac_p2','img_RgH_rc_mass_p1','img_RgH_rc_mass_p2' ...
    'img_RgH_fo_spac_p1','img_RgH_fo_spac_p2','img_RgH_fo_mass_p1','img_RgH_fo_mass_p2' ...
    %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
    } ...
    } ...
    };
  
  if allowRecallSynonyms
    ana.trl_expr = { ...
      { ...
      { ...
      sprintf('eventNumber == %d & i_catNum == 1 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))) ...
      sprintf('eventNumber == %d & i_catNum == 2 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))) ...
      } ...
      { ...
      %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      } ...
      } ...
      };
  else
    ana.trl_expr = { ...
      { ...
      { ...
      sprintf('eventNumber == %d & i_catNum == 1 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))) ...
      sprintf('eventNumber == %d & i_catNum == 2 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))) ...
      } ...
      { ...
      %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      } ...
      } ...
      };
  end
else
  ana.eventValues = {{'multistudy_image'}};
  ana.eventValuesSplit = { ...
    { ...
    { ...
    %'img_onePres' ...
    'img_RgH_rc_spac_p1','img_RgH_rc_spac_p2','img_RgH_rc_mass_p1','img_RgH_rc_mass_p2' ...
    'img_RgH_fo_spac_p1','img_RgH_fo_spac_p2','img_RgH_fo_mass_p1','img_RgH_fo_mass_p2' ...
    %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
    } ...
    } ...
    };
  
  if allowRecallSynonyms
    ana.trl_expr = { ...
      { ...
      { ...
      %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      } ...
      } ...
      };
  else
    ana.trl_expr = { ...
      { ...
      { ...
      %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
      } ...
      } ...
      };
  end
end

%% load in the subject data

keeptrials = true;
[data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

%% decide who to kick out based on trial counts

% Subjects with bad behavior
% exper.badBehSub = {{}};
exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE030','SPACE039'}};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

%% set up the hilbert transform

parameter = 'powspctrm';
parameter_orig = 'trial';

% freqs = [2 4; 4 8; 8 12; 12 30; 30 80];
% freqs = [2 4; 4 8; 8 12; 12 30; 30 50];
freqs = [4 8; 8 12; 12 30; 30 50];
freq_str = sprintf('%dfreq%dto%d',size(freqs,1),freqs(1,1),freqs(end,end));

% bandpass filter the trials at a given frequency bin
sampleRate = ft_findcfg(data_tla.(exper.sesStr{1}).(sprintf('%s_p1',dataTypes{1})).sub(1).data.cfg,'fsample');
if isempty(sampleRate)
  sampleRate = 250;
end
bpFiltType = 'but';
bpFiltOrd = 4;

% compute the absolute value (complex vector length) of the Hilbert
% transform of the bandpass-filtered signal
hilbType = 'abs';

%% initialize to store results

similarity_all = cell(length(exper.subjects),length(exper.sesStr),length(dataTypes),size(latencies,1));
similarity_ntrials = nan(length(exper.subjects),length(exper.sesStr),length(dataTypes),size(latencies,1));

%% run similarity comparison

for sub = 1:length(exper.subjects)
  
  for ses = 1:length(exper.sesStr)
    sesStr = exper.sesStr{ses};
    
    if ~exper.badSub(sub,ses)
      fprintf('\t%s %s...\n',exper.subjects{sub},exper.sesStr{ses});
      
      if accurateClassifSelect
        % equate the training categories
        trlIndTrain = cell(length(dataTypes_train),1);
        if equateTrainTrials
          nTrainTrial = nan(length(dataTypes_train),1);
          for dt = 1:length(dataTypes_train)
            nTrainTrial(dt) = size(data_tla.(sesStr).(dataTypes_train{dt}).sub(sub).data.(parameter_orig),1);
          end
          fprintf('\tEquating training categories to have %d trials.\n',min(nTrainTrial));
          for dt = 1:length(dataTypes_train)
            trlInd = randperm(nTrainTrial(dt));
            trlIndTrain{dt,1} = sort(trlInd(1:min(nTrainTrial)));
          end
        else
          fprintf('\tNot equating training category trial counts.\n');
          for dt = 1:length(dataTypes_train)
            trlIndTrain{dt,1} = 'all';
          end
        end
      end
      
      % convert to frequency data using hilbert transform before processing
      if accurateClassifSelect
        for dt = 1:length(dataTypes_train)
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.freq = [];
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.dimord = 'rpt_chan_freq_time';
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.(parameter) = data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.(parameter_orig);
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data,parameter_orig);
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data,'avg');
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data,'var');
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data,'dof');
          
          allFreq = nan( ...
            size(data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.(parameter),1), ...
            size(data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.(parameter),2), ...
            size(freqs,1), ...
            size(data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.(parameter),3));
          
          for f = 1:size(freqs,1)
            data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.freq = cat(2,data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.freq,mean(freqs(f,:)));
            for tr = 1:size(data_tla.(sesStr).(dataTypes_train{dt}).sub(sub).data.(parameter),1)
              thisFreq = ft_preproc_bandpassfilter(squeeze(data_tla.(sesStr).(dataTypes_train{dt}).sub(sub).data.(parameter)(tr,:,:)),sampleRate,freqs(f,:),bpFiltOrd,bpFiltType);
              thisFreq = ft_preproc_hilbert(thisFreq,hilbType);
              allFreq(tr,:,f,:) = thisFreq;
            end
          end
          data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data.(parameter) = allFreq;
          clear allFreq
        end
      end
      
      presStr = {'p1','p2'};
      for dt = 1:length(dataTypes)
        for p = 1:length(presStr)
          dataType = sprintf('%s_%s',dataTypes{dt},presStr{p});
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.freq = [];
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.dimord = 'rpt_chan_freq_time';
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.(parameter) = data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.(parameter_orig);
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data,parameter_orig);
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data,'avg');
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data,'var');
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data = rmfield(data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data,'dof');
          allFreq = nan( ...
            size(data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.(parameter),1), ...
            size(data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.(parameter),2), ...
            size(freqs,1), ...
            size(data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.(parameter),3));
          for f = 1:size(freqs,1)
            data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.freq = cat(2,data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.freq,mean(freqs(f,:)));
            for tr = 1:size(data_tla.(sesStr).(dataType).sub(sub).data.(parameter),1)
              thisFreq = ft_preproc_bandpassfilter(squeeze(data_tla.(sesStr).(dataType).sub(sub).data.(parameter)(tr,:,:)),sampleRate,freqs(f,:),bpFiltOrd,bpFiltType);
              thisFreq = ft_preproc_hilbert(thisFreq,hilbType);
              allFreq(tr,:,f,:) = thisFreq;
            end
            %for tr = 1:length(p1_ind)
            %  thisFreq = ft_preproc_bandpassfilter(squeeze(data_tla.(sesStr).(dataType).sub(sub).data.(parameter)(p1_ind(tr),:,:)),sampleRate,freqs(f,:),bpFiltOrd,bpFiltType);
            %  thisFreq = ft_preproc_hilbert(thisFreq,hilbType);
            %  allFreq(p1_ind(tr),:,f,:) = thisFreq;
            %end
          end
          data_tla.(exper.sesStr{ses}).(dataType).sub(sub).data.(parameter) = allFreq;
          clear allFreq
        end
      end
      
      % train for a given latency and set of electrodes
      for lat = 1:size(latencies,1)
        cfg_sel.latency = latencies(lat,:);
        
        if accurateClassifSelect
          % select the training data
          data_train = struct;
          
          % select the data
          for dt = 1:length(dataTypes_train)
            cfg_sel.trials = trlIndTrain{dt};
            data_train.(dataTypes_train{dt}) = ft_selectdata_new(cfg_sel,data_tla.(exper.sesStr{ses}).(dataTypes_train{dt}).sub(sub).data);
          end
          
          % get the category number for each training image
          imageCategory_train = data_train.(dataTypes_train{1}).trialinfo(:,categNumCol);
          for dt = 2:length(dataTypes_train)
            imageCategory_train = cat(1,imageCategory_train,data_train.(dataTypes_train{dt}).trialinfo(:,categNumCol));
          end
          
          % concatenate the data
          dat_train = data_train.(dataTypes_train{1}).(parameter);
          for dt = 2:length(dataTypes_train)
            dat_train = cat(1,dat_train,data_train.(dataTypes_train{dt}).(parameter));
          end
          
          dim = size(dat_train);
          dat_train = reshape(dat_train, dim(1), prod(dim(2:end)));
          
          if standardizeTrain
            fprintf('\t\tStandardizing the training data...');
            
            m = dml.standardizer;
            m = m.train(dat_train);
            dat_train = m.test(dat_train);
            fprintf('Done.\n');
          end
          
          fprintf('\t\tTraining classifier...');
          %facehouse = {dml.standardizer dml.enet('family','binomial','alpha',alpha)};
          facehouse = dml.enet('family','binomial','alpha',alpha);
          facehouse = facehouse.train(dat_train,imageCategory_train);
          %facehouse_svm = dml.svm;
          %facehouse_svm = facehouse_svm.train(dat_train,imageCategory_train);
          fprintf('Done.\n');
        end
        
        data_p1_append_str = [];
        data_p2_append_str = [];
        dTypes_p1 = [];
        dTypes_p2 = [];
        
        for d = 1:length(dataTypes)
          data_p1_append_str = cat(2,data_p1_append_str,sprintf(',data_tla.%s.%s_p1.sub(%d).data',sesStr,dataTypes{d},sub));
          data_p2_append_str = cat(2,data_p2_append_str,sprintf(',data_tla.%s.%s_p2.sub(%d).data',sesStr,dataTypes{d},sub));
          
          dTypes_p1 = cat(2,dTypes_p1,d*ones(1,size(data_tla.(sesStr).(sprintf('%s_p1',dataTypes{d})).sub(sub).data.(parameter),1)));
          dTypes_p2 = cat(2,dTypes_p2,d*ones(1,size(data_tla.(sesStr).(sprintf('%s_p2',dataTypes{d})).sub(sub).data.(parameter),1)));
        end
        
        cfg_ad = [];
        cfg_ad.parameter = parameter;
        data_p1 = eval(sprintf('ft_appendfreq(cfg_ad%s);',data_p1_append_str));
        data_p2 = eval(sprintf('ft_appendfreq(cfg_ad%s);',data_p2_append_str));
        
        % find where p1 and p2 both exist for each trial
        p1_ind = [];
        p2_ind = [];
        if accurateClassifSelect
          imageCategory_test = []; % 1=face, 2=house
        end
        
        for p = 1:size(data_p1.trialinfo,1)
          p1_trlInd = p;
          p1_phaseCount = data_p1.trialinfo(p1_trlInd,phaseCountCol);
          p1_stimNum = data_p1.trialinfo(p1_trlInd,stimNumCol);
          p1_categNum = data_p1.trialinfo(p1_trlInd,categNumCol);
          
          p2_trlInd = find(...
            data_p2.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
            data_p2.trialinfo(:,stimNumCol) == p1_stimNum & ...
            data_p2.trialinfo(:,categNumCol) == p1_categNum);
          
          if ~isempty(p2_trlInd)
            if length(p2_trlInd) == 1
              p1_ind = cat(2,p1_ind,p1_trlInd);
              p2_ind = cat(2,p2_ind,p2_trlInd);
              if accurateClassifSelect
                imageCategory_test = cat(2,imageCategory_test,p1_categNum);
              end
            else
              warning('more than one p2 trial found');
              keyboard
            end
          end
        end
        
        % select the requested events, times, channels, etc.
        cfg_sel.trials = p1_ind;
        data_p1 = ft_selectdata_new(cfg_sel,data_p1);
        cfg_sel.trials = p2_ind;
        data_p2 = ft_selectdata_new(cfg_sel,data_p2);
        
        dat_p1 = data_p1.(parameter);
        dat_p2 = data_p2.(parameter);
        
        clear data_p1 data_p2
        
        dTypes_p1 = dTypes_p1(p1_ind);
        dTypes_p2 = dTypes_p2(p2_ind);
        
        if accurateClassifSelect
          % reshape first, test all trials at once
          %
          % unroll data for each trial in the second dimension
          dim_p1 = size(dat_p1);
          dim_p2 = size(dat_p2);
          dat_p1_p2_z = zscore(cat(1,reshape(dat_p1, dim_p1(1), prod(dim_p1(2:end))),reshape(dat_p2, dim_p2(1), prod(dim_p2(2:end)))));
          
          % attempt to classifiy study trials
          
          if classifRequireP1
            Z = facehouse.test(dat_p1_p2_z(1:length(p1_ind),:));
            probabilityClassP1 = Z;
            
            [Y,I] = max(Z,[],2);
            
            correctClassP1 = I == imageCategory_test';
          else
            correctClassP1 = true(size(p1_ind));
          end
          
          if classifRequireP2
            Z = facehouse.test(dat_p1_p2_z(length(p1_ind)+1:end,:));
            probabilityClassP2 = Z;
            
            [Y,I] = max(Z,[],2);
            
            correctClassP2 = I == imageCategory_test';
          else
            correctClassP2 = true(size(p2_ind));
          end
          
          % only compare these trials
          p1_ind = find(correctClassP1 & correctClassP2);
          p2_ind = find(correctClassP1 & correctClassP2);
          
          dat_p1 = dat_p1(p1_ind,:,:,:);
          dat_p2 = dat_p2(p2_ind,:,:,:);
          
          dTypes_p1 = dTypes_p1(p1_ind);
          dTypes_p2 = dTypes_p2(p2_ind);
        end
        
%         if accurateClassifSelect
%           
%           % attempt to classifiy study trials
%           probabilityClassP1 = nan(length(p1_ind),2);
%           correctClassP1 = true(size(p1_ind));
%           if classifRequireP1
%             % TODO: reshape first, test all trials at once
%             for p = 1:length(p1_ind)
%               dat1 = dat_p1(p,:,:,:);
%               dim = size(dat1);
%               dat1 = reshape(dat1, dim(1), prod(dim(2:end)));
%               
%               Z = facehouse.test(zscore(dat1));
%               probabilityClassP1(p,:) = Z;
%               
%               [Y,I] = max(Z,[],2);
%               
%               correctClassP1(p) = I == imageCategory_test(p);
%             end
%           end
%           
%           probabilityClassP2 = nan(length(p1_ind),2);
%           correctClassP2 = true(size(p2_ind));
%           if classifRequireP2
%             for p = 1:length(p2_ind)
%               dat2 = dat_p2(p,:,:,:);
%               dim = size(dat2);
%               dat2 = reshape(dat2, dim(1), prod(dim(2:end)));
%               
%               Z = facehouse.test(zscore(dat2));
%               probabilityClassP2(p,:) = Z;
%               
%               [Y,I] = max(Z,[],2);
%               
%               correctClassP2(p) = I == imageCategory_test(p);
%             end
%           end
%           
%           % only compare these trials
%           p1_ind = find(correctClassP1 & correctClassP2);
%           p2_ind = find(correctClassP1 & correctClassP2);
%           
%           dat_p1 = dat_p1(p1_ind,:,:,:);
%           dat_p2 = dat_p2(p2_ind,:,:,:);
%           
%           dTypes_p1 = dTypes_p1(p1_ind);
%           dTypes_p2 = dTypes_p2(p2_ind);
%         end
        
        if ~isempty(p1_ind) && ~isempty(p2_ind)
          % unroll data for each trial in the second dimension
          dim_p1 = size(dat_p1);
          dim_p2 = size(dat_p2);
          dat_p1_p2 = cat(1,reshape(dat_p1, dim_p1(1), prod(dim_p1(2:end))),reshape(dat_p2, dim_p2(1), prod(dim_p2(2:end))));
          
          clear dat_p1 dat_p2
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Compute similarity
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          % variables: columns = electrode x time (unrolled)
          % observations/instances: rows = events
          
          % apply PCA to data
          if exist('pca','file')
            [evec_p1_p2, data_pcaspace, eval_p1_p2] = pca(zscore(dat_p1_p2), 'Economy', true);
          else
            [evec_p1_p2, data_pcaspace, eval_p1_p2] = princomp(zscore(dat_p1_p2), 'econ');
          end
          
          if strcmp(eig_criterion,'kaiser')
            crit_eig = eval_p1_p2 >= 1;
          elseif strcmp(eig_criterion,'analytic')
            % analytic: keep PC if percent variance explained is above
            % 100 / number of variables
            
            % convert to percent variance explained
            eval_PVE = (eval_p1_p2 ./ sum(eval_p1_p2)) .* 100;
            crit_eig = eval_PVE > (100 / size(eval_PVE,1));
          elseif strcmp(eig_criterion,'CV85')
            % Cumulative variance 85%: keep PCs that explain at least 85%
            % of variance
            
            % convert to percent variance explained
            eval_PVE = (eval_p1_p2 ./ sum(eval_p1_p2)) .* 100;
            eval_CV = cumsum(eval_PVE);
            cutoff_eval = find(eval_CV>=85,1,'first');
            
            crit_eig = false(size(eval_p1_p2));
            crit_eig(1:cutoff_eval) = true;
          elseif strcmp(eig_criterion,'none')
            crit_eig = true(length(eval_p1_p2),1);
          end
          % remove features with eigenvalues that didn't pass criterion
          %
          % evec_p1_p2 (coeff) lets you map from PCA space to original
          % feature space
          evec_p1_p2_crit = evec_p1_p2(:, crit_eig);
          feature_vectors = data_pcaspace(:, crit_eig);
          
          %%%%%%
          % more feature selection done here? use dummy selection for now
          %%%%%%
          select_inds = true(1, size(feature_vectors, 2));
          
          evec_p1_p2_final = evec_p1_p2_crit(:, select_inds);
          feature_vectors = feature_vectors(:, select_inds);
          
          % normalize the vector lengths of each event
          feature_vectors = feature_vectors ./ repmat(sqrt(sum(feature_vectors.^2, 2)), 1, size(feature_vectors, 2));
          
          % compute the similarities between each pair of events
          similarities = 1 - squareform(pdist(feature_vectors, sim_method));
          
          % add it to the full set
          for d = 1:length(dataTypes)
            similarity_all{sub,ses,d,lat} = similarities(cat(2,dTypes_p1,dTypes_p2) == d,cat(2,dTypes_p1,dTypes_p2) == d);
            similarity_ntrials(sub,ses,d,lat) = sum(cat(2,dTypes_p1,dTypes_p2) == d);
            %similarity_all{sub,ses,d,lat} = similarities;
            %similarity_ntrials(sub,ses,d,lat) = length(p1_ind);
          end
          
        end % ~isempty
        
      end % lat
      
      %end % d
    end % ~badSub
  end % ses
end % sub

saveFile = fullfile(dirs.saveDirProc,sprintf('RSA_PCA_hilbert_%s_%s_%s_%s_%dlat_%s_%sAvgT_%s.mat',sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,cfg_sel.avgovertime,date));
fprintf('Saving: %s\n',saveFile);
save(saveFile,'exper','dataTypes','thisROI','cfg_sel','eig_criterion','sim_method','classif_str','freqs','latencies','similarity_all','similarity_ntrials');
fprintf('Done.\n');

%% load

% analysisDate = '13-Jun-2014';
analysisDate = '07-Jul-2014';

% thisROI = {'center109'};
thisROI = {'LPI2','LPS','LT','RPI2','RPS','RT'};
% thisROI = {'LPS','RPS'};
% thisROI = {'LT','RT'};

% analysisDate = '10-Jul-2014';
% thisROI = {'LPI2','RPI2'};

if iscell(thisROI)
  roi_str = sprintf(repmat('%s',1,length(thisROI)),thisROI{:});
elseif ischar(thisROI)
  roi_str = thisROI;
end

latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
  0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
  0 0.3; 0.3 0.6; 0.6 0.9; ...
  0 0.5; 0.5 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0];

origDataType = 'hilbert';

freqs = [4 8; 8 12; 12 30; 30 50];
% freqs = [4 8];
% freqs = [8 12];
% freqs = [12 30];
% freqs = [30 50];
freq_str = sprintf('%dfreq%dto%d',size(freqs,1),freqs(1,1),freqs(end,end));

avgovertime = 'yes';
% avgovertime = 'no';

sim_method = 'cosine';
% sim_method = 'correlation';
% sim_method = 'spearman';

% accurateClassifSelect = true;
accurateClassifSelect = false;
if accurateClassifSelect
  classif_str = 'classif';
else
  classif_str = 'noClassif';
end

% eig_criterion = 'CV85';
eig_criterion = 'kaiser';
% eig_criterion = 'analytic';

saveDirProc = '~/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
rsaFile = fullfile(saveDirProc,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%dlat_%s_%sAvgT_%s_cluster.mat',origDataType,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,avgovertime,analysisDate));
if exist(rsaFile,'file')
  fprintf('loading: %s...',rsaFile);
  load(rsaFile);
  fprintf('Done.\n');
else
  error('does not exist: %s',rsaFile);
end

%% stats

nTrialThresh = 8; % 31
% nTrialThresh = 14; % 30
% nTrialThresh = 16; % 28
% nTrialThresh = 17; % 24

plotit = false;

noNans = true(length(exper.subjects),length(exper.sesStr));

passTrlThresh = true(length(exper.subjects),length(exper.sesStr));

mean_similarity = struct;
for d = 1:length(dataTypes)
  mean_similarity.(dataTypes{d}) = nan(length(exper.subjects),length(exper.sesStr),size(latencies,1));
  for lat = 1:size(latencies,1)
    
    for sub = 1:length(exper.subjects)
      for ses = 1:length(exper.sesStr)
        
        % Average Pres1--Pres2 similarity
        if ~isempty(diag(similarity_all{sub,ses,d,lat}))
          if similarity_ntrials(sub,ses,d,lat) < nTrialThresh
            passTrlThresh(sub,ses) = false;
          end
          
          mean_similarity.(dataTypes{d})(sub,ses,lat) = mean(diag(similarity_all{sub,ses,d,lat},size(similarity_all{sub,ses,d,lat},1) / 2));
          %mean_similarity.(dataTypes{d}) = cat(1,mean_similarity.(dataTypes{d}),mean(diag(similarity_all{sub,ses,d,lat},size(similarity_all{sub,ses,d,lat},1) / 2)));
          
          if plotit
            figure
            imagesc(similarity_all{sub,ses,d,lat});
            colorbar;
            axis square;
            title(sprintf('%s, %.2f to %.2f',strrep(dataTypes{d},'_','-'),latencies(lat,1),latencies(lat,2)));
          end
        else
          noNans(sub,ses) = false;
        end
      end
    end
    
  end
end

% disp(mean_similarity);

fprintf('Threshold: >= %d trials. Including %d subjects.\n',nTrialThresh,sum(noNans & passTrlThresh));

%% ANOVA: factors: spaced/massed, recalled/forgotten, latency

% latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
%   0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
%   0 0.3; 0.3 0.6; 0.6 0.9; ...
%   0 0.5; 0.5 1.0; ...
%   0.3 0.8; ...
%   0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
%   0 0.8; 0.1 0.9; 0.2 1.0;
%   0 1.0];

% % 0 to 1, in 200 ms chunks
% latInd = [1:5];

% % 0.1 to 0.9, in 200 ms chunks
% latInd = [6:9];

% % 0-0.3, 0.3-0.6, 0.6-0.9
% latInd = [10:12];

% % 0-0.5, 0.5-1 *****
latInd = [13:14];

% % 0 to 1, in 600 ms chunks
% latInd = [16:20];

% % 0 to 1 in 800 ms chunks
% latInd = [21:23];

% % 0 to 1
% latInd = 24;

% =================================================

spacings = {'spac','mass'};
memConds = {'rc', 'fo'};
latency = cell(1,length(latInd));
for i = 1:length(latency)
  latency{i} = sprintf('%.1f-%.1f',latencies(latInd(1)+i-1,1),latencies(latInd(1)+i-1,2));
end
latStr = sprintf(repmat(' %s',1,length(latency)),latency{:});

factorNames = {'spacings', 'memConds', 'latency'};

nVariables = nan(size(factorNames));
keepTheseFactors = false(size(factorNames));
levelNames_teg = cell(size(factorNames)); % TEG
for c = 1:length(factorNames)
  nVariables(c) = length(eval(factorNames{c}));
  levelNames_teg{c} = eval(factorNames{c}); % TEG
  if length(eval(factorNames{c})) > 1
    keepTheseFactors(c) = true;
  end
end

variableNames = cell(1,prod(nVariables));
levelNames = cell(prod(nVariables),length(factorNames));

ses=1;
nSub = sum(passTrlThresh(:,ses));
anovaData = nan(nSub,prod(nVariables));
rmaov_data_teg = nan(nSub*prod(nVariables),length(factorNames) + 2);

fprintf('Collecting ANOVA data for %d subjects:\n\t',nSub);
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(spacings)),spacings{:}),factorNames{1});
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(memConds)),memConds{:}),factorNames{2});
fprintf('%s (%s),',latStr,factorNames{3});
fprintf('\n\tROI: %s, Freq: %s, Eig: %s, Sim: %s...',roi_str,freq_str,eig_criterion,sim_method);

lnDone = false;
vnDone = false;
subCount = 0;
rmCount = 0;
for sub = 1:length(exper.subjects)
  if all(noNans(sub,:)) && all(passTrlThresh(sub,:))
    subCount = subCount + 1;
    for ses = 1:length(exper.sesStr)
      lnCount = 0;
      vnCount = 0;
      
      for sp = 1:length(spacings)
        for mc = 1:length(memConds)
          cond_str = sprintf('img_RgH_%s_%s',memConds{mc},spacings{sp});
          
          for lat = 1:length(latInd)
            if ~lnDone
              lnCount = lnCount + 1;
              levelNames{lnCount,1} = spacings{sp};
              levelNames{lnCount,2} = memConds{mc};
              levelNames{lnCount,3} = latency{lat};
            end
            
            vnCount = vnCount + 1;
            if ~vnDone
              variableNames{vnCount} = sprintf('Y%d',vnCount);
            end
            
            anovaData(subCount,vnCount) = mean_similarity.(cond_str)(sub,ses,latInd(lat));
            
            rmCount = rmCount + 1;
            rmaov_data_teg(rmCount,:) = [mean_similarity.(cond_str)(sub,ses,latInd(lat)) sp mc lat sub];
          end
          
        end
      end
      
      lnDone = true;
      vnDone = true;
    end
  end
end

if any(~keepTheseFactors)
  factorNames = factorNames(keepTheseFactors);
  levelNames = levelNames(:,keepTheseFactors);
  nVariables = nVariables(keepTheseFactors);
  levelNames_teg = levelNames_teg(keepTheseFactors); % TEG
  
  rmaov_data_teg = rmaov_data_teg(:,[1 (find(keepTheseFactors) + 1) size(rmaov_data_teg,2)]); % TEG
  fprintf('\n\tOnly keeping factors:%s...',sprintf(repmat(' %s',1,length(factorNames)),factorNames{:}));
end
fprintf('Done.\n');

%% TEG RM ANOVA

fprintf('=======================================\n');
fprintf('This ANOVA: ROI: %s, Latency:%s\n\n',roi_str,latStr);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: ROI: %s, Latency:%s\n',roi_str,latStr);
fprintf('=======================================\n');

%% Matlab RM ANOVA

t = array2table(anovaData,'VariableNames',variableNames);

within = cell2table(levelNames,'VariableNames',factorNames);

% rm = fitrm(t,'Y1-Y8~1','WithinDesign',within);
rm = fitrm(t,sprintf('%s-%s~1',variableNames{1},variableNames{end}),'WithinDesign',within);

fprintf('=======================================\n');
fprintf('This ANOVA: ROI: %s, Latency:%s\n\n',roi_str,latStr);

margmean(rm,factorNames)
% grpstats(rm,factorNames)

MODELSPEC = sprintf(repmat('*%s',1,length(factorNames)),factorNames{:});
MODELSPEC = MODELSPEC(2:end);

% Perform repeated measures analysis of variance.
if any(nVariables) > 2
  [ranovatbl,A,C,D] = ranova(rm, 'WithinModel',MODELSPEC)
  %Show epsilon values
  %I assume that HF epsilon values > 1 (in R) are truncated to 1 by epsilon.m
  [ranovatbl,A,C,D] = ranova(rm, 'WithinModel',MODELSPEC);
  for cn = 1:length(C)
    tbl = epsilon(rm, C(cn))
  end
else
  [ranovatbl] = ranova(rm, 'WithinModel',MODELSPEC)
end

fprintf('Prev ANOVA: ROI: %s, Latency:%s\n',roi_str,latStr);
fprintf('=======================================\n');

% multcompare(rm,'spacings','By','memConds')
% multcompare(rm,'spacings','By','oldnew')
% multcompare(rm,'memConds','By','oldnew')
% % multcompare(rm,'spacings','By','oldnew','By','memConds')

pairwiseComps = nchoosek(1:length(factorNames),2);
for i = 1:size(pairwiseComps,1)
  multcompare(rm,factorNames{pairwiseComps(i,1)},'By',factorNames{pairwiseComps(i,2)})
end

%% plot RSA spacing x subsequent memory interaction

latInd = [13 14];
theseSub = noNans & passTrlThresh;

ses=1;

rc_spaced_mean = mean(mean(mean_similarity.img_RgH_rc_spac(theseSub,ses,latInd),1));
rc_massed_mean = mean(mean(mean_similarity.img_RgH_rc_mass(theseSub,ses,latInd),1));
fo_spaced_mean = mean(mean(mean_similarity.img_RgH_fo_spac(theseSub,ses,latInd),1));
fo_massed_mean = mean(mean(mean_similarity.img_RgH_fo_mass(theseSub,ses,latInd),1));

figure
hold on

s_mark = 'ko';
m_mark = 'rx';

% recalled
plot(1.9*ones(sum(theseSub(:,ses)),1), mean(mean_similarity.img_RgH_rc_spac(theseSub,ses,latInd),3),s_mark,'LineWidth',1);
plot(2.1*ones(sum(theseSub(:,ses)),1), mean(mean_similarity.img_RgH_rc_mass(theseSub,ses,latInd),3),m_mark,'LineWidth',1);

% forgotten
plot(0.9*ones(sum(theseSub(:,ses)),1), mean(mean_similarity.img_RgH_fo_spac(theseSub,ses,latInd),3),s_mark,'LineWidth',1);
plot(1.1*ones(sum(theseSub(:,ses)),1), mean(mean_similarity.img_RgH_fo_mass(theseSub,ses,latInd),3),m_mark,'LineWidth',1);

meanSizeR = 20;
meanSizeF = 20;

% recalled
hs = plot(2, mean(rc_spaced_mean),s_mark,'LineWidth',3,'MarkerSize',meanSizeR);
hm = plot(2, mean(rc_massed_mean),m_mark,'LineWidth',3,'MarkerSize',meanSizeR);

% forgotten
plot(1, mean(fo_spaced_mean),s_mark,'LineWidth',3,'MarkerSize',meanSizeF);
plot(1, mean(fo_massed_mean),m_mark,'LineWidth',3,'MarkerSize',meanSizeF);

plot([-3 3], [0 0],'k--','LineWidth',2);

hold off
axis square
% axis([0.75 2.25 75 110]);
xlim([0.75 2.25]);
ylim([-0.25 0.25]);

set(gca,'XTick', [1 2]);

set(gca,'XTickLabel',{'Forgot','Recalled'});
ylabel('Neural Similarity');

title('Spacing \times Subsequent Memory');
legend([hs, hm],{'Spaced','Massed'},'Location','North');

publishfig(gcf,0,[],[],[]);

% print(gcf,'-depsc2','~/Desktop/similarity_spacXsm.eps');

%% plot RSA spacing x time interaction

latInd = [13 14];
theseSub = noNans & passTrlThresh;

ses=1;

% rc_spaced_mean = mean([mean(D.word_RgH_rc_spac.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_rc_spac.dissim(~exper.badSub,ses,:),3)],2);
% rc_massed_mean = mean([mean(D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_rc_mass.dissim(~exper.badSub,ses,:),3)],2);
% fo_spaced_mean = mean([mean(D.word_RgH_fo_spac.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_fo_spac.dissim(~exper.badSub,ses,:),3)],2);
% fo_massed_mean = mean([mean(D.word_RgH_fo_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_fo_mass.dissim(~exper.badSub,ses,:),3)],2);

rc_spaced_mean = squeeze(mean_similarity.img_RgH_rc_spac(theseSub,ses,latInd));
rc_massed_mean = squeeze(mean_similarity.img_RgH_rc_mass(theseSub,ses,latInd));
fo_spaced_mean = squeeze(mean_similarity.img_RgH_fo_spac(theseSub,ses,latInd));
fo_massed_mean = squeeze(mean_similarity.img_RgH_fo_mass(theseSub,ses,latInd));

spac_early = mean([rc_spaced_mean(:,1) fo_spaced_mean(:,1)],2);
spac_late = mean([rc_spaced_mean(:,2) fo_spaced_mean(:,2)],2);
mass_early = mean([rc_massed_mean(:,1) fo_massed_mean(:,1)],2);
mass_late = mean([rc_massed_mean(:,2) fo_massed_mean(:,2)],2);

figure
hold on

s_mark = 'ko';
m_mark = 'rx';

% early
plot(0.9*ones(sum(theseSub(:,ses)),1), spac_early,s_mark,'LineWidth',1);
plot(1.1*ones(sum(theseSub(:,ses)),1), mass_early,m_mark,'LineWidth',1);

% late
plot(1.9*ones(sum(theseSub(:,ses)),1), spac_late,s_mark,'LineWidth',1);
plot(2.1*ones(sum(theseSub(:,ses)),1), mass_late,m_mark,'LineWidth',1);

meanSizeR = 20;
meanSizeF = 20;

% early
plot(1, mean(spac_early),s_mark,'LineWidth',3,'MarkerSize',meanSizeF);
plot(1, mean(spac_late),m_mark,'LineWidth',3,'MarkerSize',meanSizeF);

% late
hs = plot(2, mean(mass_early),s_mark,'LineWidth',3,'MarkerSize',meanSizeR);
hm = plot(2, mean(mass_late),m_mark,'LineWidth',3,'MarkerSize',meanSizeR);

% horiz
plot([-3 3], [0 0],'k--','LineWidth',2);

hold off
axis square
% axis([0.75 2.25 75 110]);
xlim([0.75 2.25]);
ylim([-0.25 0.25]);

set(gca,'XTick', [1 2]);

set(gca,'XTickLabel',{'0-0.5 sec','0.5-1.0 sec'});
ylabel('Neural Similarity');

title('Spacing \times Time');
legend([hs, hm],{'Spaced','Massed'},'Location','North');

publishfig(gcf,0,[],[],[]);

% print(gcf,'-depsc2','~/Desktop/similarity_spacXtime.eps');

% %% RMANOVA
% 
% % dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass', ...
% %   'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'};
% 
% % latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
% %   0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
% %   0 0.3; 0.3 0.6; 0.6 0.9; ...
% %   0 0.5; 0.5 1.0; ...
% %   0.3 0.8; ...
% %   0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
% %   0 0.8; 0.1 0.9; 0.2 1.0];
% 
% % % 0 to 1, in 200 ms chunks
% % latInd = [1 5];
% 
% % % 0.1 to 0.9, in 200 ms chunks
% % latInd = [6 9];
% 
% % % 0-0.3, 0.3-0.6, 0.6-0.9
% % latInd = [10 12];
% 
% % % 0-0.5, 0.5-1 *****
% % latInd = [13 14];
% 
% % % 0 to 1, in 600 ms chunks
% % latInd = [16 20];
% 
% % % 0 to 1 in 800 ms chunks
% % latInd = [21 23];
% 
% fprintf('=======================================\n');
% fprintf('Latency: %.1f-%.1f\n\n',latencies(latInd(1),1),latencies(latInd(2),2));
% 
% anovaData = [];
% 
% for sub = 1:length(exper.subjects)
%   if all(noNans(sub,:)) && all(passTrlThresh(sub,:))
%     for ses = 1:length(exper.sesStr)
%       theseData = [];
%       
%       for d = 1:length(dataTypes)
%         for lat = latInd(1):latInd(2)
%           theseData = cat(2,theseData,mean_similarity.(dataTypes{d})(sub,ses,lat));
%         end
%       end
%     end
%     anovaData = cat(1,anovaData,theseData);
%   end
% end
% 
% latStr = cell(1,length(latInd(1):latInd(2)));
% for i = 1:length(latStr)
%   latStr{i} = sprintf('%.1f-%.1f',latencies(latInd(1)+i-1,1),latencies(latInd(1)+i-1,2));
% end
% 
% % levelnames = {{'img','word'}, {'rc', 'fo'}, {'spac','mass'}, latStr};
% % varnames = {'stimType','subseqMem','spacing','time'};
% % O = teg_repeated_measures_ANOVA(anovaData, [2 2 2 length(latInd(1):latInd(2))], varnames,[],[],[],[],[],[],levelnames);
% 
% levelnames = {{'rc', 'fo'}, {'spac','mass'}, latStr};
% varnames = {'subseqMem','spacing','time'};
% O = teg_repeated_measures_ANOVA(anovaData, [2 2 length(latInd(1):latInd(2))], varnames,[],[],[],[],[],[],levelnames);
% 
% fprintf('Latency: %.1f-%.1f\n',latencies(latInd(1),1),latencies(latInd(2),2));
% fprintf('=======================================\n\n');

% %% RMANOVA - no time dimension
% 
% % dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass', ...
% %   'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'};
% 
% % latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
% %   0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
% %   0 0.3; 0.3 0.6; 0.6 0.9; ...
% %   0 0.5; 0.5 1.0; ...
% %   0.3 0.8; ...
% %   0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
% %   0 0.8; 0.1 0.9; 0.2 1.0];
% 
% % 0-0.5
% % lat = 13;
% % % 0.5-1.0
% % lat = 14;
% 
% % % 0.3-0.8
% % lat = 15;
% 
% % % 0-0.6
% % lat = 16;
% % % 0.1-0.7
% % lat = 17;
% % % 0.2-0.8 ***
% % lat = 18;
% % % 0.3-0.9
% % lat = 19;
% % % 0.4-1.0 **
% % lat = 20;
% 
% % % 0-0.8
% % lat = 21;
% % % 0.1-0.9
% % lat = 22;
% % % 0.2-1.0
% % lat = 23;
% 
% % % 0-1.0
% % lat = 24;
% 
% fprintf('=======================================\n');
% fprintf('Latency: %.1f-%.1f\n\n',latencies(lat,:));
% 
% anovaData = [];
% 
% for sub = 1:length(exper.subjects)
%   if all(noNans(sub,:)) && all(passTrlThresh(sub,:))
% %   if all(noNans(sub,:))
%     for ses = 1:length(exper.sesStr)
%       theseData = [];
%       
%       for d = 1:length(dataTypes)
%           theseData = cat(2,theseData,mean_similarity.(dataTypes{d})(sub,ses,lat));
%       end
%     end
%     anovaData = cat(1,anovaData,theseData);
%   end
% end
% 
% % no time dimension
% 
% % varnames = {'stimType','subseqMem','spacing'};
% % levelnames = {{'img','word'}, {'rc', 'fo'}, {'spac','mass'}};
% % O = teg_repeated_measures_ANOVA(anovaData, [2 2 2], varnames,[],[],[],[],[],[],levelnames);
% 
% varnames = {'subseqMem','spacing'};
% levelnames = {{'rc', 'fo'}, {'spac','mass'}};
% O = teg_repeated_measures_ANOVA(anovaData, [2 2], varnames,[],[],[],[],[],[],levelnames);
% 
% fprintf('Latency: %.1f-%.1f\n',latencies(lat,:));
% fprintf('=======================================\n\n');
