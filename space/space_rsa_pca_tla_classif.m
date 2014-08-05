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
  %'SPACE017'; % really noisy EEG, half of ICA components rejected
  'SPACE018';
  %'SPACE019'; % low trial counts
  'SPACE020';
  'SPACE021';
  'SPACE022';
  'SPACE027';
  'SPACE029';
  'SPACE037';
  'SPACE039'; % noisy EEG; original EEG analyses stopped here
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

% do we want to lowpass filter
lpfilt = true;
% do we want to subselect time?
subselect_eeg = true;

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

dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass'};

% dataTypes = {'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'};

% dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass', ...
%   'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'};

parameter = 'trial';

% latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
%   0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
%   0 0.3; 0.3 0.6; 0.6 0.9; ...
%   0 0.5; 0.5 1.0; ...
%   0.3 0.8; ...
%   0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
%   0 0.8; 0.1 0.9; 0.2 1.0;
%   0 1.0];

latencies = [0.0 0.2; 0.22 0.4; 0.42 0.6; 0.62 0.8; 0.82 1.0; ...
  0.1 0.3; 0.32 0.5; 0.52 0.7; 0.72 0.9; ...
  0 0.3; 0.32 0.6; 0.62 0.9; ...
  0 0.5; 0.52 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0];

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
% cfg_sel.avgovertime = 'yes';
cfg_sel.avgovertime = 'no';

cfg_sel.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,thisROI)}));

sim_method = 'cosine';
% sim_method = 'correlation';
% sim_method = 'spearman';

% eigenvalue options

% % keep components with eigenvalue >= 1
eig_criterion = 'kaiser';

% % compute the percent explained variance expected from each component if
% % all events are uncorrelated with each other; keep it if above this level.
% % So, each component would explain 100/n, where n is the number of
% % events/components.
% eig_criterion = 'analytic';

% keep components that cumulatively explain at least 85% of the variance
% eig_criterion = 'CV85';

%% classif: train on expo images, test on multistudy images; no classif: just calculate similarity

sesNum = 1;

% ana.trl_order.multistudy_image = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'catNum', 'targ', 'spaced', 'lag', 'presNum', 'pairOrd', 'pairNum', 'cr_recog_acc', 'cr_recall_resp', 'cr_recall_spellCorr'};

if accurateClassifSelect
  % classifier needs expo and multistudy; can include targ==-1 because
  % those are simply buffers for multistudy
  if all(ismember({'img_RgH_rc_spac', 'img_RgH_rc_mass', 'img_RgH_fo_spac', 'img_RgH_fo_mass', 'word_RgH_rc_spac', 'word_RgH_rc_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass'},dataTypes))
    data_str = 'img_word';
    ana.eventValues = {{'expo_stim','multistudy_image','multistudy_word'}};
    ana.eventValuesSplit = { ...
      {{'Face','House'} ...
      { ...
      %'img_onePres' ...
      'img_RgH_rc_spac_p1' ,'img_RgH_rc_spac_p2' ,'img_RgH_rc_mass_p1' ,'img_RgH_rc_mass_p2' ...
      'img_RgH_fo_spac_p1' ,'img_RgH_fo_spac_p2' ,'img_RgH_fo_mass_p1' ,'img_RgH_fo_mass_p2' ...
      %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
      } ...
      { ...
      %'word_onePres' ...
      'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
      'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
      %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
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
        { ...
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
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
        { ...
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        } ...
        } ...
        };
    end
  elseif all(ismember({'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass'},dataTypes))
    data_str = 'img';
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
  elseif all(ismember({'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'},dataTypes))
    data_str = 'word';
    ana.eventValues = {{'expo_stim','multistudy_word'}};
    ana.eventValuesSplit = { ...
      {{'Face','House'} ...
      { ...
      %'word_onePres' ...
      'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
      'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
      %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
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
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
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
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        } ...
        } ...
        };
    end
  end
else
  if all(ismember({'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass', 'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'},dataTypes))
    data_str = 'img_word';
    ana.eventValues = {{'multistudy_image','multistudy_word'}};
    ana.eventValuesSplit = { ...
      { ...
      { ...
      %'img_onePres' ...
      'img_RgH_rc_spac_p1' ,'img_RgH_rc_spac_p2' ,'img_RgH_rc_mass_p1' ,'img_RgH_rc_mass_p2' ...
      'img_RgH_fo_spac_p1' ,'img_RgH_fo_spac_p2' ,'img_RgH_fo_mass_p1' ,'img_RgH_fo_mass_p2' ...
      %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
      } ...
      { ...
      %'word_onePres' ...
      'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
      'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
      %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
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
        { ...
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
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
        { ...
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        } ...
        } ...
        };
    end
  elseif all(ismember({'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass'},dataTypes))
    data_str = 'img';
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
  elseif all(ismember({'word_RgH_rc_spac', 'word_RgH_rc_mass','word_RgH_fo_spac', 'word_RgH_fo_mass'},dataTypes))
    data_str = 'word';
    ana.eventValues = {{'multistudy_word'}};
    ana.eventValuesSplit = { ...
      { ...
      { ...
      %'word_onePres' ...
      'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
      'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
      %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
      } ...
      } ...
      };
    
    if allowRecallSynonyms
      ana.trl_expr = { ...
        { ...
        { ...
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        } ...
        } ...
        };
    else
      ana.trl_expr = { ...
        { ...
        { ...
        %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
        } ...
        } ...
        };
    end
  else
    error('Did not find correct event setup!');
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
% exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE030','SPACE039'}};
exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE030'}};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

%% lowpass filter and smaller segments

if lpfilt
  sampleRate = exper.sampleRate;
  lpfreq = 40;
  lofiltord = 3;
  lpfilttype = 'but';
  lpfilt_str = sprintf('lpfilt%d',lpfreq);
else
  lpfilt_str = 'lpfiltNo';
end

if subselect_eeg
  cfg_sseeg = [];
  cfg_sseeg.latency = [min(latencies(:)) max(latencies(:))];
end

if lpfilt || subselect_eeg
  for ses = 1:length(exper.sesStr)
    for typ = 1:length(ana.eventValues{ses})
      for evVal = 1:length(ana.eventValues{ses}{typ})
        for sub = 1:length(exper.subjects)
          fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
          
          if lpfilt
            % do the lowpass filter
            for i = 1:size(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(parameter),1)
              data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(parameter)(i,:,:) = ft_preproc_lowpassfilter( ...
                squeeze(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(parameter)(i,:,:)), ...
                sampleRate,lpfreq,lofiltord,lpfilttype);
            end
          end
          
          if subselect_eeg
            % select a more limited time range
            data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_selectdata(cfg_sseeg,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
          end
        end
      end
    end
  end
end

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
            nTrainTrial(dt) = size(data_tla.(sesStr).(dataTypes_train{dt}).sub(sub).data.(parameter),1);
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
        data_p1 = eval(sprintf('ft_appenddata(cfg_ad%s);',data_p1_append_str));
        data_p2 = eval(sprintf('ft_appenddata(cfg_ad%s);',data_p2_append_str));
        
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
        
        % put it in non-raw format (ft_appenddata makes it raw)
        cfg_t = [];
        cfg_t.keeptrials = 'yes';
        data_p1 = ft_timelockanalysis(cfg_t,data_p1);
        data_p2 = ft_timelockanalysis(cfg_t,data_p2);
        
        % select the requested events, times, channels, etc.
        cfg_sel.trials = p1_ind;
        data_p1 = ft_selectdata_new(cfg_sel,data_p1);
        cfg_sel.trials = p2_ind;
        data_p2 = ft_selectdata_new(cfg_sel,data_p2);
        
        %dat_p1 = data_p1.(parameter);
        %dat_p2 = data_p2.(parameter);
        
        dTypes_p1 = dTypes_p1(p1_ind);
        dTypes_p2 = dTypes_p2(p2_ind);
        
        if accurateClassifSelect
          cfg_t = [];
          cfg_t.keeptrials = 'yes';
          if isfield(cfg_sel,'trials')
            cfg_sel = rmfield(cfg_sel,'trials');
          end
          
          % attempt to classifiy study trials
          probabilityClassP1 = nan(length(p1_ind),2);
          correctClassP1 = true(size(p1_ind));
          if classifRequireP1
            for p = 1:length(p1_ind)
              dat1 = data_p1.(parameter)(p,:,:);
              dim = size(dat1);
              dat1 = reshape(dat1, dim(1), prod(dim(2:end)));
              
              Z = facehouse.test(zscore(dat1));
              probabilityClassP1(p,:) = Z;
              
              [Y,I] = max(Z,[],2);
              
              correctClassP1(p) = I == imageCategory_test(p);
            end
          end
          
          probabilityClassP2 = nan(length(p1_ind),2);
          correctClassP2 = true(size(p2_ind));
          if classifRequireP2
            for p = 1:length(p2_ind)
              dat2 = data_p2.(parameter)(p,:,:);
              dim = size(dat2);
              dat2 = reshape(dat2, dim(1), prod(dim(2:end)));
              
              Z = facehouse.test(zscore(dat2));
              probabilityClassP2(p,:) = Z;
              
              [Y,I] = max(Z,[],2);
              
              correctClassP2(p) = I == imageCategory_test(p);
            end
          end
          
          % only compare these trials
          p1_ind = find(correctClassP1 & correctClassP2);
          p2_ind = find(correctClassP1 & correctClassP2);
          
          cfg_sel.trials = p1_ind;
          data_p1 = ft_selectdata_new(cfg_sel,data_p1);
          cfg_sel.trials = p2_ind;
          data_p2 = ft_selectdata_new(cfg_sel,data_p2);
          
          dTypes_p1 = dTypes_p1(p1_ind);
          dTypes_p2 = dTypes_p2(p2_ind);
        end
        
        if ~isempty(p1_ind) && ~isempty(p2_ind)
          % unroll data for each trial in the second dimension
          dim1 = size(data_p1.(parameter));
          dim2 = size(data_p2.(parameter));
          dat_p1_p2 = cat(1,reshape(data_p1.(parameter), dim1(1), prod(dim1(2:end))),reshape(data_p2.(parameter), dim2(1), prod(dim2(2:end))));
          
          clear data_p1 data_p2
          
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
          
          % zscore the similarity values across all trials for this subject
          similarities = zscore(similarities);
          
          %             feature_vectors = data_pcaspace(:, crit_eig);
          %
          %             %%%%%%
          %             % more feature selection done here? use dummy selection for now
          %             %%%%%%
          %             select_inds = true(1, size(feature_vectors, 2));
          %
          %             evec_p1_p2_final = evec_p1_p2_crit(:, select_inds);
          %             feature_vectors = feature_vectors(:, select_inds);
          %
          %             feature_vectors = zscore(feature_vectors,0,2);
          %
          %             % normalize the vector lengths of each event
          %             feature_vectors = feature_vectors ./ repmat(sqrt(sum(feature_vectors.^2, 2)), 1, size(feature_vectors, 2));
          %
          %             % compute the similarities between each pair of events
          %             similarities_zpre = 1 - squareform(pdist(feature_vectors, 'cosine'));
          %
          %
          %             feature_vectors = data_pcaspace(:, crit_eig);
          %
          %             %%%%%%
          %             % more feature selection done here? use dummy selection for now
          %             %%%%%%
          %             select_inds = true(1, size(feature_vectors, 2));
          %
          %             evec_p1_p2_final = evec_p1_p2_crit(:, select_inds);
          %             feature_vectors = feature_vectors(:, select_inds);
          %
          %
          %             % normalize the vector lengths of each event
          %             feature_vectors = feature_vectors ./ repmat(sqrt(sum(feature_vectors.^2, 2)), 1, size(feature_vectors, 2));
          %
          %             feature_vectors = zscore(feature_vectors,0,2);
          %             % compute the similarities between each pair of events
          %             similarities_zpost = 1 - squareform(pdist(feature_vectors, 'cosine'));
          
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

saveFile = fullfile(dirs.saveDirProc,sprintf('RSA_PCA_tla_%s_%s_%s_%s_%s_%dlat_%sAvgT_%s_%s.mat',data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),cfg_sel.avgovertime,lpfilt_str,date));
fprintf('Saving: %s\n',saveFile);
save(saveFile,'exper','dataTypes','thisROI','cfg_sel','eig_criterion','sim_method','classif_str','lpfilt','latencies','similarity_all','similarity_ntrials');
fprintf('Done.\n');

%% load

analysisDate = '01-Aug-2014';

data_str = 'word';
% data_str = 'img';
% data_str = 'img_word';

% thisROI = {'center109'};
% thisROI = {'LPI2','LPS','LT','RPI2','RPS','RT'};
thisROI = {'LPS','RPS'}; % *spacing, mem x latency
% thisROI = {'LT','RT'};
% thisROI = {'LPI2','RPI2'}; % not interesting
% thisROI = {'LAS2','FS','RAS2'}; % spacing x latency
% thisROI = {'LFP','FC','RFP'};

if iscell(thisROI)
  roi_str = sprintf(repmat('%s',1,length(thisROI)),thisROI{:});
elseif ischar(thisROI)
  roi_str = thisROI;
end

lpfilt = true;
if lpfilt
  lpfreq = 40;
  %lpfilttype = 'but';
  lpfilt_str = sprintf('lpfilt%d',lpfreq);
else
  lpfilt_str = 'lpfiltNo';
end

latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
  0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
  0 0.3; 0.3 0.6; 0.6 0.9; ...
  0 0.5; 0.5 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0];

origDataType = 'tla';

% avgovertime = 'yes';
avgovertime = 'no';

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

% procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');
procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');

if ~isempty(data_str)
  rsaFile = fullfile(procDir,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%s_%dlat_%sAvgT_%s_%s_cluster.mat',origDataType,data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),avgovertime,lpfilt_str,analysisDate));
else
  rsaFile = fullfile(procDir,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%dlat_%sAvgT_%s_%s_cluster.mat',origDataType,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),avgovertime,lpfilt_str,analysisDate));
end
if exist(rsaFile,'file')
  fprintf('loading: %s...',rsaFile);
  load(rsaFile);
  fprintf('Done.\n');
else
  error('does not exist: %s',rsaFile);
end

%% stats

% nTrialThresh = 10; % 30
nTrialThresh = 12; % 27
% nTrialThresh = 14; % 25
% nTrialThresh = 16; % 22
% nTrialThresh = 17; % 21

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

% datatyp = 'word';
% datayp = 'img';

% =================================================

spacings = {'spac','mass'};
memConds = {'rc', 'fo'};
latency = cell(1,length(latInd));
latencySec = cell(1,length(latInd));
for i = 1:length(latency)
  latency{i} = sprintf('%dto%d',latencies(latInd(1)+i-1,1)*1000,latencies(latInd(1)+i-1,2)*1000);
  latencySec{i} = sprintf('%.1f-%.1f',latencies(latInd(1)+i-1,1),latencies(latInd(1)+i-1,2));
end
latStr = sprintf(repmat('_%s',1,length(latency)),latency{:});
latStr = latStr(2:end);

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
nSub = sum(noNans & passTrlThresh(:,ses));
anovaData = nan(nSub,prod(nVariables));
rmaov_data_teg = nan(nSub*prod(nVariables),length(factorNames) + 2);

fprintf('Collecting %s ANOVA data for %d subjects:\n\t',origDataType,nSub);
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(spacings)),spacings{:}),factorNames{1});
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(memConds)),memConds{:}),factorNames{2});
fprintf('%s (%s),',latStr,factorNames{3});
fprintf('\n\tROI: %s, Eig: %s, Sim: %s...',roi_str,eig_criterion,sim_method);

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
          cond_str = sprintf('%s_RgH_%s_%s',data_str,memConds{mc},spacings{sp});
          
          for lat = 1:length(latInd)
            if ~lnDone
              lnCount = lnCount + 1;
              levelNames{lnCount,1} = spacings{sp};
              levelNames{lnCount,2} = memConds{mc};
              levelNames{lnCount,3} = latency{lat};
            end
            
            vnCount = vnCount + 1;
            if ~vnDone
              variableNames{vnCount} = sprintf('y%d',vnCount);
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
fprintf('This ANOVA: %s, ROI: %s, Latency:%s, %s, %s\n\n',origDataType,roi_str,latStr,eig_criterion,sim_method);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: %s, ROI: %s, Latency:%s, %s, %s\n',origDataType,roi_str,latStr,eig_criterion,sim_method);
fprintf('=======================================\n');

%% Matlab RM ANOVA

t = array2table(anovaData,'VariableNames',variableNames);

within = cell2table(levelNames,'VariableNames',factorNames);

% rm = fitrm(t,'Y1-Y8~1','WithinDesign',within);
rm = fitrm(t,sprintf('%s-%s~1',variableNames{1},variableNames{end}),'WithinDesign',within);

fprintf('=======================================\n');
fprintf('This ANOVA: ROI: %s, Latency:%s\n\n',roi_str,latStr);

margmean(rm,factorNames)
% % grpstats(rm,factorNames)

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
% multcompare(rm,'memConds','By','spacings')
% multcompare(rm,'spacings','By','oldnew')
% multcompare(rm,'oldnew','By','spacings')
% multcompare(rm,'memConds','By','oldnew')
% multcompare(rm,'oldnew','By','memConds')

pairwiseComps = nchoosek(1:length(factorNames),2);
for i = 1:size(pairwiseComps,1)
  multcompare(rm,factorNames{pairwiseComps(i,1)},'By',factorNames{pairwiseComps(i,2)})
  multcompare(rm,factorNames{pairwiseComps(i,2)},'By',factorNames{pairwiseComps(i,1)})
end

% % these calculate the same thing
% manova(rm)
% manova(rm, 'WithinModel','separatemeans')
% 
% % set up a specific contrast
% manova(rm, 'WithinModel',[1 -1 0 0 0 0 0 0]')
% 
% % this tests the 3-way interaction
% manova(rm, 'WithinModel',[-1 1 1 -1 1 -1 -1 1]')
% 
% % MODELSPEC
% % manova(rm, 'WithinModel',sprintf('%s-%s~1',variableNames{1},variableNames{end}))
% % tests intercept
% manova(rm, 'WithinModel','1')
% 
% % manova(rm,'By','memConds') % doesn't work, between only
% 
% 
% % coeftest?
% % coeftest(rm,eye(8),[1 0 0 0 0 0 0 -1]')
% coeftest(rm,[1 1 1 1 -1 -1 -1 -1]',[1 0 0 0 0 0 0 -1]')

%% Suppose we want to compare latencies for each combination of levels of spacing and memory

% 1. Convert factors to categorical.
within2 = within;
within2.latency = categorical(within2.latency);
within2.spacings = categorical(within2.spacings);
within2.memConds = categorical(within2.memConds);

% 2. Create an interaction factor capturing each combination of levels
% of TestCond and TMS.
% within2.TestCond_TMS = within2.TestCond .* within2.TMS;
within2.spacings_memConds = within2.spacings .* within2.memConds;
within2.latency_memConds = within2.latency .* within2.memConds;
within2.spacings_latency = within2.spacings .* within2.latency;
% within2.spacings_memConds_latency = within2.spacings .* within2.memConds .* within2.latency;

% 3. Call fitrm with the modified within design.
rm2 = fitrm(t,'y1-y8~1','WithinDesign',within2);
ranovatbl2 = ranova(rm2, 'WithinModel','spacings*memConds*latency')

% 4. Use interaction factor TestCond_TMS as the 'By' variable in multcompare.
multcompare(rm2,'latency','By','spacings_memConds')
multcompare(rm2,'spacings_memConds','By','latency')


multcompare(rm2,'memConds','By','spacings_latency')
multcompare(rm2,'memConds','By','spacings_latency')
multcompare(rm2,'spacings','By','latency_memConds')

% not correct
multcompare(rm2,'spacings_memConds')
multcompare(rm2,'latency_memConds')
multcompare(rm2,'spacings_latency')

% do not correct for multiple comparisons
multcompare(rm2,'spacings_memConds_latency','ComparisonType','lsd')

%% plot RSA spacing x subsequent memory interaction

theseSub = noNans & passTrlThresh;

ses=1;

img_RgH_rc_spac = squeeze(mean_similarity.img_RgH_rc_spac(theseSub,ses,latInd));
img_RgH_fo_spac = squeeze(mean_similarity.img_RgH_fo_spac(theseSub,ses,latInd));
img_RgH_rc_mass = squeeze(mean_similarity.img_RgH_rc_mass(theseSub,ses,latInd));
img_RgH_fo_mass = squeeze(mean_similarity.img_RgH_fo_mass(theseSub,ses,latInd));

% mean across time
img_RgH_rc_spac_sub = mean(img_RgH_rc_spac,2);
img_RgH_fo_spac_sub = mean(img_RgH_fo_spac,2);
img_RgH_rc_mass_sub = mean(img_RgH_rc_mass,2);
img_RgH_fo_mass_sub = mean(img_RgH_fo_mass,2);

plotMeanLine = false;
plotSub = true;

if plotMeanLine
  s_mark = 'ko-';
  m_mark = 'rx:';
else
  s_mark = 'ko';
  m_mark = 'rx';
end

if plotSub
  subSpacing = 0.1;
  s_mark_sub = 'ko';
  m_mark_sub = 'rx';
end

meanSizeR = 20;
meanSizeF = 20;

figure
hold on

if plotSub
  % forgotten
  plot((1-subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_fo_spac_sub,s_mark_sub,'LineWidth',1);
  plot((1+subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_fo_mass_sub,m_mark_sub,'LineWidth',1);
  
  % recalled
  plot((2-subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_rc_spac_sub,s_mark_sub,'LineWidth',1);
  plot((2+subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_rc_mass_sub,m_mark_sub,'LineWidth',1);
end

if plotMeanLine
  hs = plot([mean(img_RgH_fo_spac_sub,1) mean(img_RgH_rc_spac_sub,1)],s_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  hm = plot([mean(img_RgH_fo_mass_sub,1) mean(img_RgH_rc_mass_sub,1)],m_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  
%   % recalled
%   plot(2, ,s_mark,'LineWidth',3,'MarkerSize',meanSizeR);
%   plot(2, ,m_mark,'LineWidth',3,'MarkerSize',meanSizeR);
else
  % forgotten
  plot(1, mean(img_RgH_fo_spac_sub,1),s_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  plot(1, mean(img_RgH_fo_mass_sub,1),m_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  
  % recalled
  hs = plot(2, mean(img_RgH_rc_spac_sub,1),s_mark,'LineWidth',3,'MarkerSize',meanSizeR);
  hm = plot(2, mean(img_RgH_rc_mass_sub,1),m_mark,'LineWidth',3,'MarkerSize',meanSizeR);
end

plot([-3 3], [0 0],'k--','LineWidth',2);

hold off
axis square
% axis([0.75 2.25 75 110]);
xlim([0.75 2.25]);
ylim([-0.61 0.61]);

set(gca,'XTick', [1 2]);
set(gca,'XTickLabel',{'Forgot','Recalled'});

ylabel('Neural Similarity');

title('Spacing \times Subsequent Memory');
legend([hs, hm],{'Spaced','Massed'},'Location','North');

% ticFontSize = 20;
ticFontSize = 18;
publishfig(gcf,0,ticFontSize,[],[]);

% print(gcf,'-depsc2',sprintf('~/Desktop/similarity_spacXmem_%s_%s_%s_%s_%s.eps',origDataType,roi_str,latStr,eig_criterion,sim_method));

%% plot RSA spacing x time interaction

theseSub = noNans & passTrlThresh;

ses=1;

img_RgH_rc_spac = squeeze(mean_similarity.img_RgH_rc_spac(theseSub,ses,latInd));
img_RgH_fo_spac = squeeze(mean_similarity.img_RgH_fo_spac(theseSub,ses,latInd));
img_RgH_rc_mass = squeeze(mean_similarity.img_RgH_rc_mass(theseSub,ses,latInd));
img_RgH_fo_mass = squeeze(mean_similarity.img_RgH_fo_mass(theseSub,ses,latInd));

% mean across rc/fo
img_RgH_spac = mean(cat(3,img_RgH_rc_spac,img_RgH_fo_spac),3);
img_RgH_mass = mean(cat(3,img_RgH_rc_mass,img_RgH_fo_mass),3);

plotMeanLines = true;
plotSub = true;
plotSubLines = false;

if plotMeanLines
  s_mark = 'ko-';
  m_mark = 'rx--';
else
  s_mark = 'ko';
  m_mark = 'rx';
end

meanSize = 20;

if plotSub
  subSpacing = 0.1;
  if plotSubLines
    s_mark_sub = 'ko-';
    m_mark_sub = 'rx--';
  else
    s_mark_sub = 'ko';
    m_mark_sub = 'rx';
  end
end

figure
hold on

if plotSub
  if plotSubLines
    for s = 1:sum(theseSub)
      plot(img_RgH_spac(s,:),s_mark_sub,'LineWidth',1);
      plot(img_RgH_mass(s,:),m_mark_sub,'LineWidth',1);
    end
  else
    for t = 1:length(latInd)
      if plotSub
        plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_spac(:,t),s_mark_sub,'LineWidth',1);
        plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_mass(:,t),m_mark_sub,'LineWidth',1);
      end
    end
  end
end
if plotMeanLines
  hs = plot(mean(img_RgH_spac,1),s_mark,'LineWidth',3,'MarkerSize',meanSize);
  hm = plot(mean(img_RgH_mass,1),m_mark,'LineWidth',3,'MarkerSize',meanSize);
else
  for t = 1:length(latInd)
    hs = plot(t, mean(img_RgH_spac(:,t),1),s_mark,'LineWidth',3,'MarkerSize',meanSize);
    hm = plot(t, mean(img_RgH_mass(:,t),1),m_mark,'LineWidth',3,'MarkerSize',meanSize);
  end
end

% horiz
plot([-length(latInd)-1, length(latInd)+1], [0 0],'k--','LineWidth',2);

hold off
axis square
xlim([0.75 length(latInd)+0.25]);
ylim([-0.65 0.65]);

set(gca,'XTick', 1:length(latInd));
set(gca,'XTickLabel',latencySec);
xlabel('Time (Sec)');

ylabel('Neural Similarity');

title('Spacing \times Time');
legend([hs, hm],{'Spaced','Massed'},'Location','North');

% ticFontSize = 20;
ticFontSize = 18;
publishfig(gcf,0,ticFontSize,[],[]);

% print(gcf,'-depsc2',sprintf('~/Desktop/similarity_spacXtime_%s_%s_%s_%s_%s.eps',origDataType,roi_str,latStr,eig_criterion,sim_method));

%% plot RSA spacing x memory x time interaction

theseSub = noNans & passTrlThresh;

ses=1;

img_RgH_rc_spac = squeeze(mean_similarity.img_RgH_rc_spac(theseSub,ses,latInd));
img_RgH_fo_spac = squeeze(mean_similarity.img_RgH_fo_spac(theseSub,ses,latInd));
img_RgH_rc_mass = squeeze(mean_similarity.img_RgH_rc_mass(theseSub,ses,latInd));
img_RgH_fo_mass = squeeze(mean_similarity.img_RgH_fo_mass(theseSub,ses,latInd));

plotMeanLines = true;
plotSub = false;
plotSubLines = false;

if plotMeanLines
  s_rc_mark = 'ko-';
  s_fo_mark = 'ko--';
  m_rc_mark = 'rx-';
  m_fo_mark = 'rx--';
else
  s_rc_mark = 'ko';
  s_fo_mark = 'ro';
  m_rc_mark = 'kx';
  m_fo_mark = 'rx';
end

meanSizeS = 20;
meanSizeM = 20;

if plotSub
  subSpacing = 0.2;
  if plotSubLines
    s_rc_mark_sub = 'ko-';
    s_fo_mark_sub = 'ko--';
    m_rc_mark_sub = 'rx-';
    m_fo_mark_sub = 'rx--';
  else
    s_rc_mark_sub = 'ko';
    s_fo_mark_sub = 'ro';
    m_rc_mark_sub = 'kx';
    m_fo_mark_sub = 'rx';
  end
end

figure
hold on

if plotSub
  if plotSubLines
    plotSub = false;
    for s = 1:sum(theseSub)
      % spaced
      plot(img_RgH_rc_spac(s,:),s_rc_mark_sub,'LineWidth',1);
      plot(img_RgH_fo_spac(s,:),s_fo_mark_sub,'LineWidth',1);
      % massed
      plot(img_RgH_rc_mass(s,:),m_rc_mark_sub,'LineWidth',1);
      plot(img_RgH_fo_mass(s,:),m_fo_mark_sub,'LineWidth',1);
    end
  else
    for t = 1:length(latInd)
      % spaced
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_rc_spac(:,t),s_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_fo_spac(:,t),s_fo_mark_sub,'LineWidth',1);
      % massed
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_rc_mass(:,t),m_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), img_RgH_fo_mass(:,t),m_fo_mark_sub,'LineWidth',1);
    end
  end
end
if plotMeanLines
  % spaced
  hsr = plot(mean(img_RgH_rc_spac,1),s_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  hsf = plot(mean(img_RgH_fo_spac,1),s_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  % massed
  hmr = plot(mean(img_RgH_rc_mass,1),m_rc_mark,'LineWidth',3,'MarkerSize',meanSizeM);
  hmf = plot(mean(img_RgH_fo_mass,1),m_fo_mark,'LineWidth',3,'MarkerSize',meanSizeM);
else
  for t = 1:length(latInd)
    % spaced
    hsr = plot(t,mean(img_RgH_rc_spac(:,t),1),s_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    hsf = plot(t,mean(img_RgH_fo_spac(:,t),1),s_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    % massed
    hmr = plot(t,mean(img_RgH_rc_mass(:,t),1),m_rc_mark,'LineWidth',3,'MarkerSize',meanSizeM);
    hmf = plot(t,mean(img_RgH_fo_mass(:,t),1),m_fo_mark,'LineWidth',3,'MarkerSize',meanSizeM);
  end
end

% horiz
plot([-length(latInd)-1, length(latInd)+1], [0 0],'k--','LineWidth',2);

hold off
axis square
xlim([0.75 length(latInd)+0.25]);
ylim([-0.35 0.45]);

set(gca,'XTick', 1:length(latInd));
set(gca,'XTickLabel',latencySec);
xlabel('Time (Sec)');

ylabel('Neural Similarity');

title('Spacing \times Memory \times Time');
legend([hsr, hsf, hmr, hmf],{'Spaced Rc','Spaced Fo','Massed Rc','Massed Fo'},'Location','North');

% ticFontSize = 20;
ticFontSize = 18;
publishfig(gcf,0,ticFontSize,[],[]);

% print(gcf,'-depsc2',sprintf('~/Desktop/similarity_spacXmemXtime_%s_%s_%s_%s_%s.eps',origDataType,roi_str,latStr,eig_criterion,sim_method));
