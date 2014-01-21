% space_classifier

%% load the analysis details

procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';

subjects = {
  'SPACE001';
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
  'SPACE017'; % really noisy EEG, half of ICA components rejected
  'SPACE018';
  'SPACE019';
  };

% only one cell, with all session names
sesNames = {'session_1'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = false;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% list the event values to analyze; specific to each experiment

% this is useful for when there are multiple types of event values, for
% example, hits and CRs in two conditions. You don't have to enter anything
% if you just want all events from exper.eventValues together in a single
% cell because it will get set to {exper.eventValues}, but it needs to be a
% cell containing a cell of eventValue strings

% this is only used by mm_ft_checkCondComps to create pairwise combinations
% either within event types {'all_within_types'} or across all event types
% {'all_across_types'}; mm_ft_checkCondComps is called within subsequent
% analysis functions

%% expo

% can include targ==-1 because those are simply buffers for multistudy

sesNum = 1;

ana.eventValues = {{'expo_stim'}};
ana.eventValuesSplit = {{{'Face','House'}}};
ana.trl_expr = {...
  {{sprintf('eventNumber == %d & i_catNum == 1 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
  sprintf('eventNumber == %d & i_catNum == 2 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim')))}}};

% ana.eventValues = {{'expo_stim'}};
% ana.eventValuesSplit = {{{'Face_VU','Face_SU','Face_SA','Face_VA','House_VU','House_SU','House_SA','House_VA',}}};
% ana.trl_expr = {...
%   {{sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 1 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 2 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 3 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 4 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 1 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 2 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 3 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 4 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim')))}}};

%% multistudy events

sesNum = 1;

% ana.trl_order.multistudy_image = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'catNum', 'targ', 'spaced', 'lag', 'presNum', 'pairOrd', 'pairNum', 'cr_recog_acc', 'cr_recall_resp', 'cr_recall_spellCorr'};

ana.eventValues = {{'multistudy_image','multistudy_word'}};
ana.eventValuesSplit = { ...
  {{'img_onePres' ...
  'img_RgH_rc_spac_p1','img_RgH_rc_spac_p2','img_RgH_rc_mass_p1','img_RgH_rc_mass_p2' ...
  'img_RgH_fo_spac_p1','img_RgH_fo_spac_p2','img_RgH_fo_mass_p1','img_RgH_fo_mass_p2' ...
  %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
  } ...
  {'word_onePres' ...
  'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
  'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
  %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
  }} ...
  };
ana.trl_expr = { ...
  {{ ...
  sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
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
  sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  }} ...
  };

% ana.eventValues = {{'multistudy_image', 'multistudy_word'}};
% ana.eventValuesSplit = { ...
%   {{'img_onePres' ...
%   'img_RgH_spac_p1','img_RgH_spac_p2','img_RgH_mass_p1','img_RgH_mass_p2' ...
%   %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
%   } ...
%   {'word_onePres' ...
%   'word_RgH_spac_p1','word_RgH_spac_p2','word_RgH_mass_p1','word_RgH_mass_p2' ...
%   %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
%   }} ...
%   };
% ana.trl_expr = { ...
%   {{ ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   } ...
%   { ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   }} ...
%   };

% ana.eventValues = {{'multistudy_image','multistudy_word'}};
% ana.eventValuesSplit = {...
%   {{'img_onePres', 'img_spac_p1','img_spac_p2','img_mass_p1','img_mass_p2'} ...
%   {'word_onePres', 'word_spac_p1','word_spac_p2','word_mass_p1','word_mass_p2'}} ...
%   };
% ana.trl_expr = {...
%   {{...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   } ...
%   {...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   }} ...
%   };

%% recognition events

% sesNum = 1;

% ana.trl_order.cued_recall_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'pairNum', 'recog_resp', 'recog_acc', 'recog_rt', 'new_resp', 'new_acc', 'new_rt', 'recall_resp', 'recall_spellCorr', 'recall_rt'};

% ana.eventValues = {{'cued_recall_stim'}};
% ana.eventValuesSplit = {{{'RgH','CR'}}};
% ana.trl_expr = {...
%   {{sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}}};

% ana.eventValues = {{'cued_recall_stim'}};
% ana.eventValuesSplit = {{{'RgH_rc_spac','RgH_rc_mass','RgH_fo_spac','RgH_fo_mass','CR'}}};
% ana.trl_expr = {...
%   {{sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}}};

%% load in the subject data

% % make sure ana.eventValues is set properly
% if ~iscell(ana.eventValues{1})
%   ana.eventValues = {ana.eventValues};
% end
% if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
%   ana.eventValues = {exper.eventValues};
% end

% [data_tla,exper] = mm_ft_loadSubjectData(exper,dirs,ana,'tla',1,'trialinfo');
[data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',1,'trialinfo');

% %% get rid of the bad channels
% 
% cfg = [];
% cfg.printRoi = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% [data_tla] = mm_rmBadChan(cfg,exper,ana,data_tla);

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {{}};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,3,[],'vert');

%% classifier - use face/house to predict image paired with P2 word

% dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass', 'word_RgH_rc_spac', 'word_RgH_rc_mass', ...
%   'img_RgH_fo_spac', 'img_RgH_fo_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass'};

thisROI = 'center91';

% latencies = [0.0 0.2; 0.1 0.3; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
latencies = [0.0 0.5; 0.5 1.0];

parameter = 'trial';

% sub = 1;
ses = 1;

% data1 = data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{1}).sub(sub).data;
% data2 = data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{2}).sub(sub).data;
% data1_trials = data_tla.(exper.sesStr{ses}).Face.sub(sub).data;
% data2_trials = data_tla.(exper.sesStr{ses}).House.sub(sub).data;

for sub = 1:length(exper.subjects)
  fprintf('Subject: %s...\n',exper.subjects{sub});
  
  eventsFile = fullfile(dirs.dataroot,dirs.behDir,exper.subjects{sub},'events','events.mat');
  if exist(eventsFile,'file')
    load(eventsFile)
  else
    error('events file not found: %s',eventsFile);
  end
  
  %dataTypes_train = {'img_RgH_rc_spac_p1', 'img_RgH_rc_mass_p1', 'img_RgH_fo_spac_p1', 'img_RgH_fo_mass_p1'};
  dataTypes_train = {'Face', 'House'};
  
  % get the category number for each training image
  catNumCol = 7;
  imageCategory_train = data_tla.(exper.sesStr{ses}).(dataTypes_train{1}).sub(sub).data.trialinfo(:,catNumCol);
  for d = 2:length(dataTypes_train)
    imageCategory_train = cat(1,imageCategory_train,data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data.trialinfo(:,catNumCol));
  end
  
  % test on p2 word
  dataTypes_test = {'word_RgH_rc_spac_p2', 'word_RgH_rc_mass_p2', 'word_RgH_fo_spac_p2', 'word_RgH_fo_mass_p2'};
  
  sesName = 'oneDay';
  phaseName = 'multistudy';
  phaseCountCol = 4;
  trialCol = 5;
  stimNumCol = 6;
  presNumCol = 11;
  pairNumCol = 13;
  pairedCategory_test = [];
  for d = 1:length(dataTypes_test)
    for i = 1:size(data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.(cfg_stat.parameter),1)
      phaseNum = data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.trialinfo(i,phaseCountCol);
      
      thisEventInd = find(...
        [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.trial] == data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.trialinfo(i,trialCol) & ...
        [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.stimNum] == data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.trialinfo(i,stimNumCol) & ...
        [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.presNum] == data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.trialinfo(i,presNumCol) & ...
        [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.pairNum] == data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.trialinfo(i,pairNumCol) ...
        );
      
      if ~isempty(thisEventInd) && length(thisEventInd) == 1
        pairCatNum = events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data(thisEventInd+1).catNum;
        pairedCategory_test = cat(1,pairedCategory_test,pairCatNum);
      elseif ~isempty(thisEventInd) && length(thisEventInd) > 1
        %keyboard
        error('Found more than one event!');
      elseif isempty(thisEventInd)
        %keyboard
        error('No events found!');
      end
    end
  end
  
  for lat = 1:size(latencies,1)
    fprintf('\t%.2fs to %.2fs...\n',latencies(lat,1),latencies(lat,2));
    
    cfg_data = [];
    cfg_data.latency = latencies(lat,:);
    cfg_data.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,thisROI)});
    cfg_data.avgoverchan = 'no';
    cfg_data.avgovertime = 'no';
    % cfg_data.avgovertime = 'yes';
    
    data_train = struct;
    
    % select the data
    for d = 1:length(dataTypes_train)
      data_train.(dataTypes_train{d}) = ft_selectdata_new(cfg_data,data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data);
    end
    
    % concatenate the data
    dat_train = data_train.(dataTypes_train{1}).(parameter);
    for d = 2:length(dataTypes_train)
      dat_train = cat(1,dat_train,data_train.(dataTypes_train{d}).(parameter));
    end
    
    dim = size(dat_train);
    dat_train = reshape(dat_train, dim(1), prod(dim(2:end)));
    %dat_train_ravel = reshape(dat_train,dim);
    
    % whether to equate the training data
    equateTrainTrials = true;
    if equateTrainTrials
      nTrainCateg = unique(imageCategory_train);
      nTrainTrial = nan(length(nTrainCateg),1);
      for i = 1:length(nTrainCateg)
        nTrainTrial(i) = sum(imageCategory_train == i);
      end
      for i = 1:length(nTrainCateg)
        if nTrainTrial(i) > min(nTrainTrial)
          thisCatInd = find(imageCategory_train == i);
          thisCatInd = thisCatInd(randperm(length(thisCatInd)));
          tr_to_keep = thisCatInd(1:min(nTrainTrial));
          
          tr_to_keep = sort(cat(1,tr_to_keep,find(imageCategory_train ~= i)));
          imageCategory_train = imageCategory_train(tr_to_keep);
          dat_train = dat_train(tr_to_keep,:);
        end
      end
    end
    
    % whether to standardize the data
    standardize = true;
    if standardize
      fprintf('Standardizing the data...');
      
      m = dml.standardizer;
      m = m.train(dat_train);
      dat_train = m.test(dat_train);
      %dat_train_ravel = reshape(dat_train,dim);
      
      %       mu = nanmean(dat);
      %       sigma = nanstd(dat);
      %
      %       dat = dat;
      %       idx = ~isnan(mu);
      %       dat(:,idx) = bsxfun(@minus,dat(:,idx),mu(idx));
      %       idx = ~isnan(sigma) & ~(sigma==0);
      %       dat(:,idx) = bsxfun(@rdivide,dat(:,idx),sigma(idx));
      %       dat_ravel = reshape(dat,dim);
      %
      %       % put it back
      %       data1.(parameter) = dat_z_ravel(1:size(data1.(parameter),1),:,:);
      %       data2.(parameter) = dat_z_ravel(size(data2.(parameter),1)+1:size(data1.(parameter),1)+size(data2.(parameter),1),:,:);
      
      fprintf('Done.\n');
    end
    
    %facehouse = {dml.standardizer dml.enet('family','binomial','alpha',0.2)};
    facehouse = dml.enet('family','binomial','alpha',0.2);
    facehouse = facehouse.train(dat_train,imageCategory_train);
    
%     cfg = [];
%     %cfg.parameter = 'trial';
%     %cfg.layout = 'GSN-HydroCel-129.sfp';
%     %cfg.method = 'crossvalidate_mm';
%     cfg.mva = dml.enet('family','binomial','alpha',0.2);
%     cfg.statistic = {'accuracy', 'binomial', 'contingency', 'confusion'};
%     cfg.compact = true;
%     cfg.verbose = true;
%     cfg.resample = false;
%     cfg.type = 'nfold';
%     cfg.nfolds = 5;
%     %cfg.type = 'split';
%     %cfg.proportion = 0.75;
%     
%     if strcmp(cfg.type,'nfold')
%       cv = dml.crossvalidator('mva',cfg.mva,'type',cfg.type,'folds',cfg.nfolds,'stat',cfg.statistic,'resample',cfg.resample,'compact',cfg.compact,'verbose',cfg.verbose);
%     elseif strcmp(cfg.type,'split')
%       cv = dml.crossvalidator('mva',cfg.mva,'type',cfg.type,'proportion',cfg.proportion,'stat',cfg.statistic,'resample',cfg.resample,'compact',cfg.compact,'verbose',cfg.verbose);
%     elseif strcmp(cfg.type,'loo') || strcmp(cfg.type,'bloo')
%       cv = dml.crossvalidator('mva',cfg.mva,'type',cfg.type,'stat',cfg.statistic,'resample',cfg.resample,'compact',cfg.compact,'verbose',cfg.verbose);
%     end
% 
%     % cfg.design = imageCategory_train';
%     % cfg.design = pairedCategory_test';
%     
%     %cv = cv.train(dat_train,imageCategory_train);
    
    % set up the test data
    data_test = struct;
    
    % select the data
    for d = 1:length(dataTypes_test)
      data_test.(dataTypes_test{d}) = ft_selectdata_new(cfg_data,data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data);
    end
    
    % concatenate the data
    dat_test = data_test.(dataTypes_test{1}).(parameter);
    for d = 2:length(dataTypes_test)
      dat_test = cat(1,dat_test,data_test.(dataTypes_test{d}).(parameter));
    end
    
    dim = size(dat_test);
    dat_test = reshape(dat_test, dim(1), prod(dim(2:end)));
    %dat_test_ravel = reshape(dat_test,dim);
    
    % whether to equate the test data
    equateTestTrials = true;
    if equateTestTrials
      nTestCateg = unique(pairedCategory_test);
      nTestTrial = nan(length(nTestCateg),1);
      for i = 1:length(nTestCateg)
        nTestTrial(i) = sum(pairedCategory_test == i);
      end
      for i = 1:length(nTestCateg)
        if nTestTrial(i) > min(nTestTrial)
          thisCatInd = find(pairedCategory_test == i);
          thisCatInd = thisCatInd(randperm(length(thisCatInd)));
          tr_to_keep = thisCatInd(1:min(nTestTrial));
          
          tr_to_keep = sort(cat(1,tr_to_keep,find(pairedCategory_test ~= i)));
          pairedCategory_test = pairedCategory_test(tr_to_keep);
          dat_test = dat_test(tr_to_keep,:);
        end
      end
    end
    
    % whether to standardize the data
    standardize = true;
    if standardize
      fprintf('Standardizing the data...');
      
      m = dml.standardizer;
      m = m.train(dat_test);
      dat_test = m.test(dat_test);
      %dat_test_ravel = reshape(dat_test,dim);
      
      %       mu = nanmean(dat);
      %       sigma = nanstd(dat);
      %
      %       dat = dat;
      %       idx = ~isnan(mu);
      %       dat(:,idx) = bsxfun(@minus,dat(:,idx),mu(idx));
      %       idx = ~isnan(sigma) & ~(sigma==0);
      %       dat(:,idx) = bsxfun(@rdivide,dat(:,idx),sigma(idx));
      %       dat_ravel = reshape(dat,dim);
      %
      %       % put it back
      %       data1.(parameter) = dat_z_ravel(1:size(data1.(parameter),1),:,:);
      %       data2.(parameter) = dat_z_ravel(size(data2.(parameter),1)+1:size(data1.(parameter),1)+size(data2.(parameter),1),:,:);
      
      fprintf('Done.\n');
    end
    
    %cv = cv.train(dat_test,pairedCategory_test);
    
    Z = facehouse.test(dat_test);
    [Y,I] = max(Z,[],2);
    
    %[I pairedCategory_test]
    fprintf('\tAccuracy: %.4f\n',mean(I == pairedCategory_test));
    
  end
end


%% classification from the tutorial

cfg_ana = [];

cfg_ana.latencies = [0.0 0.2; 0.1 0.3; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
lat = 2;
typ = 1;

sub = 1;
ses = 1;

% data1 = data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{1}).sub(sub).data;
% data2 = data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{2}).sub(sub).data;
data1 = data_tla.(exper.sesStr{ses}).img_onePres.sub(sub).data;
data2 = data_tla.(exper.sesStr{ses}).word_onePres.sub(sub).data;

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

% data1 = ft_selectdata_new(cfg,data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{1}).sub(sub).data);
% data2 = ft_selectdata_new(cfg,data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{2}).sub(sub).data);

% data1 = ft_selectdata_old(data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{1}).sub(sub).data,...
%   'param',cfg.parameter,...
%   'toilim',cfg.latency,...
%   'channel',cfg.channel,...
%   'avgoverchan',cfg.avgoverchan,...
%   'avgovertime',cfg.avgovertime);
% 
% data2 = ft_selectdata_old(data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{2}).sub(sub).data,...
%   'param',cfg.parameter,...
%   'toilim',cfg.latency,...
%   'channel',cfg.channel,...
%   'avgoverchan',cfg.avgoverchan,...
%   'avgovertime',cfg.avgovertime);


cfg.design = [ones(size(data1.(cfg.parameter),1),1); 2*ones(size(data2.(cfg.parameter),1),1)]';

cfg.mva = {dml.standardizer dml.enet('family','binomial','alpha',0.2)};
%cfg.mva = {dml.standardizer dml.enet('family','binomial','alpha',1)};
%cfg.mva = {dml.crossvalidator('mva',{dml.enet('family','binomial','alpha',0.2)},'stat',{'accuracy','binomial','contingency'},'resample',true,'verbose',true)};
% cfg.mva = dml.analysis({dml.enet('family','binomial','alpha',0.2)});
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


%     % whether to equate the training data
%     equateTrainTrials = false;
%     if equateTrainTrials
%       minEv = Inf;
%       min_d = [];
%       for d = 1:length(dataTypes_train)
%         if size(data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data.(parameter),1) < minEv
%           minEv = size(data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data.(parameter),1);
%           min_d = d;
%         end
%       end
%       for d = 1:length(dataTypes_train)
%         cfg = [];
%         if d == min_d
%           cfg.trials = 'all';
%         else
%           cfg.trials = randperm(minEv);
%         end
%         data_train.(dataTypes_train{d}) = ft_selectdata_new(cfg,data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data);
%       end
%     else
%       for d = 1:length(dataTypes_train)
%         data_train.(dataTypes_train{d}) = data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data;
%       end
%     end
    

%     % whether to equate the test data
%     equateTestTrials = false;
%     if equateTestTrials
%       minEv = Inf;
%       min_d = [];
%       for d = 1:length(dataTypes_test)
%         if size(data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.(parameter),1) < minEv
%           minEv = size(data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.(parameter),1);
%           min_d = d;
%         end
%       end
%       for d = 1:length(dataTypes_test)
%         cfg = [];
%         if d == min_d
%           cfg.trials = 'all';
%         else
%           cfg.trials = randperm(minEv);
%         end
%         data_test.(dataTypes_test{d}) = ft_selectdata_new(cfg,data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data);
%       end
%     else
%       for d = 1:length(dataTypes_test)
%         data_test.(dataTypes_test{d}) = data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data;
%       end
%     end
    
