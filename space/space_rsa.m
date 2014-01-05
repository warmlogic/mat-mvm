% space RSA

%% load the analysis details

procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_eq0_art_ftManual_ftICA/tla';

subjects = {
  'SPACE001';
  'SPACE002';
%   'SPACE003';
%   'SPACE004';
%   'SPACE005';
%   'SPACE006';
%   'SPACE007';
%   %'SPACE008';
%   'SPACE009';
%   'SPACE010';
%   'SPACE011';
%   'SPACE012';
%   'SPACE013';
%   'SPACE014';
  %'SPACE015';
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
ana.eventValuesSplit = {{'Face','House'}};
ana.trl_expr = {...
  {sprintf('eventNumber == %d & i_catNum == 1 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
  sprintf('eventNumber == %d & i_catNum == 2 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim')))}};

% ana.eventValues = {{'expo_stim'}};
% ana.eventValuesSplit = {{'Face_VU','Face_SU','Face_SA','Face_VA','House_VU','House_SU','House_SA','House_VA',}};
% ana.trl_expr = {...
%   {sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 1 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 2 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 3 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 1 & expo_response == 4 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 1 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 2 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 3 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
%   sprintf('eventNumber == %d & i_catNum == 2 & expo_response == 4 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim')))}};

%% multistudy events

sesNum = 1;

% ana.trl_order.multistudy_image = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'catNum', 'targ', 'spaced', 'lag', 'presNum', 'pairOrd', 'pairNum', 'cr_recog_acc', 'cr_recall_resp', 'cr_recall_spellCorr'};

% ana.eventValues = {{'multistudy_image'},{'multistudy_word'}};
% ana.eventValuesSplit = {...
%   {'img_onePres' ...
%   'img_RgH_rc_spac_p1','img_RgH_rc_spac_p2','img_RgH_rc_mass_p1','img_RgH_rc_mass_p2' ...
%   'img_RgH_fo_spac_p1','img_RgH_fo_spac_p2','img_RgH_fo_mass_p1','img_RgH_fo_mass_p2' ...
%   %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
%   } ...
%   {'word_onePres' ...
%   'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
%   'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
%   %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
%   } ...
%   };
% ana.trl_expr = {...
%   {...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%   } ...
%   {...
%   sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
%   } ...
%   };

ana.eventValues = {{'multistudy_image'},{'multistudy_word'}};
ana.eventValuesSplit = {...
  {'img_onePres' ...
  'img_RgH_spac_p1','img_RgH_spac_p2','img_RgH_mass_p1','img_RgH_mass_p2' ...
  %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
  } ...
  {'word_onePres' ...
  'word_RgH_spac_p1','word_RgH_spac_p2','word_RgH_mass_p1','word_RgH_mass_p2' ...
  %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
  } ...
  };
ana.trl_expr = {...
  {...
  sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
  } ...
  {...
  sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
  } ...
  };

% ana.eventValues = {{'multistudy_image'},{'multistudy_word'}};
% ana.eventValuesSplit = {...
%   {'img_onePres', 'img_spac_p1','img_spac_p2','img_mass_p1','img_mass_p2'} ...
%   {'word_onePres', 'word_spac_p1','word_spac_p2','word_mass_p1','word_mass_p2'} ...
%   };
% ana.trl_expr = {...
%   {...
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
%   } ...
%   };

%% recognition events

% sesNum = 1;

% ana.trl_order.cued_recall_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'pairNum', 'recog_resp', 'recog_acc', 'recog_rt', 'new_resp', 'new_acc', 'new_rt', 'recall_resp', 'recall_spellCorr', 'recall_rt'};

% ana.eventValues = {{'cued_recall_stim'}};
% ana.eventValuesSplit = {{'RgH','CR'}};
% ana.trl_expr = {...
%   {sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}};

% ana.eventValues = {{'cued_recall_stim'}};
% ana.eventValuesSplit = {{'RgH_rc_spac','RgH_rc_mass','RgH_fo_spac','RgH_fo_mass','CR'}};
% ana.trl_expr = {...
%   {sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%   sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}};

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

%% RSA - very basic

% dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass', 'word_RgH_rc_spac', 'word_RgH_rc_mass', ...
%   'img_RgH_fo_spac', 'img_RgH_fo_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass'};

dataTypes = {'img_RgH_spac', 'img_RgH_mass', 'word_RgH_spac', 'word_RgH_mass'};

% dataTypes = {'img_RgH_rc_spac'};

% latencies = [0 1.0];
% latencies = [0 0.5; 0.5 1.0];
latencies = [0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
% latencies = [0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
% latencies = [0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9];
% latencies = [-0.2 0];

standardize = true;

distanceMetric = 'euclidean';
% distanceMetric = 'seuclidean';
% distanceMetric = 'spearman';
% distanceMetric = 'cosine';
% distanceMetric = 'correlation';

parameter = 'trial';
plotit = false;
verbose = false;

% sub = 1;
ses = 1;

% thisROI = 'LPS';
% thisROI = 'PI';
% thisROI = {'posterior'};
% thisROI = {'LPI', 'PI', 'RPI'};
% thisROI = {'LPS', 'RPS'};
% thisROI = {'LPS', 'RPS', 'LPI', 'PI', 'RPI'};
% thisROI = {'all129'};
thisROI = {'center91'};
% thisROI = 'Cz';
% thisROI = {'E70', 'E83'};
% thisROI = {'E83'};
if all(ismember(thisROI,ana.elecGroupsStr))
  elecInd = ismember(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{1}).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,thisROI)})));
elseif ~all(ismember(thisROI,ana.elecGroupsStr)) && all(ismember(thisROI,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{1}).sub(sub).data.label))
  elecInd = ismember(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{1}).sub(sub).data.label,unique(thisROI));
else
  error('Cannot find specified electrode(s)');
end

% simAcross = 'time';
% simAcross = 'chan';

% column numbers in trialinfo
phaseCountCol = 4;
stimNumCol = 6;
categNumCol = 7;

% store the distance values
D = struct;
for d = 1:length(dataTypes)
%   fprintf('%s\n',dataTypes{d});
%   D.(dataTypes{d}).nTrial
  D.(dataTypes{d}).dissim = nan(length(exper.subjects),size(latencies,1));
  D.(dataTypes{d}).nTrial = nan(length(exper.subjects),1);
end

if plotit
  if strcmp(distanceMetric,'euclidean')
    distanceScale = [0 100];
  elseif strcmp(distanceMetric,'spearman') || strcmp(distanceMetric,'correlation') || strcmp(distanceMetric,'cosine')
    distanceScale = [0 2];
  elseif strcmp(distanceMetric,'seuclidean')
    distanceScale = [0 20];
  else
    distanceScale = [];
  end
end

for lat = 1:size(latencies,1)
  fprintf('%.2f sec to %.2f sec...\n',latencies(lat,1),latencies(lat,2));
  timeS = latencies(lat,:);
  timeInd = data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{1}).sub(sub).data.time >= timeS(1) & data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{1}).sub(sub).data.time <= timeS(2) + 0.0001;
  
  for d = 1:length(dataTypes)
    dataType = dataTypes{d};
    
    fprintf('Processing %s...\n',dataType);
    
    for sub = 1:length(exper.subjects)
      
      subD = nan(size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter),1),1);
      
      if standardize
        nTrl1 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter),1);
        nTrl2 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.(parameter),1);
        data_cat = cat(1, ...
          data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd), ...
          data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd));
        dim = size(data_cat);
        dat = reshape(data_cat, dim(1), prod(dim(2:end)));
        
        mu = nanmean(dat);
        sigma = nanstd(dat);
        
        dat_z = dat;
        idx = ~isnan(mu);
        dat_z(:,idx) = bsxfun(@minus,dat(:,idx),mu(idx));
        idx = ~isnan(sigma) & ~(sigma==0);
        dat_z(:,idx) = bsxfun(@rdivide,dat_z(:,idx),sigma(idx));
        dat_z_ravel = reshape(dat_z,dim);
        
        % put it back
        %data_tla.(fn{f}).sub(sub).ses(ses).data.(parameter) = dat_z_ravel;
        dat1 = dat_z_ravel(1:nTrl1,:,:);
        dat2 = dat_z_ravel(nTrl1+1:nTrl1+nTrl2,:,:);
      end
      
      for i = 1:size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter),1)
        %p1_trlInd = 1;
        p1_trlInd = i;
        p1_phaseCount = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trialinfo(p1_trlInd,phaseCountCol);
        p1_stimNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trialinfo(p1_trlInd,stimNumCol);
        p1_categNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trialinfo(p1_trlInd,categNumCol);
        
        p2_trlInd = find(...
          data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
          data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trialinfo(:,stimNumCol) == p1_stimNum & ...
          data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trialinfo(:,categNumCol) == p1_categNum);
        %p2_trlInd = 1;
        
        if ~isempty(p2_trlInd)
          % pdist2: rows (dim 1) are observations; columns (dim 2) are variables;
          % distances are measured between observations
          
          if sum(elecInd) == 1
            p1_data = squeeze(dat1(p1_trlInd,:,:));
            p2_data = squeeze(dat2(p2_trlInd,:,:));
          elseif sum(elecInd) > 1
            p1_data = squeeze(dat1(p1_trlInd,:,:))';
            p2_data = squeeze(dat2(p2_trlInd,:,:))';
          end
          
%           if strcmp(simAcross,'time')
%             % rows = samples; cols = channels
%             if sum(elecInd) == 1
%               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd));
%               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd));
%             elseif sum(elecInd) > 1
%               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd))';
%               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd))';
%             end
%             
%             if plotit
%               xaxis = linspace(timeS(1),timeS(2),size(p1_data,1));
%               yaxis = linspace(timeS(1),timeS(2),size(p2_data,1));
%             end
%           elseif strcmp(simAcross,'chan')
%             % rows = channels; cols = samples
%             p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd));
%             p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd));
%             
%             if plotit
%               xaxis = 1:sum(elecInd);
%               yaxis = 1:sum(elecInd);
%               
%               elecLabels_x = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.label(elecInd);
%               elecLabels_y = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.label(elecInd);
%             end
%           end
          
          subD(i) = pdist2(p1_data(:)',p2_data(:)',distanceMetric);
          
          if plotit
            figure;
            if exist('distanceScale','var') && ~isempty(distanceScale)
              imagesc(xaxis,yaxis,D,distanceScale);
            else
              imagesc(xaxis,yaxis,D);
            end
            
            if strcmp(simAcross,'chan')
              set(gca,'XTickLabel',elecLabels_x);
              set(gca,'YTickLabel',elecLabels_y);
            end
            
            title(sprintf('%s (%d vs %d): phaseCount=%d stimNum=%d categNum=%d\n',strrep(dataType,'_','-'),p1_trlInd,p2_trlInd,p1_phaseCount,p1_stimNum,p1_categNum));
            if strcmp(simAcross,'time')
              xlabel('P1: Time (s)');
              ylabel('P2: Time (s)');
            elseif strcmp(simAcross,'chan')
              xlabel('P1: Electrode');
              ylabel('P2: Electrode');
            end
            
            hc = colorbar;
            set(get(hc,'YLabel'),'string','Dissimilarity');
          end
          
          %     figure;
          %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd),2)),'b');
          %     hold on
          %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd),2)),'r');
          %     %plot(linspace(timeS(1),timeS(2),size(p2_data,1)),p2_data,'r');
          %     hold off
          %     axis([timeS(1) timeS(2) -20 20]);
          
        else
          if verbose
            fprintf('%s: No p2 found for p1 phaseCount=%d stimNum=%d categNum=%d\n',dataType,p1_phaseCount,p1_stimNum,p1_categNum);
          end
        end
      end % trials
      D.(dataType).dissim(sub,lat) = nanmean(subD);
      D.(dataType).nTrial(sub) = sum(~isnan(subD));
    end % sub
  end % dataTypes
end % lat