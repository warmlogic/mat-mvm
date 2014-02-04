% space_classifier

%% load the analysis details

allowRecallSynonyms = true;

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
  'SPACE017'; % previous assessment: really noisy EEG, half of ICA components rejected
  'SPACE018';
  'SPACE019';
  'SPACE020';
  'SPACE021';
  'SPACE022';
  'SPACE027';
  'SPACE029';
  };

% only one cell, with all session names
sesNames = {'session_1'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

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

%% classifier needs expo and multistudy

% can include targ==-1 because those are simply buffers for multistudy

sesNum = 1;

ana.eventValues = {{'expo_stim','multistudy_word'}};
ana.eventValuesSplit = { ...
  {{'Face','House'} ...
  {'word_RgH_rc_spac_p2' ,'word_RgH_rc_mass_p2' ...
  'word_RgH_fo_spac_p2' ,'word_RgH_fo_mass_p2' ...
  }} ...
  };
if allowRecallSynonyms
  ana.trl_expr = {...
    {{sprintf('eventNumber == %d & i_catNum == 1 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
    sprintf('eventNumber == %d & i_catNum == 2 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim')))} ...
    {sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word')))}} ...
    };
else
  ana.trl_expr = {...
    {{sprintf('eventNumber == %d & i_catNum == 1 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
    sprintf('eventNumber == %d & i_catNum == 2 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim')))} ...
    {sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word')))}} ...
    };
end

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

%% expo

% can include targ==-1 because those are simply buffers for multistudy

sesNum = 1;

ana.eventValues = {{'expo_stim'}};
ana.eventValuesSplit = {{{'Face','House'}}};
ana.trl_expr = {...
  {{sprintf('eventNumber == %d & i_catNum == 1 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim'))), ...
  sprintf('eventNumber == %d & i_catNum == 2 & expo_response ~= 0 & rt < 3000',find(ismember(exper.eventValues{sesNum},'expo_stim')))}}};

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
if allowRecallSynonyms
  ana.trl_expr = { ...
    {{ ...
    sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
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
    sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
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
    }} ...
    };
else
  ana.trl_expr = { ...
    {{ ...
    sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
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
    sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
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
    }} ...
    };
end

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

ana.eventValues = {{'cued_recall_stim'}};
ana.eventValuesSplit = {{{'RgH_rc_spac','RgH_rc_mass','RgH_fo_spac','RgH_fo_mass','CR'}}};
if allowRecallSynonyms
  ana.trl_expr = {...
    {{sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr > 0 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr > 0 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}}};
else
  ana.trl_expr = {...
    {{sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr < 1 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr < 1 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
    sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}}};
end

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
[exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

%% lowpass filter and segment for ERPs

data_tla_backup = data_tla;

lpfilt = false;

if lpfilt
  sampleRate = 250;
  lpfreq = 40;
  lofiltord = 4;
  lpfilttype = 'but';
end

cfg_sel = [];
cfg_sel.latency = [-0.2 1.0];

for ses = 1:length(exper.sesStr)
  for typ = 1:length(ana.eventValues{ses})
    for evVal = 1:length(ana.eventValues{ses}{typ})
      for sub = 1:length(exper.subjects)
        fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
        
        if lpfilt
          for i = 1:size(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial,1)
            data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial(i,:,:) = ft_preproc_lowpassfilter( ...
              squeeze(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial(i,:,:)), ...
              sampleRate,lpfreq,lofiltord,lpfilttype);
          end
        end
        
        % select
        data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_selectdata_new(cfg_sel,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
        
      end
    end
  end
end

%% classifier - use face/house to predict image paired with P2 word

% dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass', 'word_RgH_rc_spac', 'word_RgH_rc_mass', ...
%   'img_RgH_fo_spac', 'img_RgH_fo_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass'};

thisROI = 'center91';

% latencies = [0.0 0.2; 0.1 0.3; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
latencies = [0.0 0.5; 0.5 1.0];

parameter = 'trial';

if allowRecallSynonyms
  synStr = 'syn';
else
  synStr = 'nosyn';
end

alpha = 0.2;
% alpha = 0.8;
% % L1
% alpha = 1;
% % L2
% alpha = 0;

% sub = 1;
% ses = 1;

%dataTypes_train = {'img_RgH_rc_spac_p1', 'img_RgH_rc_mass_p1', 'img_RgH_fo_spac_p1', 'img_RgH_fo_mass_p1'};

% train on expo faces and houses
dataTypes_train = {'Face', 'House'};
equateTrainTrials = true;
standardizeTrain = true;

% test on p2 word
dataTypes_test = {'word_RgH_rc_spac_p2', 'word_RgH_rc_mass_p2', 'word_RgH_fo_spac_p2', 'word_RgH_fo_mass_p2'};
equateTestTrials = false;
standardizeTest = true;
standardizeTestSeparately = false;

% data1 = data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{1}).sub(sub).data;
% data2 = data_tla.(exper.sesStr{ses}).(ana.eventValues{typ}{2}).sub(sub).data;
% data1_trials = data_tla.(exper.sesStr{ses}).Face.sub(sub).data;
% data2_trials = data_tla.(exper.sesStr{ses}).House.sub(sub).data;

testAcc = nan(length(exper.subjects),length(exper.sesStr),size(latencies,1),length(dataTypes_test));
continTab = cell(length(exper.subjects),length(exper.sesStr),size(latencies,1),length(dataTypes_test));
trainLambda = nan(length(exper.subjects),length(exper.sesStr),size(latencies,1));
trainWeights = cell(length(exper.subjects),length(exper.sesStr),size(latencies,1));

for sub = 1:length(exper.subjects)
  fprintf('%s...',exper.subjects{sub});
  for ses = 1:length(exper.sesStr)
    fprintf('%s...',exper.sesStr{ses});
    
    if ~exper.badSub(sub,ses)
      eventsFile = fullfile(dirs.dataroot,dirs.behDir,exper.subjects{sub},'events','events.mat');
      if exist(eventsFile,'file')
        fprintf('Loading events file...');
        load(eventsFile)
        fprintf('Done.\n');
      else
        error('events file not found: %s',eventsFile);
      end
      
      % equate the training categories
      trlIndTrain = cell(length(dataTypes_train),1);
      if equateTrainTrials
        nTrainTrial = nan(length(dataTypes_train),1);
        for d = 1:length(dataTypes_train)
          nTrainTrial(d) = size(data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data.(parameter),1);
        end
        fprintf('\tEquating training categories to have %d trials.\n',min(nTrainTrial));
        for d = 1:length(dataTypes_train)
          trlInd = randperm(nTrainTrial(d));
          trlIndTrain{d,1} = sort(trlInd(1:min(nTrainTrial)));
        end
      else
        fprintf('\tNot equating training category trial counts.\n');
        for d = 1:length(dataTypes_train)
          trlIndTrain{d,1} = 'all';
        end
      end
      
      % equate the testing categories
      trlIndTest = cell(length(dataTypes_test),1);
      if equateTestTrials
        nTestTrial = nan(length(dataTypes_test),1);
        for d = 1:length(dataTypes_test)
          nTestTrial(d) = size(data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data.(parameter),1);
        end
        fprintf('Equating testing categories to have %d trials.\n',min(nTestTrial));
        for d = 1:length(dataTypes_test)
          trlInd = randperm(nTestTrial(d));
          trlIndTest{d,1} = sort(trlInd(1:min(nTestTrial)));
        end
      else
        fprintf('\tNot equating testing category trial counts.\n');
        for d = 1:length(dataTypes_test)
          trlIndTest{d,1} = 'all';
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
          cfg_data.trials = trlIndTrain{d};
          data_train.(dataTypes_train{d}) = ft_selectdata_new(cfg_data,data_tla.(exper.sesStr{ses}).(dataTypes_train{d}).sub(sub).data);
        end
        
        % get the category number for each training image
        catNumCol = 7;
        imageCategory_train = data_train.(dataTypes_train{1}).trialinfo(:,catNumCol);
        for d = 2:length(dataTypes_train)
          imageCategory_train = cat(1,imageCategory_train,data_train.(dataTypes_train{d}).trialinfo(:,catNumCol));
        end
        
        % concatenate the data
        dat_train = data_train.(dataTypes_train{1}).(parameter);
        for d = 2:length(dataTypes_train)
          dat_train = cat(1,dat_train,data_train.(dataTypes_train{d}).(parameter));
        end
        
        dim = size(dat_train);
        dat_train = reshape(dat_train, dim(1), prod(dim(2:end)));
        %dat_train_ravel = reshape(dat_train,dim);
        
        % whether to standardize the training data
        if standardizeTrain
          fprintf('\t\tStandardizing the training data...');
          
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
        
        fprintf('\t\tTraining classifier...');
        %facehouse = {dml.standardizer dml.enet('family','binomial','alpha',alpha)};
        facehouse = dml.enet('family','binomial','alpha',alpha);
        facehouse = facehouse.train(dat_train,imageCategory_train);
        fprintf('Done.\n');
        trainLambda(sub,ses,lat) = facehouse.lambda;
        trainWeights{sub,ses,lat} = facehouse.weights;
        
%         cfg = [];
%         %cfg.parameter = 'trial';
%         %cfg.layout = 'GSN-HydroCel-129.sfp';
%         %cfg.method = 'crossvalidate_mm';
%         cfg.mva = dml.enet('family','binomial','alpha',alpha);
%         cfg.statistic = {'accuracy', 'binomial', 'contingency', 'confusion'};
%         cfg.compact = true;
%         cfg.verbose = true;
%         cfg.resample = false;
%         cfg.type = 'nfold';
%         cfg.nfolds = 5;
%         %cfg.type = 'split';
%         %cfg.proportion = 0.75;
%         
%         % cfg.design = imageCategory_train';
%         % cfg.design = pairedCategory_test';
%         
%         if strcmp(cfg.type,'nfold')
%           facehouse2 = dml.crossvalidator('mva',cfg.mva,'type',cfg.type,'folds',cfg.nfolds,'stat',cfg.statistic,'resample',cfg.resample,'compact',cfg.compact,'verbose',cfg.verbose);
%         elseif strcmp(cfg.type,'split')
%           facehouse2 = dml.crossvalidator('mva',cfg.mva,'type',cfg.type,'proportion',cfg.proportion,'stat',cfg.statistic,'resample',cfg.resample,'compact',cfg.compact,'verbose',cfg.verbose);
%         elseif strcmp(cfg.type,'loo') || strcmp(cfg.type,'bloo')
%           facehouse2 = dml.crossvalidator('mva',cfg.mva,'type',cfg.type,'stat',cfg.statistic,'resample',cfg.resample,'compact',cfg.compact,'verbose',cfg.verbose);
%         end
%         facehouse2 = facehouse2.train(dat_train,imageCategory_train);
        
%         keyboard
%         trainClass = struct;
%         trainClass.model = facehouse.model;
%         %trainClass.obj = facehouse;
%         
%         if length(trainClass.model) > 1
%           fn = fieldnames(trainClass.model{1});
%           for i=1:length(stat.model)
%           for k=1:length(fn)
%             if numel(trainClass.model{i}.(fn{k}))==prod(dim(2:end))
%               trainClass.model{i}.(fn{k}) = squeeze(reshape(trainClass.model{i}.(fn{k}),dim(2:end)));
%             end
%           end
%           end
%         elseif length(trainClass.model) == 1
%           fn = fieldnames(trainClass.model);
%           for k=1:length(fn)
%             if numel(trainClass.model.(fn{k}))==prod(dim(2:end))
%               trainClass.model.(fn{k}) = squeeze(reshape(trainClass.model.(fn{k}),dim(2:end)));
%             end
%           end
%         end
%         trainClass.trial = [];
%         
%         trainClass.dimord = 'chan';
%         trainClass.label = data_train.(dataTypes_train{1}).label;
%         
%         if length(trainClass.model) > 1
%           trainClass.mymodel = trainClass.model{1}.weights;
%         else
%           trainClass.mymodel = trainClass.model.weights;
%         end
%         cfg = [];
%         cfg.parameter = 'mymodel';
%         cfg.layout = 'GSN-HydroCel-129.sfp';
%         cfg.xlim = latencies(lat,:);
%         cfg.comments = '';
%         cfg.colorbar = 'yes';
%         cfg.interplimits= 'electrodes';
%         ft_topoplotER(cfg,trainClass);
        
        fprintf('\t\tSetting up testing data...\n');
        % set up the test data
        data_test = struct;
        
        sesName = 'oneDay';
        phaseName = 'multistudy';
        phaseCountCol = 4;
        trialCol = 5;
        stimNumCol = 6;
        presNumCol = 11;
        pairNumCol = 13;
        
        if standardizeTestSeparately
          for d = 1:length(dataTypes_test)
            % select the data
            %for d = 1:length(dataTypes_test)
            cfg_data.trials = trlIndTest{d};
            data_test.(dataTypes_test{d}) = ft_selectdata_new(cfg_data,data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data);
            %end
            
            % find whether p2 words were paired with a face or a house
            pairedCategory_test = [];
            %for d = 1:length(dataTypes_test)
            for i = 1:size(data_test.(dataTypes_test{d}).(parameter),1)
              phaseNum = data_test.(dataTypes_test{d}).trialinfo(i,phaseCountCol);
              
              thisEventInd = find(...
                ismember({events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.type},'STUDY_IMAGE') & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.trial] == data_test.(dataTypes_test{d}).trialinfo(i,trialCol) & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.stimNum] == data_test.(dataTypes_test{d}).trialinfo(i,stimNumCol) & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.presNum] == data_test.(dataTypes_test{d}).trialinfo(i,presNumCol) & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.pairNum] == data_test.(dataTypes_test{d}).trialinfo(i,pairNumCol) ...
                );
              
              if isempty(thisEventInd)
                thisEventInd = find(...
                  ismember({events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.type},'STUDY_WORD') & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.trial] == data_test.(dataTypes_test{d}).trialinfo(i,trialCol) & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.stimNum] == data_test.(dataTypes_test{d}).trialinfo(i,stimNumCol) & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.presNum] == data_test.(dataTypes_test{d}).trialinfo(i,presNumCol) & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.pairNum] == data_test.(dataTypes_test{d}).trialinfo(i,pairNumCol) ...
                  );
              end
              
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
            %end
            
            % get the test data
            dat_test = data_test.(dataTypes_test{d}).(parameter);
            
            % concatenate the data
            %dat_test = data_test.(dataTypes_test{1}).(parameter);
            %for d = 2:length(dataTypes_test)
            %  dat_test = cat(1,dat_test,data_test.(dataTypes_test{d}).(parameter));
            %end
            
            dim = size(dat_test);
            dat_test = reshape(dat_test, dim(1), prod(dim(2:end)));
            %dat_test_ravel = reshape(dat_test,dim);
            
            % whether to standardize the testing data
            if standardizeTest
              fprintf('\t\tStandardizing the testing data (%s)...',dataTypes_test{d});
              
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
            
            faceInd = pairedCategory_test == 1;
            houseInd = pairedCategory_test == 2;
            
            fprintf('\t\tTest data: %s\n',dataTypes_test{d});
            fprintf('\t\t%d trials: %d paired with faces, %d paired with houses\n',size(dat_test,1),sum(faceInd),sum(houseInd));
            fprintf('\t\tRunning classifier on testing data...');
            Z = facehouse.test(dat_test);
            fprintf('Done.\n');
            
            % find which category it fit best
            [Y,I] = max(Z,[],2);
            
            testAcc(sub,ses,lat,d) = mean(I == pairedCategory_test);
            
            %[I pairedCategory_test]
            fprintf('\t\tAccuracy (%s): %.4f\n',dataTypes_test{d},testAcc(sub,ses,lat,d));
            
            
            faceCorr = pairedCategory_test(faceInd) == I(faceInd);
            faceWrong = pairedCategory_test(faceInd) ~= I(faceInd);
            houseCorr = pairedCategory_test(houseInd) == I(houseInd);
            houseWrong = pairedCategory_test(houseInd) ~= I(houseInd);
            
            continTab{sub,ses,lat,d} = [sum(faceCorr), sum(faceWrong); sum(houseWrong), sum(houseCorr)];
            
            % contingency table
            fprintf('\n');
            fprintf('\t\tPredicted\n');
            fprintf('\t\t%s\t%s\n',dataTypes_train{1},dataTypes_train{2});
            fprintf('True\t%s\t%d\t%d\n',dataTypes_train{1},continTab{sub,ses,lat,d}(1,1),continTab{sub,ses,lat,d}(1,2));
            fprintf('\t%s\t%d\t%d\n',dataTypes_train{2},continTab{sub,ses,lat,d}(2,1),continTab{sub,ses,lat,d}(2,2));
            fprintf('\n');
          end % d
        else
          % standardize test data together
          
          % select the data
          for d = 1:length(dataTypes_test)
            cfg_data.trials = trlIndTest{d};
            data_test.(dataTypes_test{d}) = ft_selectdata_new(cfg_data,data_tla.(exper.sesStr{ses}).(dataTypes_test{d}).sub(sub).data);
          end
          
          % find whether p2 words were paired with a face or a house
          pairedCategory_test = [];
          for d = 1:length(dataTypes_test)
            for i = 1:size(data_test.(dataTypes_test{d}).(parameter),1)
              phaseNum = data_test.(dataTypes_test{d}).trialinfo(i,phaseCountCol);
              
              thisEventInd = find(...
                ismember({events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.type},'STUDY_IMAGE') & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.trial] == data_test.(dataTypes_test{d}).trialinfo(i,trialCol) & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.stimNum] == data_test.(dataTypes_test{d}).trialinfo(i,stimNumCol) & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.presNum] == data_test.(dataTypes_test{d}).trialinfo(i,presNumCol) & ...
                [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.pairNum] == data_test.(dataTypes_test{d}).trialinfo(i,pairNumCol) ...
                );
              
              if isempty(thisEventInd)
                thisEventInd = find(...
                  ismember({events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.type},'STUDY_WORD') & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.trial] == data_test.(dataTypes_test{d}).trialinfo(i,trialCol) & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.stimNum] == data_test.(dataTypes_test{d}).trialinfo(i,stimNumCol) & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.presNum] == data_test.(dataTypes_test{d}).trialinfo(i,presNumCol) & ...
                  [events.(sesName).(sprintf('%s_%d',phaseName,phaseNum)).data.pairNum] == data_test.(dataTypes_test{d}).trialinfo(i,pairNumCol) ...
                  );
              end
              
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
          
          % concatenate the data
          dat_test = data_test.(dataTypes_test{1}).(parameter);
          testDataIndices = ones(size(data_test.(dataTypes_test{1}).(parameter),1),1);
          for d = 2:length(dataTypes_test)
            dat_test = cat(1,dat_test,data_test.(dataTypes_test{d}).(parameter));
            testDataIndices = cat(1,testDataIndices,d * ones(size(data_test.(dataTypes_test{d}).(parameter),1),1));
          end
          
          dim = size(dat_test);
          dat_test = reshape(dat_test, dim(1), prod(dim(2:end)));
          %dat_test_ravel = reshape(dat_test,dim);
          
          % whether to standardize the testing data
          if standardizeTest
            fprintf('\t\tStandardizing the testing data (%s)...',dataTypes_test{d});
            
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
          
          for d = 1:length(dataTypes_test)
            datInd = testDataIndices == d;
            thisPairedCategory_test = pairedCategory_test(datInd,:);
            faceInd = thisPairedCategory_test == 1;
            houseInd = thisPairedCategory_test == 2;
            
            fprintf('\t\tTest data: %s\n',dataTypes_test{d});
            fprintf('\t\t%d trials: %d paired with faces, %d paired with houses\n',size(dat_test(datInd,:),1),sum(faceInd),sum(houseInd));
            fprintf('\t\tRunning classifier on testing data...');
            Z = facehouse.test(dat_test(datInd,:));
            fprintf('Done.\n');
            
            % find which category it fit best
            [Y,I] = max(Z,[],2);
            
            testAcc(sub,ses,lat,d) = mean(I == pairedCategory_test(datInd,:));
            
            %[I pairedCategory_test]
            fprintf('\t\tAccuracy (%s): %.4f\n',dataTypes_test{d},testAcc(sub,ses,lat,d));
            
            faceCorr = thisPairedCategory_test(faceInd) == I(faceInd);
            faceWrong = thisPairedCategory_test(faceInd) ~= I(faceInd);
            houseCorr = thisPairedCategory_test(houseInd) == I(houseInd);
            houseWrong = thisPairedCategory_test(houseInd) ~= I(houseInd);
            
            continTab{sub,ses,lat,d} = [sum(faceCorr), sum(faceWrong); sum(houseWrong), sum(houseCorr)];
            
            % contingency table
            fprintf('\n');
            fprintf('\t\tPredicted\n');
            fprintf('\t\t%s\t%s\n',dataTypes_train{1},dataTypes_train{2});
            fprintf('True\t%s\t%d\t%d\n',dataTypes_train{1},continTab{sub,ses,lat,d}(1,1),continTab{sub,ses,lat,d}(1,2));
            fprintf('\t%s\t%d\t%d\n',dataTypes_train{2},continTab{sub,ses,lat,d}(2,1),continTab{sub,ses,lat,d}(2,2));
            fprintf('\n');
          end % d
          
        end
      end % lat
    else
      fprintf('bad subject!\n');
    end % badSub
  end % ses
end % sub

save(sprintf('fhClass_100HzLP_%s_alpha%s_4Feb2014.mat',synStr,strrep(num2str(alpha),'.','')),'testAcc','continTab','trainLambda','trainWeights');

%% testAcc

% ses = 1;
% lat = 2;
% mean(testAcc(~exper.badSub(:,ses),ses,lat,1),1)
% mean(testAcc(~exper.badSub(:,ses),ses,lat,2),1)
% mean(testAcc(~exper.badSub(:,ses),ses,lat,3),1)
% mean(testAcc(~exper.badSub(:,ses),ses,lat,4),1)

alpha = 0.05;
tails = 'both';
for ses = 1:length(exper.sesStr)
  fprintf('%s...',exper.sesStr{ses});
  for lat = 1:size(latencies,1)
    fprintf('%.2fs to %.2fs...\n',latencies(lat,1),latencies(lat,2));
    
    % recalled
    fprintf('Recalled\n');
    spaced = 1;
    massed = 2;
    data1_str = dataTypes_test{spaced};
    data2_str = dataTypes_test{massed};
    data1 = testAcc(~exper.badSub(:,ses),ses,lat,spaced);
    data2 = testAcc(~exper.badSub(:,ses),ses,lat,massed);
    
%     data2_str = 'chance';
%     data2 = 0.5*ones(sum(~exper.badSub(:,ses)),1);
%     d = mm_effect_size('single',data1,0.5);
%     [h, p, ci, stats] = ttest(data1,data2,alpha,tails);
    
    d = mm_effect_size('within',data1,data2);
    [h, p, ci, stats] = ttest(data1,data2,alpha,tails);
    
    fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
      data1_str, ...
      mean(data1), ...
      std(data1) / sqrt(length(data1)), ...
      data2_str, ...
      mean(data2), ...
      std(data2) / sqrt(length(data2)), ...
      stats.df, ...
      stats.tstat, ...
      d, ...
      std(data1 - data2),...
      std(data1 - data2) / sqrt(length(data1)),...
      p);
    
    % forgotten
    fprintf('Forgotten\n');
    spaced = 3;
    massed = 4;
    data1_str = dataTypes_test{spaced};
    data2_str = dataTypes_test{massed};
    data1 = testAcc(~exper.badSub(:,ses),ses,lat,spaced);
    data2 = testAcc(~exper.badSub(:,ses),ses,lat,massed);
    
    d = mm_effect_size('within',data1,data2);
    [h, p, ci, stats] = ttest(data1,data2,alpha,tails);
    
    fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
      data1_str, ...
      mean(data1), ...
      std(data1) / sqrt(length(data1)), ...
      data2_str, ...
      mean(data2), ...
      std(data2) / sqrt(length(data2)), ...
      stats.df, ...
      stats.tstat, ...
      d, ...
      std(data1 - data2),...
      std(data1 - data2) / sqrt(length(data1)),...
      p);
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

