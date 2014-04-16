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

% procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
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
  'SPACE030'; % low trial counts
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

% ana.eventValues = {{'multistudy_image','multistudy_word'}};
ana.eventValues = {{'multistudy_word'}};
ana.eventValuesSplit = { ...
  { ...
%   { ...
%   %'img_onePres' ...
%   'img_RgH_rc_spac_p1','img_RgH_rc_spac_p2','img_RgH_rc_mass_p1','img_RgH_rc_mass_p2' ...
%   'img_RgH_fo_spac_p1','img_RgH_fo_spac_p2','img_RgH_fo_mass_p1','img_RgH_fo_mass_p2' ...
%   %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
%   } ...
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
%     { ...
%     %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     } ...
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
%     { ...
%     %sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     } ...
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
% if allowRecallSynonyms
%   ana.trl_expr = {...
%     {{sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr > 0 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr > 0 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 0 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}}};
% else
%   ana.trl_expr = {...
%     {{sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr == 1 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr < 1 & spaced == 1 & lag > 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 1 & recog_resp == 1 & recog_acc == 1 & recog_rt < 3000 & recall_spellCorr < 1 & spaced == 0 & lag == 0',find(ismember(exper.eventValues{sesNum},'cued_recall_stim'))), ...
%     sprintf('eventNumber == %d & targ == 0 & recog_resp == 2 & recog_acc == 1 & recog_rt < 3000 & new_resp ~= 0 & new_acc == 1',find(ismember(exper.eventValues{sesNum},'cued_recall_stim')))}}};
% end

%% load in the subject data

% % make sure ana.eventValues is set properly
% if ~iscell(ana.eventValues{1})
%   ana.eventValues = {ana.eventValues};
% end
% if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
%   ana.eventValues = {exper.eventValues};
% end

keeptrials = true;
% [data_tla,exper] = mm_ft_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');
[data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');

% %% get rid of the bad channels
% 
% cfg = [];
% cfg.printRoi = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% [data_tla] = mm_rmBadChan(cfg,exper,ana,data_tla);

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

%% decide who to kick out based on trial counts

% Subjects with bad behavior
% exper.badBehSub = {{}};
exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE039'}};

% SPACE019 has particularly low distance (high similarity) values

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

%% convert to scalp current density

if ~exist('data_tla_backup','var')
  data_tla_backup = data_tla;
else
  error('data_tla_backup already exists, don''t want to overwrite it');
end

cfg_sel = [];
cfg_sel.latency = [-0.2 1.0];

cfg_scd = [];
cfg_scd.elec = ana.elec;
cfg_scd.method = 'spline';

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    for typ = 1:length(ana.eventValues{ses})
      for evVal = 1:length(ana.eventValues{ses}{typ})
        if ~exper.badSub(sub)
          fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
          
          % select
          data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_selectdata_new(cfg_sel,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
          
          % scalp current density
          data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_scalpcurrentdensity(cfg_scd,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
        end
      end
    end
  end
end

% cfg_tla = [];
% cfg_tla.keeptrials = 'no';
% c = ft_timelockanalysis(cfg_tla,b);
% 
% cfg_t = [];
% cfg_t.xlim = [0.1 0.3];
% % cfg_t.xlim = [0.5 0.8];
% cfg_t.layout = ana.elec;
% ft_topoplotER(cfg_t,c);

% fn = fieldnames(data_tla.session_1);
% for i = 1:length(fn)
%   if ~ismember(fn{i},ana.eventValues{1}{1})
%     data_tla.session_1 = rmfield(data_tla.session_1,fn{i});
%   end
% end

%% RSA - very basic

dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass', 'word_RgH_rc_spac', 'word_RgH_rc_mass', ...
  'img_RgH_fo_spac', 'img_RgH_fo_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass'};

% dataTypes = {'word_RgH_rc_spac', 'word_RgH_rc_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass'};

% % testing out SCD
% dataTypes = {'word_RgH_rc_spac', 'word_RgH_rc_mass'};

% dataTypes = {'img_RgH_rc_spac'};

% dataTypes = {'img_RgH_spac', 'img_RgH_mass', 'word_RgH_spac', 'word_RgH_mass'};
% dataTypes = {'img_RgH_spac', 'img_RgH_mass'};
% dataTypes = {'Face', 'House'};

% compare the different types of stimuli
% dataTypes = {{'img_RgH_spac', 'word_RgH_spac'}, {'img_RgH_mass', 'word_RgH_mass'}};



% latencies = [0 1.0];
% latencies = [0 0.5; 0.5 1.0];
% latencies = [0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0];
% latencies = [-0.2 -0.1; -0.1 0; 0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
% latencies = [-0.2 -0.1; -0.1 0; 0 0.1; 0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5; 0.5 0.6; 0.6 0.7; 0.7 0.8; 0.8 0.9; 0.9 1.0];
latencies = [0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9];
% latencies = [-0.2 0];

standardize = true;

distanceMetric = 'euclidean';
% distanceMetric = 'seuclidean';
% distanceMetric = 'spearman';
% distanceMetric = 'cosine';
% distanceMetric = 'correlation';

if strcmp(distanceMetric,'correlation')
  warning('need to do a Fischer Z-transform before t-test/ANOVA');
  %http://www.mathworks.com/matlabcentral/fileexchange/25367-homogeneity-test-for-multiple-correlation-coefficients/content/fisherz.m
end

parameter = 'trial';
plotit = false;
verbose = false;

sub = 1;
ses = 1;
evVal = 1;

% thisROI = {'center91'};
thisROI = {'center109'};
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
if all(ismember(thisROI,ana.elecGroupsStr))
  elecInd = ismember(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}{1}).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,thisROI)})));
elseif ~all(ismember(thisROI,ana.elecGroupsStr)) && all(ismember(thisROI,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}{1}).sub(sub).data.label))
  elecInd = ismember(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}{1}).sub(sub).data.label,unique(thisROI));
else
  error('Cannot find specified electrode(s)');
end

% simAcross = 'time';
% simAcross = 'chan';

% column numbers in trialinfo
% trialNumCol = 5;
phaseCountCol = 4;
stimNumCol = 6;
categNumCol = 7;
pairNumCol = 13;

% if plotit
if strcmp(distanceMetric,'euclidean')
  distanceScale = [0 100];
elseif strcmp(distanceMetric,'spearman') || strcmp(distanceMetric,'correlation') || strcmp(distanceMetric,'cosine')
  distanceScale = [0 2];
elseif strcmp(distanceMetric,'seuclidean')
  distanceScale = [0 20];
else
  distanceScale = [];
end
% end

%% run it for p1 vs p2 comparing across time

% initialize to store the distance values
D = struct;
for d = 1:length(dataTypes)
  %   fprintf('%s\n',dataTypes{d});
  %   D.(dataTypes{d}).nTrial
  D.(dataTypes{d}).dissim = nan(length(exper.subjects),length(exper.sessions),size(latencies,1),size(latencies,1));
  D.(dataTypes{d}).nTrial = nan(length(exper.subjects),length(exper.sessions));
end

for d = 1:length(dataTypes)
  dataType = dataTypes{d};
  
  fprintf('Processing %s...\n',dataType);
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      if ~exper.badSub(sub,ses)
        
      fprintf('\t%s %s...\n',exper.subjects{sub},exper.sesStr{ses});
      
      subD = nan(size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter),1),size(latencies,1),size(latencies,1));
      
      for lat1 = 1:size(latencies,1)
        if verbose
          fprintf('\t%.2f sec to %.2f sec...\n',latencies(lat1,1),latencies(lat1,2));
        end
        timeS = latencies(lat1,:);
        timeInd1 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.time <= timeS(2) + 0.0001;
        
        for lat2 = 1:size(latencies,1)
          if verbose
            fprintf('\t\t%.2f sec to %.2f sec...\n',latencies(lat2,1),latencies(lat2,2));
          end
          timeS = latencies(lat2,:);
          timeInd2 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.time <= timeS(2) + 0.0001;
          
          if standardize
            nTrl1 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter),1);
            nTrl2 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.(parameter),1);
            data_cat = cat(1, ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd1), ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd2));
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
          else
            dat1 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd1);
            dat2 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd2);
          end
          
          for i = 1:size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter),1)
            p1_trlInd = i;
            p1_phaseCount = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trialinfo(p1_trlInd,phaseCountCol);
            p1_stimNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trialinfo(p1_trlInd,stimNumCol);
            p1_categNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trialinfo(p1_trlInd,categNumCol);
            
            p2_trlInd = find(...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trialinfo(:,stimNumCol) == p1_stimNum & ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trialinfo(:,categNumCol) == p1_categNum);
            
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
              %               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1));
              %               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2));
              %             elseif sum(elecInd) > 1
              %               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1))';
              %               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2))';
              %             end
              %
              %             if plotit
              %               xaxis = linspace(timeS(1),timeS(2),size(p1_data,1));
              %               yaxis = linspace(timeS(1),timeS(2),size(p2_data,1));
              %             end
              %           elseif strcmp(simAcross,'chan')
              %             % rows = channels; cols = samples
              %             p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1));
              %             p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2));
              %
              %             if plotit
              %               xaxis = 1:sum(elecInd);
              %               yaxis = 1:sum(elecInd);
              %
              %               elecLabels_x = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.label(elecInd);
              %               elecLabels_y = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.label(elecInd);
              %             end
              %           end
              
              subD(i,lat1,lat2) = pdist2(p1_data(:)',p2_data(:)',distanceMetric);
              
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
              %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1),2)),'b');
              %     hold on
              %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd2),2)),'r');
              %     %plot(linspace(timeS(1),timeS(2),size(p2_data,1)),p2_data,'r');
              %     hold off
              %     axis([timeS(1) timeS(2) -20 20]);
              
            else
              if verbose
                fprintf('%s: No p2 found for p1 phaseCount=%d stimNum=%d categNum=%d\n',dataType,p1_phaseCount,p1_stimNum,p1_categNum);
              end
            end
          end % trials
          D.(dataType).dissim(sub,ses,lat1,lat2) = nanmean(subD(:,lat1,lat2),1);
        end % lat2
      end % lat1
      D.(dataType).nTrial(sub,ses) = sum(~isnan(subD(:,1,1)));
      end
    end
  end % sub
end % dataTypes

%% stats

% data1_str = 'word_RgH_rc_spac';
% data2_str = 'word_RgH_rc_mass';
% data1_str = 'img_RgH_rc_spac';
% data2_str = 'img_RgH_rc_mass';
% data1_str = 'word_RgH_fo_spac';
% data2_str = 'word_RgH_fo_mass';
% data1_str = 'img_RgH_fo_spac';
% data2_str = 'img_RgH_fo_mass';

comparisons = {...
  {'word_RgH_rc_spac', 'word_RgH_rc_mass'}, ...
  {'img_RgH_rc_spac', 'img_RgH_rc_mass'}, ...
  {'word_RgH_fo_spac', 'word_RgH_fo_mass'}, ...
  {'img_RgH_fo_spac', 'img_RgH_fo_mass'}};

% comparisons = {...
%   {'word_RgH_spac', 'word_RgH_mass'}, ...
%   {'img_RgH_spac', 'img_RgH_mass'}};

% comparisons = {{'img_RgH_spac', 'img_RgH_mass'}};

alpha = 0.05;
tails = 'both';

% trial count threshold - need n or more trials in both comparison conds
nThresh = 5;

for ses = 1:length(exper.sessions)
  for lat1 = 1:size(latencies,1)
    fprintf('\n');
    fprintf('%.2f sec to %.2f sec...\n',latencies(lat1,1),latencies(lat1,2));
    for lat2 = 1:size(latencies,1)
      fprintf('\t%.2f sec to %.2f sec...\n',latencies(lat2,1),latencies(lat2,2));
      
      for cmp = 1:length(comparisons)
        data1_str = comparisons{cmp}{1};
        data2_str = comparisons{cmp}{2};
        
        threshSub = D.(data1_str).nTrial(:,ses) >= nThresh & D.(data2_str).nTrial(:,ses) >= nThresh;
        
        data1 = D.(data1_str).dissim(:,ses,lat1,lat2);
        data1 = data1(threshSub);
        data2 = D.(data2_str).dissim(:,ses,lat1,lat2);
        data2 = data2(threshSub);
        
        if strcmp(distanceMetric,'correlation')
          % need to do a Fischer Z-transform before t-test/ANOVA
          %
          % http://www.mathworks.com/matlabcentral/fileexchange/25367-homogeneity-test-for-multiple-correlation-coefficients/content/fisherz.m
          data1 = data1(:);
          data1 = 0.5 .* log((1+data1) ./ (1-data1));
          data2 = data2(:);
          data2 = 0.5 .* log((1+data2) ./ (1-data2));
        end
        
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
      end % cmp
    end % lat2
  end % lat2
end % ses

%% run it for p1 vs p2 comparing across stim types (p1 image, p2 word)

% looking for evidence of the image (defined by p1) during word p2

dataTypes = {{'img_RgH_rc_spac', 'word_RgH_rc_spac'}, {'img_RgH_rc_mass', 'word_RgH_rc_mass'}, ...
  {'img_RgH_fo_spac', 'word_RgH_fo_spac'}, {'img_RgH_fo_mass', 'word_RgH_fo_mass'}};

% initialize to store the distance values
D = struct;
for d = 1:length(dataTypes)
  %   fprintf('%s\n',dataTypes{d});
  %   D.(dataTypes{d}).nTrial
  datstr = dataTypes{d}{1};
  if length(dataTypes{d}) > 1
    for di = 2:length(dataTypes{d})
      datstr = cat(2,datstr,'_',dataTypes{d}{di});
    end
  end
  D.(datstr).dissim = nan(length(exper.subjects),length(exper.sessions),size(latencies,1));
  D.(datstr).nTrial = nan(length(exper.subjects),length(exper.sessions));
end

for d = 1:length(dataTypes)
  datstr = dataTypes{d}{1};
  if length(dataTypes{d}) > 1
    for di = 2:length(dataTypes{d})
      datstr = cat(2,datstr,'_',dataTypes{d}{di});
    end
  end
  
  dataType = dataTypes{d};
  
  fprintf('Processing %s vs %s...\n',dataType{1},dataType{2});
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      if ~exper.badSub(sub,ses)
      fprintf('\t%s %s...\n',exper.subjects{sub},exper.sesStr{ses});
      
      subD = nan(size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter),1),size(latencies,1));
      
      for lat1 = 1:size(latencies,1)
        if verbose
          fprintf('\t%.2f sec to %.2f sec...\n',latencies(lat1,1),latencies(lat1,2));
        end
        timeS = latencies(lat1,:);
        timeInd1 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.time <= timeS(2) + 0.0001;
          timeInd2 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.time <= timeS(2) + 0.0001;
          
          if standardize
            nTrl1 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter),1);
            nTrl2 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.(parameter),1);
            data_cat = cat(1, ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter)(:,elecInd,timeInd1), ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.(parameter)(:,elecInd,timeInd2));
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
          else
            dat1 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter)(:,elecInd,timeInd1);
            dat2 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.(parameter)(:,elecInd,timeInd2);
          end
          
          for i = 1:size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter),1)
            p1_trlInd = i;
            %p1_trialNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,trialNumCol);
            p1_phaseCount = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,phaseCountCol);
            p1_stimNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,stimNumCol);
            p1_categNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,categNumCol);
            p1_pairNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,pairNumCol);
            
%             p2_trlInd = find(...
%               data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
%               data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,stimNumCol) == p1_stimNum & ...
%               data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,categNumCol) == p1_categNum);
            
            p2_trlInd = find(...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,pairNumCol) == p1_pairNum & ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,categNumCol) ~= p1_categNum);
            
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
              %               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1));
              %               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2));
              %             elseif sum(elecInd) > 1
              %               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1))';
              %               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2))';
              %             end
              %
              %             if plotit
              %               xaxis = linspace(timeS(1),timeS(2),size(p1_data,1));
              %               yaxis = linspace(timeS(1),timeS(2),size(p2_data,1));
              %             end
              %           elseif strcmp(simAcross,'chan')
              %             % rows = channels; cols = samples
              %             p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1));
              %             p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2));
              %
              %             if plotit
              %               xaxis = 1:sum(elecInd);
              %               yaxis = 1:sum(elecInd);
              %
              %               elecLabels_x = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.label(elecInd);
              %               elecLabels_y = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.label(elecInd);
              %             end
              %           end
              
              subD(i,lat1) = pdist2(p1_data(:)',p2_data(:)',distanceMetric);
              
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
                
                title(sprintf('%s vs %s (%d vs %d): phaseCount=%d stimNum=%d categNum=%d\n',strrep(dataType{1},'_','-'),strrep(dataType{2},'_','-'),p1_trlInd,p2_trlInd,p1_phaseCount,p1_stimNum,p1_categNum));
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
              %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1),2)),'b');
              %     hold on
              %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd2),2)),'r');
              %     %plot(linspace(timeS(1),timeS(2),size(p2_data,1)),p2_data,'r');
              %     hold off
              %     axis([timeS(1) timeS(2) -20 20]);
              
            else
              if verbose
                fprintf('%s vs %s: No p2 found for p1 phaseCount=%d stimNum=%d categNum=%d\n',dataType{1},dataType{2},p1_phaseCount,p1_stimNum,p1_categNum);
              end
            end
          end % trials
          D.(datstr).dissim(sub,ses,lat1) = nanmean(subD(:,lat1));
      end % lat1
      D.(datstr).nTrial(sub,ses) = sum(~isnan(subD(:,1)));
      end
    end
  end % sub
end % dataTypes

%% RMANOVA - p1 vs p2 comparing stim types

cnds = {'img_RgH_rc_spac_word_RgH_rc_spac', 'img_RgH_fo_spac_word_RgH_fo_spac', 'img_RgH_rc_mass_word_RgH_rc_mass', 'img_RgH_fo_mass_word_RgH_fo_mass'};

nThresh = 1;

goodSub = ones(length(exper.subjects),1);

for ses = 1:length(exper.sesStr)
  for cnd = 1:length(cnds)
    goodSub = goodSub .* D.(cnds{cnd}).nTrial(:,ses) >= nThresh;
  end
end


anovaData = [];

for sub = 1:length(exper.subjects)
  
  %if ~exper.badSub(sub)
  if goodSub(sub)
    for ses = 1:length(exper.sesStr)
      
      theseData = [];
      
      for cnd = 1:length(cnds)
        
        for t = 1:size(D.(cnds{cnd}).dissim,3)
          
          theseData = cat(2,theseData,D.(cnds{cnd}).dissim(sub,ses,t));
          
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
  
end

varnames = {'spacing','subseqMem','time'};
O = teg_repeated_measures_ANOVA(anovaData, [2 2 4], varnames);

%% run it for p1 vs p2 comparing across time, comparing across stim types (p1 image, p2 word)

dataTypes = {{'img_RgH_rc_spac', 'word_RgH_rc_spac'}, {'img_RgH_rc_mass', 'word_RgH_rc_mass'}, ...
  {'img_RgH_fo_spac', 'word_RgH_fo_spac'}, {'img_RgH_fo_mass', 'word_RgH_fo_mass'}};

% initialize to store the distance values
D = struct;
for d = 1:length(dataTypes)
  %   fprintf('%s\n',dataTypes{d});
  %   D.(dataTypes{d}).nTrial
  datstr = dataTypes{d}{1};
  if length(dataTypes{d}) > 1
    for di = 2:length(dataTypes{d})
      datstr = cat(2,datstr,'_',dataTypes{d}{di});
    end
  end
  D.(datstr).dissim = nan(length(exper.subjects),length(exper.sessions),size(latencies,1),size(latencies,1));
  D.(datstr).nTrial = nan(length(exper.subjects),length(exper.sessions));
end

for d = 1:length(dataTypes)
  datstr = dataTypes{d}{1};
  if length(dataTypes{d}) > 1
    for di = 2:length(dataTypes{d})
      datstr = cat(2,datstr,'_',dataTypes{d}{di});
    end
  end
  
  dataType = dataTypes{d};
  
  fprintf('Processing %s vs %s...\n',dataType{1},dataType{2});
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      if ~exper.badSub(sub,ses)
      fprintf('\t%s %s...\n',exper.subjects{sub},exper.sesStr{ses});
      
      subD = nan(size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter),1),size(latencies,1),size(latencies,1));
      
      for lat1 = 1:size(latencies,1)
        if verbose
          fprintf('\t%.2f sec to %.2f sec...\n',latencies(lat1,1),latencies(lat1,2));
        end
        timeS = latencies(lat1,:);
        timeInd1 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.time <= timeS(2) + 0.0001;
        
        for lat2 = 1:size(latencies,1)
          if verbose
            fprintf('\t\t%.2f sec to %.2f sec...\n',latencies(lat2,1),latencies(lat2,2));
          end
          timeS = latencies(lat2,:);
          timeInd2 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.time <= timeS(2) + 0.0001;
          
          if standardize
            nTrl1 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter),1);
            nTrl2 = size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.(parameter),1);
            data_cat = cat(1, ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter)(:,elecInd,timeInd1), ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.(parameter)(:,elecInd,timeInd2));
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
          else
            dat1 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter)(:,elecInd,timeInd1);
            dat2 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.(parameter)(:,elecInd,timeInd2);
          end
          
          for i = 1:size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.(parameter),1)
            p1_trlInd = i;
            %p1_trialNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,trialNumCol);
            p1_phaseCount = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,phaseCountCol);
            p1_stimNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,stimNumCol);
            p1_categNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,categNumCol);
            p1_pairNum = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType{1})).sub(sub).data.trialinfo(p1_trlInd,pairNumCol);
            
%             p2_trlInd = find(...
%               data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
%               data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,stimNumCol) == p1_stimNum & ...
%               data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,categNumCol) == p1_categNum);
            
            p2_trlInd = find(...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,phaseCountCol) == p1_phaseCount & ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,pairNumCol) == p1_pairNum & ...
              data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType{2})).sub(sub).data.trialinfo(:,categNumCol) ~= p1_categNum);
            
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
              %               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1));
              %               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2));
              %             elseif sum(elecInd) > 1
              %               p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1))';
              %               p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2))';
              %             end
              %
              %             if plotit
              %               xaxis = linspace(timeS(1),timeS(2),size(p1_data,1));
              %               yaxis = linspace(timeS(1),timeS(2),size(p2_data,1));
              %             end
              %           elseif strcmp(simAcross,'chan')
              %             % rows = channels; cols = samples
              %             p1_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1));
              %             p2_data = squeeze(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p2_trlInd,elecInd,timeInd2));
              %
              %             if plotit
              %               xaxis = 1:sum(elecInd);
              %               yaxis = 1:sum(elecInd);
              %
              %               elecLabels_x = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.label(elecInd);
              %               elecLabels_y = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.label(elecInd);
              %             end
              %           end
              
              subD(i,lat1,lat2) = pdist2(p1_data(:)',p2_data(:)',distanceMetric);
              
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
                
                title(sprintf('%s vs %s (%d vs %d): phaseCount=%d stimNum=%d categNum=%d\n',strrep(dataType{1},'_','-'),strrep(dataType{2},'_','-'),p1_trlInd,p2_trlInd,p1_phaseCount,p1_stimNum,p1_categNum));
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
              %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd1),2)),'b');
              %     hold on
              %     plot(linspace(timeS(1),timeS(2),size(p1_data,1)),squeeze(mean(data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.trial(p1_trlInd,elecInd,timeInd2),2)),'r');
              %     %plot(linspace(timeS(1),timeS(2),size(p2_data,1)),p2_data,'r');
              %     hold off
              %     axis([timeS(1) timeS(2) -20 20]);
              
            else
              if verbose
                fprintf('%s vs %s: No p2 found for p1 phaseCount=%d stimNum=%d categNum=%d\n',dataType{1},dataType{2},p1_phaseCount,p1_stimNum,p1_categNum);
              end
            end
          end % trials
          D.(datstr).dissim(sub,ses,lat1,lat2) = nanmean(subD(:,lat1,lat2),1);
        end % lat2
      end % lat1
      D.(datstr).nTrial(sub,ses) = sum(~isnan(subD(:,1,1)));
      end
    end
  end % sub
end % dataTypes

%% stats

% data1_str = 'word_RgH_rc_spac';
% data2_str = 'word_RgH_rc_mass';
% data1_str = 'img_RgH_rc_spac';
% data2_str = 'img_RgH_rc_mass';
% data1_str = 'word_RgH_fo_spac';
% data2_str = 'word_RgH_fo_mass';
% data1_str = 'img_RgH_fo_spac';
% data2_str = 'img_RgH_fo_mass';

% comparisons = {...
%   {'word_RgH_rc_spac', 'word_RgH_rc_mass'}, ...
%   {'img_RgH_rc_spac', 'img_RgH_rc_mass'}, ...
%   {'word_RgH_fo_spac', 'word_RgH_fo_mass'}, ...
%   {'img_RgH_fo_spac', 'img_RgH_fo_mass'}};

% comparisons = {...
%   {'word_RgH_spac', 'word_RgH_mass'}, ...
%   {'img_RgH_spac', 'img_RgH_mass'}};

% comparisons = {{'img_RgH_spac_word_RgH_spac','img_RgH_mass_word_RgH_mass'}};
comparisons = {{'img_RgH_rc_spac_word_RgH_rc_spac','img_RgH_rc_mass_word_RgH_rc_mass'},{'img_RgH_fo_spac_word_RgH_fo_spac','img_RgH_fo_mass_word_RgH_fo_mass'}};

% comparisons = {{'img_RgH_spac', 'img_RgH_mass'}};

alpha = 0.05;
tails = 'both';

% trial count threshold - need n or more trials in both comparison conds
nThresh = 5;

for ses = 1:length(exper.sessions)
  for lat1 = 1:size(latencies,1)
    fprintf('\n');
    fprintf('%.2f sec to %.2f sec...\n',latencies(lat1,1),latencies(lat1,2));
    for lat2 = 1:size(latencies,1)
      fprintf('\t%.2f sec to %.2f sec...\n',latencies(lat2,1),latencies(lat2,2));
      
      for cmp = 1:length(comparisons)
        data1_str = comparisons{cmp}{1};
        data2_str = comparisons{cmp}{2};
        
        threshSub = D.(data1_str).nTrial(:,ses) >= nThresh & D.(data2_str).nTrial(:,ses) >= nThresh;
        
        data1 = D.(data1_str).dissim(:,ses,lat1,lat2);
        data1 = data1(threshSub);
        data2 = D.(data2_str).dissim(:,ses,lat1,lat2);
        data2 = data2(threshSub);
        
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
      end % cmp
    end % lat2
  end % lat2
end % ses

%% run it for p1 vs p2

ses = 1;
evVal = 1;

% initialize to store the distance values
D = struct;
for d = 1:length(dataTypes)
%   fprintf('%s\n',dataTypes{d});
%   D.(dataTypes{d}).nTrial
  D.(dataTypes{d}).dissim = nan(length(exper.subjects),length(exper.sessions),size(latencies,1));
  D.(dataTypes{d}).nTrial = nan(length(exper.subjects),length(exper.sessions));
end

for d = 1:length(dataTypes)
  dataType = dataTypes{d};
  
  fprintf('Processing %s...\n',dataType);
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      if ~exper.badSub(sub,ses)
      fprintf('\t%s %s...\n',exper.subjects{sub},exper.sesStr{ses});
      
      subD = nan(size(data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter),1),1);
      
      for lat = 1:size(latencies,1)
        if verbose
          fprintf('\t%.2f sec to %.2f sec...\n',latencies(lat,1),latencies(lat,2));
        end
        timeS = latencies(lat,:);
        timeInd = data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}{1}).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}{1}).sub(sub).data.time <= timeS(2) + 0.0001;
        
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
        else
          dat1 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p1',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd);
          dat2 = data_tla.(exper.sesStr{ses}).(sprintf('%s_p2',dataType)).sub(sub).data.(parameter)(:,elecInd,timeInd);
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
        D.(dataType).dissim(sub,ses,lat) = nanmean(subD);
      end % lat
      D.(dataType).nTrial(sub,ses) = sum(~isnan(subD));
      end
    end
  end % sub
end % dataTypes

%% bar plot


% bar(D.img_RgH_rc_spac.dissim(~exper.badSub,1,1));

measure = 'recall_hr';

figure
hold on

plot(results.oneDay.cued_recall.spaced.recall.(measure)(~exper.badSub)',D.img_RgH_rc_spac.dissim(~exper.badSub,1,1),'ko');
plot(results.oneDay.cued_recall.spaced.recall.(measure)(~exper.badSub)',D.img_RgH_rc_mass.dissim(~exper.badSub,1,1),'rx');

hold off

%% correlation

corrType = 'Pearson';
% corrType = 'Spearman';

measure = 'recall_hr';
% measure = 'recall_rt';
% measure = 'recall_rt_hit';

% distcases = {'img_RgH_rc_spac', 'img_RgH_rc_mass', 'word_RgH_rc_spac', 'word_RgH_rc_mass', ...
%   'img_RgH_fo_spac', 'img_RgH_fo_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass'};

distcases = {'word_RgH_rc_spac', 'word_RgH_rc_mass'};

for t = 1:size(D.(distcases{1}).dissim,3)
  fprintf('\n');
  latStr = sprintf('%.2f s to %.2f s',latencies(t,1),latencies(t,2));
  fprintf('%s...\n',latStr);
  for d = 1:length(distcases)
    if ~isempty(strfind(distcases{d},'spac'))
      cond = 'spaced';
    elseif ~isempty(strfind(distcases{d},'mass'))
      cond = 'massed';
    end
    
    %x = results.oneDay.cued_recall.(cond).recall.(measure)(~exper.badSub);
    %y = D.(distcases{d}).dissim(~exper.badSub,1,t);
    
    x = D.(distcases{d}).dissim(~exper.badSub,1,t);
    y = results.oneDay.cued_recall.(cond).recall.(measure)(~exper.badSub);
    
    stats = regstats(x,y,'linear','cookd');
    outliers = stats.cookd > 3*mean(stats.cookd);
    
    [rho ,p] = corr(x(~outliers),y(~outliers),'type',corrType);
    fprintf('\t%s vs %s (%d),\trho=%.3f, p=%.5f\n',measure,distcases{d},sum(~outliers),rho,p);
    
    if p < .05
      figure
      hold on
      
      plot(x(~outliers),y(~outliers),'ko');
      plot(x(outliers),y(outliers),'rx');
      
      hold off
      
      title(sprintf('%s, %s: rho=%.3f, p=%.5f',latStr,strrep(distcases{d},'_','-'),rho,p));
      %xlabel(strrep(measure,'_','-'));
      %ylabel('distance');
      
      xlabel('distance');
      ylabel(sprintf('%s %s',cond,strrep(measure,'_','-')));
      %ylim([75 110]);
    end
  end
end


% mdl = fitlm(x,y);
% 
% plotDiagnostics(mdl,'cookd')

%% RM ANOVA - p1 vs p2: spac x stimType x memory x time

cnds = {'word_RgH_rc_spac', 'img_RgH_rc_spac', 'word_RgH_fo_spac', 'img_RgH_fo_spac', 'word_RgH_rc_mass', 'img_RgH_rc_mass', 'word_RgH_fo_mass', 'img_RgH_fo_mass'};

nThresh = 1;

goodSub = ones(length(exper.subjects),1);

for ses = 1:length(exper.sesStr)
  for cnd = 1:length(cnds)
    goodSub = goodSub .* D.(cnds{cnd}).nTrial(:,ses) >= nThresh;
  end
end


anovaData = [];

for sub = 1:length(exper.subjects)
  
  if ~exper.badSub(sub)
  %if goodSub(sub)
    for ses = 1:length(exper.sesStr)
      
      theseData = [];
      
      for cnd = 1:length(cnds)
        
        for t = 1:size(D.(cnds{cnd}).dissim,3)
          
          theseData = cat(2,theseData,D.(cnds{cnd}).dissim(sub,ses,t));
          
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
  
end

varnames = {'spacing','subseqMem','stimType','time'};
O = teg_repeated_measures_ANOVA(anovaData, [2 2 2 4], varnames);

%% RM ANOVA - p1 vs p2: spac x memory x time

cnds = {'word_RgH_rc_spac', 'word_RgH_fo_spac', 'word_RgH_rc_mass', 'word_RgH_fo_mass'};

nThresh = 1;

goodSub = ones(length(exper.subjects),1);

for ses = 1:length(exper.sesStr)
  for cnd = 1:length(cnds)
    goodSub = goodSub .* D.(cnds{cnd}).nTrial(:,ses) >= nThresh;
  end
end


anovaData = [];

for sub = 1:length(exper.subjects)
  
  if ~exper.badSub(sub)
  %if goodSub(sub)
    for ses = 1:length(exper.sesStr)
      
      theseData = [];
      
      for cnd = 1:length(cnds)
        
        for t = 1:size(D.(cnds{cnd}).dissim,3)
          
          theseData = cat(2,theseData,D.(cnds{cnd}).dissim(sub,ses,t));
          
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
  
end

varnames = {'spacing','subseqMem','time'};
O = teg_repeated_measures_ANOVA(anovaData, [2 2 size(D.(cnds{cnd}).dissim,3)], varnames);

%% RM ANOVA - p1 vs p2: spac x time (trying with SCD)

cnds = {'word_RgH_rc_spac', 'word_RgH_rc_mass'};

nThresh = 1;

goodSub = ones(length(exper.subjects),1);

for ses = 1:length(exper.sesStr)
  for cnd = 1:length(cnds)
    goodSub = goodSub .* D.(cnds{cnd}).nTrial(:,ses) >= nThresh;
  end
end


anovaData = [];

for sub = 1:length(exper.subjects)
  
  if ~exper.badSub(sub)
  %if goodSub(sub)
    for ses = 1:length(exper.sesStr)
      
      theseData = [];
      
      for cnd = 1:length(cnds)
        
        for t = 1:size(D.(cnds{cnd}).dissim,3)
          
          theseData = cat(2,theseData,D.(cnds{cnd}).dissim(sub,ses,t));
          
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
  
end

varnames = {'spacing','time'};
O = teg_repeated_measures_ANOVA(anovaData, [2 4], varnames);

%% plot RSA spacing x subsequent memory interaction

ses=1;

% rc_spaced_mean = mean([mean(D.word_RgH_rc_spac.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_rc_spac.dissim(~exper.badSub,ses,:),3)],2);
% rc_massed_mean = mean([mean(D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_rc_mass.dissim(~exper.badSub,ses,:),3)],2);
% fo_spaced_mean = mean([mean(D.word_RgH_fo_spac.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_fo_spac.dissim(~exper.badSub,ses,:),3)],2);
% fo_massed_mean = mean([mean(D.word_RgH_fo_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_fo_mass.dissim(~exper.badSub,ses,:),3)],2);

rc_spaced_mean = mean(D.word_RgH_rc_spac.dissim(~exper.badSub,ses,:),3);
rc_massed_mean = mean(D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:),3);
fo_spaced_mean = mean(D.word_RgH_fo_spac.dissim(~exper.badSub,ses,:),3);
fo_massed_mean = mean(D.word_RgH_fo_mass.dissim(~exper.badSub,ses,:),3);

figure
hold on

rc_mark = 'ko';
fo_mark = 'rx';

% recalled
plot(0.9*ones(sum(~exper.badSub(:,ses)),1), rc_spaced_mean,rc_mark,'LineWidth',1);
plot(1.9*ones(sum(~exper.badSub(:,ses)),1), rc_massed_mean,rc_mark,'LineWidth',1);

% forgotten
plot(1.1*ones(sum(~exper.badSub(:,ses)),1), fo_spaced_mean,fo_mark,'LineWidth',1);
plot(2.1*ones(sum(~exper.badSub(:,ses)),1), fo_massed_mean,fo_mark,'LineWidth',1);

meanSizeR = 15;
meanSizeF = 20;

% recalled
hr = plot(1, mean(rc_spaced_mean),rc_mark,'LineWidth',3,'MarkerSize',meanSizeR);
plot(2, mean(rc_massed_mean),rc_mark,'LineWidth',3,'MarkerSize',meanSizeR);

% forgotten
hf = plot(1, mean(fo_spaced_mean),fo_mark,'LineWidth',3,'MarkerSize',meanSizeF);
plot(2, mean(fo_massed_mean),fo_mark,'LineWidth',3,'MarkerSize',meanSizeF);

hold off
axis square
% axis([0.75 2.25 75 110]);
xlim([0.75 2.25]);

set(gca,'XTick', [1 2]);

set(gca,'XTickLabel',{'Spaced','Massed'});
ylabel('Euclidean Distance');

title('RSA: Spacing \times Subsequent Memory');
legend([hr, hf],{'Recalled','Forgotten'},'Location','North');

publishfig(gcf,0,[],[],[]);

print(gcf,'-depsc2','~/Desktop/rsa_spacXsm.eps');

%% plot RSA spacing x time interaction

ses=1;

spaced_mean = squeeze(mean(cat(2,D.word_RgH_rc_spac.dissim(~exper.badSub,ses,:), D.img_RgH_rc_spac.dissim(~exper.badSub,ses,:), ...
  D.word_RgH_fo_spac.dissim(~exper.badSub,ses,:), D.img_RgH_fo_spac.dissim(~exper.badSub,ses,:)),2));

% spaced_mean = squeeze(D.word_RgH_rc_spac.dissim(~exper.badSub,ses,:));

spaced_sem = std(spaced_mean,0,1) / sqrt(sum(~exper.badSub));

massed_mean = squeeze(mean(cat(2,D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:), D.img_RgH_rc_mass.dissim(~exper.badSub,ses,:), ...
  D.word_RgH_fo_mass.dissim(~exper.badSub,ses,:), D.img_RgH_fo_mass.dissim(~exper.badSub,ses,:)),2));

% massed_mean = squeeze(D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:));

massed_sem = std(massed_mean,0,1) / sqrt(sum(~exper.badSub));

% rc_massed_mean = mean([mean(D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_rc_mass.dissim(~exper.badSub,ses,:),3)],2);
% fo_massed_mean = mean([mean(D.word_RgH_fo_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_fo_mass.dissim(~exper.badSub,ses,:),3)],2);

sp_mark = 'ks-';
ma_mark = 'ro-';

meanSize_spac = 15;
meanSize_mass = 15;

figure
hold on

% plot(mean(spaced_mean,1)',sp_mark,'LineWidth',2,'MarkerSize',meanSizeS);
% plot(mean(massed_mean,1)',ma_mark,'LineWidth',2,'MarkerSize',meanSizeM);
hs = errorbar(mean(spaced_mean,1)',spaced_sem,sp_mark,'LineWidth',2,'MarkerSize',meanSize_spac,'MarkerFaceColor','k');
hm = errorbar(mean(massed_mean,1)',massed_sem,ma_mark,'LineWidth',2,'MarkerSize',meanSize_mass,'MarkerFaceColor','r');

hold off
axis square
% axis([0.75 4.25 97 103]);
% axis([0.75 4.25 95 105]);
xlim([0.75 4.25]);

set(gca,'XTick', [1 2 3 4]);

set(gca,'XTickLabel',{'100-300','300-500','500-700','700-900'});
xlabel('Time (ms)');
ylabel('Euclidean Distance');

title('RSA: Spacing \times Time');
legend([hm,hs],{'Massed','Spaced'},'Location','SouthWest');

publishfig(gcf,0,[],[],[]);

% print(gcf,'-depsc2','~/Desktop/rsa_spacXtime.eps');

%% plot RSA spacing x stimType x time interaction

ses=1;

spaced_w_mean = squeeze(mean(cat(2,D.word_RgH_rc_spac.dissim(~exper.badSub,ses,:), ...
  D.word_RgH_fo_spac.dissim(~exper.badSub,ses,:)),2));

spaced_w_sem = std(spaced_w_mean,0,1) / sqrt(sum(~exper.badSub));

spaced_i_mean = squeeze(mean(cat(2,D.img_RgH_rc_spac.dissim(~exper.badSub,ses,:), ...
  D.img_RgH_fo_spac.dissim(~exper.badSub,ses,:)),2));

spaced_i_sem = std(spaced_i_mean,0,1) / sqrt(sum(~exper.badSub));

massed_w_mean = squeeze(mean(cat(2,D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:), ...
  D.word_RgH_fo_mass.dissim(~exper.badSub,ses,:)),2));

massed_w_sem = std(massed_w_mean,0,1) / sqrt(sum(~exper.badSub));

massed_i_mean = squeeze(mean(cat(2,D.img_RgH_rc_mass.dissim(~exper.badSub,ses,:), ...
  D.img_RgH_fo_mass.dissim(~exper.badSub,ses,:)),2));

massed_i_sem = std(massed_i_mean,0,1) / sqrt(sum(~exper.badSub));


% rc_massed_mean = mean([mean(D.word_RgH_rc_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_rc_mass.dissim(~exper.badSub,ses,:),3)],2);
% fo_massed_mean = mean([mean(D.word_RgH_fo_mass.dissim(~exper.badSub,ses,:),3) mean(D.img_RgH_fo_mass.dissim(~exper.badSub,ses,:),3)],2);

spac_w_mark = 'ks-';
spac_i_mark = 'ko--';
mass_w_mark = 'rs-';
mass_i_mark = 'ro--';

meanSize_spac = 15;
meanSize_mass = 15;

figure
hold on

% plot(mean(spaced_mean,1)',sp_mark,'LineWidth',2,'MarkerSize',meanSizeS);
% plot(mean(massed_mean,1)',ma_mark,'LineWidth',2,'MarkerSize',meanSizeM);
hsw = errorbar(mean(spaced_w_mean,1)',spaced_w_sem,spac_w_mark,'LineWidth',2,'MarkerSize',meanSize_spac,'MarkerFaceColor','k');
hsi = errorbar(mean(spaced_i_mean,1)',spaced_i_sem,spac_i_mark,'LineWidth',2,'MarkerSize',meanSize_spac);
hmw = errorbar(mean(massed_w_mean,1)',massed_w_sem,mass_w_mark,'LineWidth',2,'MarkerSize',meanSize_mass,'MarkerFaceColor','r');
hmi = errorbar(mean(massed_i_mean,1)',massed_i_sem,mass_i_mark,'LineWidth',2,'MarkerSize',meanSize_mass);

hold off
axis square
%axis([0.75 4.25 95 105]);
xlim([0.75 4.25]);

set(gca,'XTick', [1 2 3 4]);

set(gca,'XTickLabel',{'100-300','300-500','500-700','700-900'});
xlabel('Time (ms)');
ylabel('Euclidean Distance');

title('RSA: Spacing \times Stim Type \times Time');
legend([hmw, hmi, hsw, hsi],{'Massed Word','Massed Image','Spaced Word','Spaced Image'},'Location','SouthWest');

publishfig(gcf,0,[],[],[]);

% print(gcf,'-depsc2','~/Desktop/rsa_spacXstimXtime.eps');

%% stats - ttest

% data1_str = 'word_RgH_rc_spac';
% data2_str = 'word_RgH_rc_mass';
% data1_str = 'img_RgH_rc_spac';
% data2_str = 'img_RgH_rc_mass';
% data1_str = 'word_RgH_fo_spac';
% data2_str = 'word_RgH_fo_mass';
% data1_str = 'img_RgH_fo_spac';
% data2_str = 'img_RgH_fo_mass';

% comparisons = {...
%   {'word_RgH_rc_spac', 'word_RgH_rc_mass'}, ...
%   {'img_RgH_rc_spac', 'img_RgH_rc_mass'}, ...
%   {'word_RgH_fo_spac', 'word_RgH_fo_mass'}, ...
%   {'img_RgH_fo_spac', 'img_RgH_fo_mass'}};

comparisons = {...
  {'word_RgH_rc_spac', 'word_RgH_rc_mass'}, ...
  {'word_RgH_fo_spac', 'word_RgH_fo_mass'}};

% comparisons = {...
%   {'word_RgH_rc_spac', 'word_RgH_rc_mass'}};

% comparisons = {...
%   {'word_RgH_spac', 'word_RgH_mass'}, ...
%   {'img_RgH_spac', 'img_RgH_mass'}};

% comparisons = {{'img_RgH_spac', 'img_RgH_mass'}};

alpha = 0.05;
tails = 'both';

% trial count threshold - need n or more trials in both comparison conds
nThresh = 1;

for ses = 1:length(exper.sessions)
  for lat = 1:size(latencies,1)
    fprintf('\n');
    fprintf('%.2f sec to %.2f sec...\n',latencies(lat,1),latencies(lat,2));
    
    for cmp = 1:length(comparisons)
      data1_str = comparisons{cmp}{1};
      data2_str = comparisons{cmp}{2};
      
      threshSub = D.(data1_str).nTrial(:,ses) >= nThresh & D.(data2_str).nTrial(:,ses) >= nThresh;
      
      data1 = D.(data1_str).dissim(:,ses,lat);
      data1 = data1(threshSub);
      data2 = D.(data2_str).dissim(:,ses,lat);
      data2 = data2(threshSub);
      
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
end

%% plot

conds = {'word_RgH_rc_spac', 'word_RgH_rc_mass', ...
  'img_RgH_rc_spac', 'img_RgH_rc_mass', ...
  'word_RgH_fo_spac', 'word_RgH_fo_mass', ...
  'img_RgH_fo_spac', 'img_RgH_fo_mass'};

plotStyles = {'ro-','ko-','r^-','k^-','ro--','ko--','r^--','k^--'};
% red = spaced
% black = massed
% solid = recalled
% dashed = forgotten
% circles = words
% triangles = images
figure;
for i = 1:length(conds)
  plot(squeeze(nanmean(D.(conds{i}).dissim,1)),plotStyles{i},'LineWidth',1.2);
  hold on
end
hold off
legend(strrep(conds,'_','-'),'Location','SouthWest');
axis([0.5 4.5 85 100]);
% ylim([85 100]);


%% run it for face vs house

ses = 1;
evVal = 1;

dissim = struct;

for sub = 1:length(exper.subjects)
  dissim.sub(sub).D = [];
  
  fprintf('Subject %s...\n',exper.subjects{sub});
  
  nTrl1 = size(data_tla.(exper.sesStr{ses}).(dataTypes{1}).sub(sub).data.(parameter),1);
  nTrl2 = size(data_tla.(exper.sesStr{ses}).(dataTypes{2}).sub(sub).data.(parameter),1);
  
  subD = nan(nTrl1+nTrl2,nTrl1+nTrl2);
  
  if ~exper.badSub(sub,ses)
  for lat = 1:size(latencies,1)
    fprintf('\t%.2f sec to %.2f sec...\n',latencies(lat,1),latencies(lat,2));
    timeS = latencies(lat,:);
    timeInd = data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}{1}).sub(sub).data.time >= timeS(1) - 0.0001 & data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}{1}).sub(sub).data.time <= timeS(2) + 0.0001;
    
    if standardize
      data_cat = cat(1, ...
        data_tla.(exper.sesStr{ses}).(dataTypes{1}).sub(sub).data.(parameter)(:,elecInd,timeInd), ...
        data_tla.(exper.sesStr{ses}).(dataTypes{2}).sub(sub).data.(parameter)(:,elecInd,timeInd));
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
    else
      dat1 = data_tla.(exper.sesStr{ses}).(dataTypes{1}).sub(sub).data.(parameter)(:,elecInd,timeInd);
      dat2 = data_tla.(exper.sesStr{ses}).(dataTypes{2}).sub(sub).data.(parameter)(:,elecInd,timeInd);
    end
    
    d1count = 0;
    for d = 1:length(dataTypes)
      dataType = dataTypes{d};
      
      fprintf('\t\tProcessing %s...\n',dataTypes{d});
      
      for i = 1:size(data_tla.(exper.sesStr{ses}).(dataTypes{d}).sub(sub).data.(parameter),1)
        d1count = d1count + 1;
        %dt1_data = squeeze(dat1(i,:,:));
        dt1_data = squeeze(dat_z_ravel(d1count,:,:));
        
        d2count = 0;
        for d2 = 1:length(dataTypes)
          for j = 1:size(data_tla.(exper.sesStr{ses}).(dataTypes{d2}).sub(sub).data.(parameter),1)
            d2count = d2count + 1;
            
            %dt2_data = squeeze(dat2(j,:,:));
            dt2_data = squeeze(dat_z_ravel(d2count,:,:));
            
            subD(d1count,d2count) = pdist2(dt1_data(:)',dt2_data(:)',distanceMetric);
            
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
            
            %         else
            %           if verbose
            %             fprintf('%s: No p2 found for p1 phaseCount=%d stimNum=%d categNum=%d\n',dataType,p1_phaseCount,p1_stimNum,p1_categNum);
            %           end
            %         end
          end
        end
        
      end % trials
      %D.(dataType).dissim(sub,lat) = nanmean(subD);
      %D.(dataType).nTrial(sub) = sum(~isnan(subD));
    end % dataTypes
    dissim.sub(sub).D = cat(3,dissim.sub(sub).D,subD);
  end % lat
  end
end % sub

%% plot it

% distanceScale = [0 150];
% distanceScale = [0 100];
% distanceScale = [0 75];

distanceScale = [0 1.5];

for sub = 1:length(exper.subjects)
  for lat = 1:size(latencies,1)
    figure;
    imagesc(dissim.sub(sub).D(:,:,lat),distanceScale);
    axis square;
    title(sprintf('\t%.2f sec to %.2f sec\n',latencies(lat,1),latencies(lat,2)));
  end
end
