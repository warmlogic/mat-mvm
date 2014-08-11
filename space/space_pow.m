%% new analysis details loading method

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
% procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');
procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');

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

replaceDatatype = {'tla','pow'};
[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot,replaceDatatype);

% %% load the analysis details

% adFile = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla/analysisDetails.mat';
% 
% mm_splitAnalysisDetails(adFile,[],true);

% % adFile = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla/analysisDetails.mat';
% % adFile_other = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla/analysisDetails_9_14.mat';
% % [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails(adFile,adFile_other,true,true,true);
% 
% % % server_adFile = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla/analysisDetails.mat';
% % if exist(server_adFile,'file')
% %   mm_mergeAnalysisDetails(adFile,server_adFile,true,false,false);
% % end
% 
% % [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
% [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,false);

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

ana.freq = mm_freqSet('ndtools');

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
% ana.eventValues = {{'multistudy_word'}};
ana.eventValuesSplit = { ...
  { ...
  { ...
  'img_onePres' ...
  'img_RgH_rc_spac_p1','img_RgH_rc_spac_p2','img_RgH_rc_mass_p1','img_RgH_rc_mass_p2' ...
  'img_RgH_fo_spac_p1','img_RgH_fo_spac_p2','img_RgH_fo_mass_p1','img_RgH_fo_mass_p2' ...
%   %'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
  } ...
  { ...
  'word_onePres' ...
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
    sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr > 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
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
    } ...
    } ...
    };
else
  ana.trl_expr = { ...
    { ...
    { ...
    sprintf('eventNumber == %d & targ == 1 & spaced == 0 & lag == -1 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr == 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 1 & cr_recall_spellCorr < 1 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
%     %sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
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
    } ...
    } ...
    };
end

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

%% use in conjunction with post-loading baseline

% keeptrials = false;
% [data_pow_lsd,exper] = mm_loadSubjectData(exper,dirs,ana,'pow',keeptrials,'trialinfo',true);
% 
% % overwrite ana.eventValues with the new split events
% ana.eventValues = ana.eventValuesSplit;

%% new loading workflow - pow

cfg = [];

% cfg.loadMethod = 'seg';
cfg.loadMethod = 'trialinfo';
cfg.latency = 'all';
cfg.frequency = 'all';

cfg.keeptrials = 'no';
%cfg.keeptrials = 'yes';
% cfg.equatetrials = 'no';
% %cfg.equatetrials = 'yes';

cfg.rmPreviousCfg = true;

% type of input (used in the filename to load)
cfg.ftype = 'pow';
% cfg.ftype = 'fourier';

% type of output: 'pow', 'coh', 'phase'
cfg.output = 'pow';

% transformation: 'log10', 'log', 'vec'
cfg.transform = '';
% cfg.transform = 'log10';
% cfg.transform = 'vec';

% normalization of single or average trials
% cfg.norm_trials = 'single'; % Grandchamp & Delorme (2011)
cfg.norm_trials = 'average';

% baseline type
% % 'zscore', 'absolute', 'relchange', 'relative', 'db'
cfg.baseline_type = 'zscore';
% cfg.baseline_type = 'db';
% cfg.baseline_type = 'absolute';
% cfg.baseline_type = 'relchange';
% cfg.baseline_type = 'relative';

% % baseline using all events
% cfg.baseline_events = 'all';

% baseline only using word (since image immedpately follows word)
cfg.baseline_events = { ...
  'word_onePres' ...
  'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
  'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
  %'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
  };

% baseline period
%cfg.baseline_time = [-0.2 0];
% cfg.baseline_time = [-0.3 0];
cfg.baseline_time = [-0.3 -0.1];

% at what data stage should it be baseline corrected?
cfg.baseline_data = 'pow';
% mod is not an option
% % cfg.baseline_data = 'mod';

%cfg.saveFile = true;
cfg.saveFile = false;

% only keep induced data by removing evoked?
cfg.rmevoked = 'no';
cfg.rmevokedfourier = 'no';
cfg.rmevokedpow = 'no';
% cfg.rmevoked = 'yes';
% cfg.rmevokedfourier = 'yes';
% cfg.rmevokedpow = 'no';
% if strcmp(cfg.rmevoked,'yes') && ~exist('data_evoked','var')
%   % load('/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0_art_zeroVar/tla_-1000_2000_avg/data_evoked.mat');
%   
%   % local testing
%   %load('/Users/matt/data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar_badChanManual_badChanEP/tla_-1000_2000_avg/data_evoked.mat');
% end

% if isfield(cfg,'equatetrials') && strcmp(cfg.equatetrials,'yes')
%   eq_str = '_eq';
% else
%   eq_str = '';
% end
% if isfield(cfg,'keeptrials') && strcmp(cfg.keeptrials,'yes')
%   kt_str = '_trials';
% else
%   kt_str = '_avg';
% end
% if isfield(cfg,'rmevoked') && strcmp(cfg.rmevoked,'yes')
%   indu_str = '_induced';
% else
%   indu_str = '_whole';
% end
% saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s%s%s%s.mat',cfg.output,eq_str,kt_str,indu_str));

[data_pow,exper] = mm_ft_loadData_multiSes2(cfg,exper,dirs,ana);

% if exist(saveFile,'file')
%   fprintf('Loading saved file: %s\n',saveFile);
%   load(saveFile);
% else
%   fprintf('Running mm_ft_loadData_multiSes\n');
%   if exist('data_evoked','var')
%     [data_pow,exper] = mm_ft_loadData_multiSes(cfg,exper,dirs,ana,data_evoked);
%   else
%     [data_pow,exper] = mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
%   end
%   if cfg.saveFile
%     fprintf('Saving %s...\n',saveFile);
%     save(saveFile,sprintf('data_%s',cfg.output),'exper','cfg');
%   end
% end
fprintf('Done.\n');

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

% %% OLD - load in the subject data
% 
% % % make sure ana.eventValues is set properly
% % if ~iscell(ana.eventValues{1})
% %   ana.eventValues = {ana.eventValues};
% % end
% % if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
% %   ana.eventValues = {exper.eventValues};
% % end
% 
% % keeptrials = true;
% % % [data_tla,exper] = mm_ft_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');
% % [data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');
% 
% keeptrials = false;
% [data_pow,exper] = mm_loadSubjectData(exper,dirs,ana,'pow',keeptrials,'trialinfo');
% 
% % %% get rid of the bad channels
% % 
% % cfg = [];
% % cfg.printRoi = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% % [data_tla] = mm_rmBadChan(cfg,exper,ana,data_tla);
% 
% % overwrite ana.eventValues with the new split events
% ana.eventValues = ana.eventValuesSplit;

%% save

% saveDir = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
% saveDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
saveDir = dirs.saveDirProc;
% save(fullfile(saveDir,'space_word_data_ga_pow.mat'),'data_pow','ga_pow','exper','ana','dirs','files','-v7.3');
% save(fullfile(saveDir,'space_word_img_data_ga_pow.mat'),'data_pow','ga_pow','exper','ana','dirs','files','-v7.3');

% save(fullfile(saveDir,'space_word_img_data_pow.mat'),'data_pow','exper','ana','dirs','files','-v7.3');
save(fullfile(saveDir,'space_word_img_data_pow_wordBL.mat'),'data_pow','exper','ana','dirs','files','-v7.3');
% clear data_pow

%% load

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
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
else
  error('Data directory not found.');
end

loadDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/pow');
% loadDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow');

% load(fullfile(loadDir,'space_word_img_data_pow.mat'));
load(fullfile(loadDir,'space_word_img_data_pow_wordBL.mat'));

[dirs] = mm_checkDirs(dirs);

ana.freq = mm_freqSet('ndtools');

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {{}};
% exper.badBehSub = {{'SPACE001','SPACE008','SPACE011','SPACE017','SPACE019','SPACE039'}};
% exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE039'}};
exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE030'}};

% SPACE001 - low trial counts
% SPACE008 - did not do task correctly

% SPACE011 - noisy TF data

evToCheck = { ...
  { ...
  { ...
  'img_onePres' ...
  'img_RgH_rc_spac_p2','img_RgH_rc_mass_p2','img_RgH_fo_spac_p2','img_RgH_fo_mass_p2' ...
  %'img_RgH_rc_spac_p1','img_RgH_rc_mass_p1','img_RgH_fo_spac_p1','img_RgH_fo_mass_p1' ...
  %'img_RgM_spac_p1,'img_RgM_mass_p1'','img_RgM_spac_p2','img_RgM_mass_p2' ...
  } ...
  { ...
  'word_onePres' ...
  'word_RgH_rc_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_spac_p2','word_RgH_fo_mass_p2' ...
  %'word_RgH_rc_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_spac_p1','word_RgH_fo_mass_p1' ...
  %'word_RgM_spac_p1','word_RgM_mass_p1','word_RgM_spac_p2','word_RgM_mass_p2' ...
  } ...
  } ...
  };

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,10,[],'vert',evToCheck);

%% Test plots to make sure data look ok

cfg_ft = [];
% cfg_ft.baseline = [-0.3 -0.1];
% cfg_ft.baselinetype = 'absolute';
% if strcmp(cfg_ft.baselinetype,'absolute')
%   %cfg_ft.zlim = [-400 400];
%   cfg_ft.zlim = [-2 2];
% elseif strcmp(cfg_ft.baselinetype,'relative')
%   cfg_ft.zlim = [0 2.0];
% end
% cfg_ft.parameter = 'powspctrm';
% cfg_ft.ylim = [3 9];
% cfg_ft.ylim = [3 14];
cfg_ft.ylim = [3 32];
% cfg_ft.showlabels = 'yes';
% cfg_ft.fontsize = 12;
% cfg_ft.colorbar = 'yes';
% cfg_ft.interactive = 'yes';
% cfg_ft.layout = ft_prepare_layout([],ana);
% sub=6;
% ses=1;
% for i = 1:length(ana.eventValues{1})
%   figure
%   ft_multiplotTFR(cfg_ft,data_freq.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
%   title(ana.eventValues{1}{i});
% end

% cfg_ft.channel = {'E124'};
% cfg_ft.channel = {'E117'};
cfg_ft.channel = {'E55'};

% cfg_ft.baseline = [-0.4 -0.1];
% cfg_ft.baselinetype = 'absolute';
% % if strcmp(cfg_ft.baselinetype,'absolute')
% %   %cfg_ft.zlim = [-2 2];
% %   cfg_ft.zlim = [-500 500];
% % elseif strcmp(cfg_ft.baselinetype,'relative')
% %   cfg_ft.zlim = [0 1.5];
% % end
% cfg_ft.showlabels = 'yes';
% cfg_ft.colorbar = 'yes';
% cfg_ft.ylim = [4 8];

%cfg_ft.zlim = [-150 150];
%cfg_ft.zlim = [-300 300];

% cfg_ft.zlim = [-.6 .6];
%cfg_ft.zlim = [-.30 .30];
% cfg_ft.zlim = [-.15 .15];

%cfg_ft.zlim = [-2 2];
% cfg_ft.zlim = [-1 1];
cfg_ft.zlim = [-3 3];


% cfg_ft.showlabels = 'yes';
% cfg_ft.fontsize = 12;
% cfg_ft.colorbar = 'yes';
% cfg_ft.layout = ft_prepare_layout([],ana);

sub = 24;
ses = 1;
typ = 1;
for evVal = 1:length(ana.eventValues{ses}{typ})
  figure
  ft_singleplotTFR(cfg_ft,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
  set(gcf,'Name',sprintf('%s',strrep(ana.eventValues{ses}{typ}{evVal},'_','-')))
  title(strrep(ana.eventValues{ses}{typ}{evVal},'_','-'));
  
%   figure
%   ft_multiplotTFR(cfg_ft,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
%   set(gcf,'Name',sprintf('%s',strrep(ana.eventValues{ses}{typ}{evVal},'_','-')))
%   title(strrep(ana.eventValues{ses}{typ}{evVal},'_','-'));
end

%% baseline

% % log power
% data_pow_log = data_pow;
% for sub = 1:length(exper.subjects)
%   for cnd = 1:length(conds)
%     data_pow_log.(exper.sesStr{ses}).(conds{cnd}).sub(sub).data.powspctrm = log10(data_pow.(exper.sesStr{ses}).(conds{cnd}).sub(sub).data.powspctrm);
%   end
% end

% % fieldtrip's baseline method
% cfg_bl = [];
% cfg_bl.baseline = [-0.3 -0.1];
% cfg_bl.baselinetype = 'db';
% for ses = 1:length(exper.sesStr)
%   for typ = 1:length(ana.eventValues{ses})
%     for evVal = 1:length(ana.eventValues{ses}{typ})
%       for sub = 1:length(exper.subjects)
%         data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_freqbaseline(cfg_bl,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
%       end
%     end
%   end
% end

% % my zscore method - within an event value
% cfg = [];
% cfg.baseline_time = [-0.3 -0.1];
% cfg.baseline_type = 'zscore';
% cfg.param = 'powspctrm';
% for ses = 1:length(exper.sesStr)
%   for typ = 1:length(ana.eventValues{ses})
%     for evVal = 1:length(ana.eventValues{ses}{typ})
%       for sub = 1:length(exper.subjects)
%         subSesEvData = data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data;
%         
%         % get the baseline time indices
%         blt = subSesEvData.time >= cfg.baseline_time(1) & subSesEvData.time <= cfg.baseline_time(2);
%         % mean across baseline period
%         
%         % individual trial
%         %blm = nanmean(subSesEvData.(cfg.param)(:,:,:,blt),4);
%         % average
%         blm = nanmean(subSesEvData.(cfg.param)(:,:,blt),3);
%         
%         if strcmp(cfg.baseline_type,'zscore')
%           fprintf('Z-transforming %s relative to mean([%.2f %.2f] pre-stimulus).\n',ana.eventValues{ses}{typ}{evVal},cfg.baseline_time(1),cfg.baseline_time(2));
%           % std across time, then avg across events (lower freqs often get smaller std)
%           
%           % individual trial
%           %blstd = nanmean(nanstd(subSesEvData.(cell2mat(fn)).(cfg.param)(:,:,:,blt),0,4),1);
%           % average
%           blstd = nanstd(subSesEvData.(cfg.param)(:,:,blt),0,3);
%           
%           % % avg across time, then std across events (higher freqs often get smaller std)
%           % blstd = nanstd(nanmean(subSesEvData.(cfg.param)(:,:,:,blt),4),0,1);
%           
%           % % concatenate all times of all events and std (equivalent std across freqs)
%           %blstd = shiftdim(squeeze(std(double(reshape(shiftdim(subSesEvData.(cfg.param)(:,:,:,blt),3),...
%           %  size(subSesEvData.(cfg.param)(:,:,:,blt),1)*size(subSesEvData.(cfg.param)(:,:,:,blt),4),...
%           %  size(subSesEvData.(cfg.param)(:,:,:,blt),2),size(subSesEvData.(cfg.param)(:,:,:,blt),3))),0,1)),-1);
%           
%           data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(cfg.param),blm),blstd);
%         end
%       end
%     end
%   end
% end

% my zscore method - combine all events first
cfg = [];
cfg.baseline_time = [-0.3 -0.1];
cfg.baseline_type = 'zscore';
% cfg.baseline_type = 'db';
cfg.param = 'powspctrm';
for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    allData = [];
    % get the baseline time indices
    blt = false(size(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{1}{1}).sub(sub).data.time));
    blt_tbeg = nearest(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{1}{1}).sub(sub).data.time,cfg.baseline_time(1));
    blt_tend = nearest(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{1}{1}).sub(sub).data.time,cfg.baseline_time(2));
    blt(blt_tbeg:blt_tend) = true;
    
    for typ = 1:length(ana.eventValues{ses})
      for evVal = 1:length(ana.eventValues{ses}{typ})
        allData = cat(4,allData,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param)(:,:,blt));
      end
    end
    
    % average across time, then average all events together
    blm = nanmean(nanmean(allData,3),4);
    
    if strcmp(cfg.baseline_type,'zscore')
      % std across time, then avg across events (lower freqs often get smaller std)
      blstd = nanmean(nanstd(allData,0,3),4);
    end
    
    % % avg across time, then std across events (higher freqs often get smaller std)
    % blstd = nanstd(nanmean(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param)(:,:,:,blt),4),0,1);
    
    % % concatenate all times of all events and std (equivalent std across freqs)
    %blstd = shiftdim(squeeze(std(double(reshape(shiftdim(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param)(:,:,:,blt),3),...
    %  size(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param)(:,:,:,blt),1)*size(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param)(:,:,:,blt),4),...
    %  size(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param)(:,:,:,blt),2),size(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param)(:,:,:,blt),3))),0,1)),-1);
    for typ = 1:length(ana.eventValues{ses})
      for evVal = 1:length(ana.eventValues{ses}{typ})
        if strcmp(cfg.baseline_type,'zscore')
          fprintf('Z-transforming %s relative to mean([%.2f %.2f] pre-stimulus).\n',ana.eventValues{ses}{typ}{evVal},cfg.baseline_time(1),cfg.baseline_time(2));
          data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param) = bsxfun(@rdivide,bsxfun(@minus,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param),blm),blstd);
        elseif strcmp(cfg.baseline_type,'db')
          fprintf('\tConverting to db relative to mean([%.2f %.2f]).\n',cfg.baseline_time(1),cfg.baseline_time(2));
          % divide by baseline mean and convert to db
          data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param) = 10*log10(bsxfun(@rdivide,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(cfg.param),blm));
        end
      end
    end
    
  end
end

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.data_str = 'data_pow';

cfg_ft = [];
cfg_ft.keepindividual = 'no';

if strcmp(cfg_ana.data_str,'data_pow')
  ga_pow = struct;
elseif strcmp(cfg_ana.data_str,'data_pow_log')
  ga_pow_log = struct;
elseif strcmp(cfg_ana.data_str,'data_coh')
  ga_coh = struct;
elseif strcmp(cfg_ana.data_str,'data_evoked')
  ga_evoked = struct;
end

for ses = 1:length(exper.sesStr)
  for typ = 1:length(ana.eventValues{ses})
    cfg_ana.conditions = ana.eventValues{ses}{typ};
    cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
    
    for evVal = 1:length(ana.eventValues{ses}{typ})
      %tic
      fprintf('Running ft_freqgrandaverage on %s...',ana.eventValues{ses}{typ}{evVal});
      if strcmp(cfg_ana.data_str,'data_pow')
        cfg_ft.parameter = 'powspctrm';
        ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
      elseif strcmp(cfg_ana.data_str,'data_pow_log')
        cfg_ft.parameter = 'powspctrm';
        ga_pow_log.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
      elseif strcmp(cfg_ana.data_str,'data_coh')
        %cfg_ft.parameter = 'plvspctrm';
        cfg_ft.parameter = 'powspctrm';
        ga_coh.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
      elseif strcmp(cfg_ana.data_str,'data_evoked')
        cfg_ft.parameter = 'powspctrm';
        ga_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
      end
      fprintf('Done.\n');
      %toc
    end
  end
end

% % turn keeptrial data into average for statistical functions because proper
% % processing of dimord is currently broken
% 
% data_pow_avg = struct;
% 
% cfg_ana = [];
% cfg_ana.is_ga = 0;
% %cfg_ana.conditions = ana.eventValues{ses};
% cfg_ana.data_str = 'data_pow';
% 
% cfg_ft = [];
% cfg_ft.keepindividual = 'yes';
% 
% for ses = 1:length(exper.sesStr)
%   for typ = 1:length(ana.eventValues{ses})
%     cfg_ana.conditions = ana.eventValues{ses}{typ};
%     cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
%     
%     for evVal = 1:length(ana.eventValues{ses}{typ})
%       % collapse subjects together
%       fprintf('%s, %s\n',exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
%       data_pow_avg.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
%       
%       % % keep subjects separate
%       %for sub = 1:length(exper.subjects)
%       %  fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
%       %  data_pow_avg.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_freqgrandaverage(cfg_ft,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
%       %end
%     end
%   end
% end

%% simple plot

chan = 73; % 73 = Pz

% zlim = [-200 400];
% zlim = [0 100];
zlim = [-4 4];

sub = 1;
ses = 1;
typ = 1;

for evVal = 1:length(ana.eventValues{ses}{typ})
  figure;
  
  surf(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.time,...
    data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.freq,...
    squeeze(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.powspctrm(chan,:,:)));
  
  shading interp;view([0,90]);axis tight;
  caxis(zlim);
  %imagesc(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.time,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.freq,squeeze(data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.powspctrm(chan,:,:)),zlim);
  %axis xy;
  set(gcf,'Name',sprintf('%s',strrep(ana.eventValues{ses}{typ}{evVal},'_','-')))
  title(strrep(ana.eventValues{ses}{typ}{evVal},'_','-'));
  colorbar
end

%% plot the conditions - simple

cfg_ft = [];
%cfg_ft.baseline = [-0.3 -0.1];
%cfg_ft.baselinetype = 'absolute';
%if strcmp(cfg_ft.baselinetype,'absolute')
%cfg_ft.xlim = [0 1];
%cfg_ft.ylim = [4 64];
%cfg_ft.ylim = [4 8];
%cfg_ft.ylim = [8 12];
%cfg_ft.ylim = [12 28];
%cfg_ft.ylim = [28 64];
%cfg_ft.zlim = [-1 1];
%cfg_ft.zlim = [-2 2];
%elseif strcmp(cfg_ft.baselinetype,'relative')
cfg_ft.zlim = [-3 3];
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
for ses = 1:length(ana.eventValues)
  for typ = 1:length(ana.eventValues{ses})
    for evVal = 1:length(ana.eventValues{ses}{typ})
      figure
      ft_multiplotTFR(cfg_ft,ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}));
      set(gcf,'Name',sprintf('%s',strrep(ana.eventValues{ses}{typ}{evVal},'_','-')))
    end
  end
end

%% subplots of each subject's power spectrum

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'FS'},{'PS'}};
% cfg_plot.rois = {{'PS'}};
% cfg_plot.rois = {{'RPI2'}};
%cfg_plot.rois = {'E124'};
%cfg_plot.rois = {'E25'};
%cfg_plot.rois = {'RAS'};
cfg_plot.rois = {{'RAS2'}};
%cfg_plot.rois = {'LPS','RPS'};
%cfg_plot.rois = {'LPS'};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
ses=1;
% cfg_plot.condByROI = repmat(ana.eventValues{ses},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{{'word_RgH_rc_spac_p2','word_RgH_rc_mass_p2'}}}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{{'word_RgH_fo_spac_p2'}}}},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{{{'img_RgH_rc_spac_p2'}}}},size(cfg_plot.rois));

cfg_plot.types = {''};


cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.zlim = [-2 2];
cfg_ft.parameter = 'powspctrm';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotTFR(cfg_ft,cfg_plot,ana,exper,data_pow,ses);
end

%% make some GA plots

files.saveFigs = false;

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
cfg_ft.xlim = [0 1.0]; % time
% cfg_ft.ylim = [3 8]; % freq
% cfg_ft.ylim = [8 12]; % freq
% cfg_ft.ylim = [8 10]; % freq
% cfg_ft.ylim = [10 12]; % freq
cfg_ft.ylim = [3 28]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.zlim = [-2 2]; % pow

cfg_ft.parameter = 'powspctrm';

cfg_plot = [];
cfg_plot.plotTitle = 1;

%cfg_plot.rois = {{'FS'},{'LAS','RAS'},{'LPS','RPS'}};
%cfg_plot.rois = {{'FS'},{'PS'}};
%cfg_plot.rois = {'E71'};
cfg_plot.rois = {'all'};

cfg_plot.is_ga = 1;
% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'word_onePres', 'word_RgH_spac_p2', 'word_RgH_mass_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'word_onePres', 'word_RgH_spac_p1', 'word_RgH_mass_p1', 'word_RgH_spac_p2', 'word_RgH_mass_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{{'word_RgH_spac_p1', 'word_RgH_mass_p1', 'word_RgH_spac_p2', 'word_RgH_mass_p2'}}}},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{{{'word_RgH_rc_spac_p2', 'word_RgH_rc_mass_p2', 'word_RgH_fo_spac_p2', 'word_RgH_fo_mass_p2'}}}},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{{{'img_RgH_rc_spac_p2', 'img_RgH_rc_mass_p2', 'img_RgH_fo_spac_p2', 'img_RgH_fo_mass_p2'}}}},size(cfg_plot.rois));

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

%cfg_plot.ftFxn = 'ft_singleplotTFR';

% cfg_plot.ftFxn = 'ft_topoplotTFR';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.comment = 'no';
% %cfg_ft.xlim = [0.5 0.8]; % time
% cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time

cfg_plot.ftFxn = 'ft_multiplotTFR';
cfg_ft.showlabels = 'yes';
cfg_ft.comment = '';

ses = 1;
for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_plotTFR(cfg_ft,cfg_plot,ana,exper,files,dirs,ga_pow,ses);
end

%% line plots - no cluster stats

files.saveFigs = true;
files.figPrintFormat = 'png';

cfg = [];
cfg.parameter = 'powspctrm';

%cfg.times = [-0.2:0.05:0.9; -0.1:0.05:1.0]';
%cfg.times = [-0.2:0.1:0.9; -0.1:0.1:1.0]';
% cfg.times = [-0.2:0.2:0.8; 0:0.2:1.0]';

cfg.times = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
% cfg.times = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap

% cfg.freqs = [3 8; 8 12; 12 28; 28 50; 50 80];
% cfg.freqs = [3 8; 8 10; 10 12];
% cfg.freqs = [3 7; 8 12; 8 10; 11 12; 13 20; 23 30; 32 47; 51 80];
% cfg.freqs = [3 7; 8 12; 13 20; 23 30; 32 47; 51 80];
% cfg.freqs = [4.1 7.7; 8.4 10.1; 11 12; 8.4 12; 13.1 20.5; 22.4 29.2; 31.9 49.7; 54.3 77.4];

cfg.freqs = [ ...
  ana.freq.theta; ...
  ana.freq.alpha; ...
  ana.freq.alpha_lower; ...
  ana.freq.alpha_upper; ...
  ana.freq.beta_lower; ...
  ana.freq.beta_upper; ...
  ana.freq.gamma_lower; ...
  ana.freq.gamma_upper];

% cfg.rois = {...
%   {'LAS'},{'FS'},{'RAS'},...
%   {'LPS'},{'PS'},{'RPS'},...
%   };

% cfg.rois = {...
%   {'LAI'},{'FI'},{'RAI'},...
%   {'LAS'},{'FS'},{'RAS'},...
%   {'LPS'},{'PS'},{'RPS'},...
%   {'LPI'},{'PI'},{'RPI'},...
%   };

cfg.rois = {...
  {'LAS2'},{'FC'},{'RAS2'},...
  {'LT'},{'C'},{'RT'},...
  {'LPS'},{'PS'},{'RPS'},...
  {'LPI2'},{'PI'},{'RPI2'},...
  {'Oz'}
  };

% ses=1;
% cfg.conditions = ana.eventValues{ses};

cfg.conditions = {{'word_onePres','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'}};
cfg.conditions = {{'img_onePres','img_RgH_rc_spac_p2','img_RgH_fo_spac_p2','img_RgH_rc_mass_p2','img_RgH_fo_mass_p2'}};

cfg.plotTitle = true;
cfg.plotLegend = true;

cfg.plotErrorBars = false;
cfg.eb_transp = true;

cfg.plotClusSig = false;
% cfg.clusAlpha = 0.1;
% %cfg.clusTimes = cfg.times;
% % cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
% %cfg.clusTimes = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
% cfg.clusTimes = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
% cfg.clusLimits = false;

cfg.linewidth = 2;
% cfg.limitlinewidth = 0.5;
% cfg.textFontSize = 10;

%cfg.ylim = [-0.6 0.6];
%cfg.ylim = [-0.5 0.2];
cfg.nCol = 3;

% % whole power
% cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
mm_ft_lineTFR(cfg,ana,exper,files,dirs,ga_pow);


% mm_ft_clusterplotTFR

% nk_ft_avgpowerbytime - see cosi2_ft_seg_pow line 1158

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
cfg_ana.rois = {{'PS'},{'FS'},{'LPS'},{'RPS'},{'PS'},{'PS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.2 0.4; 0.6 1.0; 0.5 0.8; 0.5 0.8; 0.5 0.8; 0.5 1.0];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [4 8; 4 8; 4 8; 4 8; 4 8; 8 12];

%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};
cfg_ana.conditions = {'all'};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.avgoverfreq = 'yes';
cfg_ft.parameter = 'powspctrm';
cfg_ft.correctm = 'fdr';

% line plot parameters
cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
cfg_plot.ylims = repmat([-1 1],size(cfg_ana.rois'));
%cfg_plot.ylims = repmat([-100 100],size(cfg_ana.rois'));

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_freq);
end

%% cluster statistics

cfg_ft = [];
cfg_ft.avgoverchan = 'no';
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverfreq = 'yes';
%cfg_ft.avgoverfreq = 'no';

cfg_ft.parameter = 'powspctrm';

% debugging
%cfg_ft.numrandomization = 100;

cfg_ft.numrandomization = 500;
cfg_ft.clusteralpha = .05;
cfg_ft.alpha = .025;

cfg_ana = [];
cfg_ana.roi = 'all';
%cfg_ana.roi = 'center74';
cfg_ana.avgFrq = cfg_ft.avgoverfreq;
cfg_ana.avgTime = cfg_ft.avgovertime;
%cfg_ana.conditions = {'all'};

% cfg_ana.conditions = {{'recall','no_recall'}};

% cfg_ana.conditions = {{'word_RgH_spac_p1', 'word_RgH_mass_p1'}, {'word_RgH_spac_p2', 'word_RgH_mass_p2'}, ...
%   {'img_RgH_spac_p1', 'img_RgH_mass_p1'}, {'img_RgH_spac_p2', 'img_RgH_mass_p2'}};

cfg_ana.conditions = {...
  {'img_RgH_rc_spac_p1', 'img_RgH_fo_spac_p1'} ... % Spac P1 SME
  {'img_RgH_rc_spac_p2', 'img_RgH_fo_spac_p2'} ... % Spac P2 SME
  {'img_RgH_rc_mass_p1', 'img_RgH_fo_mass_p1'} ... % Mass P1 SME
  {'img_RgH_rc_mass_p2', 'img_RgH_fo_mass_p2'} ... % Mass P2 SME
  {'img_RgH_rc_spac_p1', 'img_RgH_rc_spac_p2'} ... % Spac Rc Repetition
  {'img_RgH_fo_spac_p1', 'img_RgH_fo_spac_p2'} ... % Spac Fo Repetition
  {'img_RgH_rc_mass_p1', 'img_RgH_rc_mass_p2'} ... % Mass Rc Repetition
  {'img_RgH_fo_mass_p1', 'img_RgH_fo_mass_p2'} ... % Mass Fo Repetition
  {'img_RgH_rc_spac_p1', 'img_RgH_rc_mass_p1'} ... % P1 Rc Spacing
  {'img_RgH_rc_spac_p2', 'img_RgH_rc_mass_p2'} ... % P2 Rc Spacing
  {'img_RgH_fo_spac_p1', 'img_RgH_fo_mass_p1'} ... % P1 Fo Spacing
  {'img_RgH_fo_spac_p2', 'img_RgH_fo_mass_p2'} ... % P2 Fo Spacing
  {'img_RgH_rc_spac_p1', 'img_RgH_fo_mass_p1'} ... % P1 spac x mem
  {'img_RgH_rc_mass_p1', 'img_RgH_fo_spac_p1'} ... % P1 spac x mem
  {'img_RgH_rc_spac_p2', 'img_RgH_fo_mass_p2'} ... % P2 spac x mem
  {'img_RgH_rc_mass_p2', 'img_RgH_fo_spac_p2'} ... % P2 spac x mem
  {'img_onePres', 'img_RgH_rc_spac_p1'} ...
  {'img_onePres', 'img_RgH_rc_mass_p1'} ...
  {'img_onePres', 'img_RgH_fo_spac_p1'} ...
  {'img_onePres', 'img_RgH_fo_mass_p1'} ...
  {'img_onePres', 'img_RgH_rc_spac_p2'} ...
  {'img_onePres', 'img_RgH_rc_mass_p2'} ...
  {'img_onePres', 'img_RgH_fo_spac_p2'} ...
  {'img_onePres', 'img_RgH_fo_mass_p2'} ...
  {'word_RgH_rc_spac_p1', 'word_RgH_fo_spac_p1'} ... % Spac P1 SME
  {'word_RgH_rc_spac_p2', 'word_RgH_fo_spac_p2'} ... % Spac P2 SME
  {'word_RgH_rc_mass_p1', 'word_RgH_fo_mass_p1'} ... % Mass P1 SME
  {'word_RgH_rc_mass_p2', 'word_RgH_fo_mass_p2'} ... % Mass P2 SME
  {'word_RgH_rc_spac_p1', 'word_RgH_rc_spac_p2'} ... % Spac Rc Repetition
  {'word_RgH_fo_spac_p1', 'word_RgH_fo_spac_p2'} ... % Spac Fo Repetition
  {'word_RgH_rc_mass_p1', 'word_RgH_rc_mass_p2'} ... % Mass Rc Repetition
  {'word_RgH_fo_mass_p1', 'word_RgH_fo_mass_p2'} ... % Mass Fo Repetition
  {'word_RgH_rc_spac_p1', 'word_RgH_rc_mass_p1'} ... % P1 Rc Spacing
  {'word_RgH_rc_spac_p2', 'word_RgH_rc_mass_p2'} ... % P2 Rc Spacing
  {'word_RgH_fo_spac_p1', 'word_RgH_fo_mass_p1'} ... % P1 Fo Spacing
  {'word_RgH_fo_spac_p2', 'word_RgH_fo_mass_p2'} ... % P2 Fo Spacing
  {'word_RgH_rc_spac_p1', 'word_RgH_fo_mass_p1'} ... % P1 spac x mem
  {'word_RgH_rc_mass_p1', 'word_RgH_fo_spac_p1'} ... % P1 spac x mem
  {'word_RgH_rc_spac_p2', 'word_RgH_fo_mass_p2'} ... % P2 spac x mem
  {'word_RgH_rc_mass_p2', 'word_RgH_fo_spac_p2'} ... % P2 spac x mem
  {'word_onePres', 'word_RgH_rc_spac_p1'} ...
  {'word_onePres', 'word_RgH_rc_mass_p1'} ...
  {'word_onePres', 'word_RgH_fo_spac_p1'} ...
  {'word_onePres', 'word_RgH_fo_mass_p1'} ...
  {'word_onePres', 'word_RgH_rc_spac_p2'} ...
  {'word_onePres', 'word_RgH_rc_mass_p2'} ...
  {'word_onePres', 'word_RgH_fo_spac_p2'} ...
  {'word_onePres', 'word_RgH_fo_mass_p2'} ...
};

cfg_ana.dirStr = '';

if strcmp(cfg_ft.avgovertime,'no')
  %cfg_ana.latencies = [0 0.5];
  %cfg_ana.latencies = [0 1.0];
  %cfg_ana.latencies = [0 0.5; 0.5 1.0];
%   cfg_ana.latencies = [0.02 0.46; 0.5 0.98];
  cfg_ana.latencies = [0.02 0.5; 0.52 1.0];
elseif strcmp(cfg_ft.avgovertime,'yes')
  %cfg_ana.latencies = [-0.2:0.1:0.9; -0.1:0.1:1.0]'; % some overlap
  %cfg_ana.latencies = [-0.2:0.2:0.8; 0:0.2:1.0]'; % some overlap
  
  %cfg_ana.latencies = [-0.18:0.12:0.92; -0.1:0.12:1.0]'; % 100 no overlap
  %cfg_ana.latencies = [-0.18:0.2:0.82; -0.02:0.2:1.0]'; % 200 no overlap
  
  % new analyses
  % cfg_ana.latencies = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
  % cfg_ana.latencies = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
  % cfg_ana.latencies = [0.02:0.25:0.92; 0.25:0.25:1.0]'; % 250 no overlap
  % cfg_ana.latencies = [0.02:0.33:0.92; 0.33:0.33:1.0]'; % 330 no overlap
  cfg_ana.latencies = [0.02:0.5:0.92; 0.5:0.5:1.0]'; % 500 no overlap
  
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end

if strcmp(cfg_ft.avgoverfreq,'no')
  % cfg_ana.frequencies = [4 100];
  cfg_ana.frequencies = ana.freq.fullRange;
elseif strcmp(cfg_ft.avgoverfreq,'yes')
  % % cfg_ana.frequencies = [4 8; 8.1 14; 14.1 21; 21.1 28; 28.1 42; 42.1 64; 64.1 80];
  % % cfg_ana.frequencies = [4 8; 8 12; 12 28; 28 50; 50 100];
  %   cfg_ana.frequencies = [3 8; 9 12; 9 10; 11 12; 13 28; 29 50; 51 80];
  %   cfg_ana.frequencies = [3 7; 8 12; 8 10; 11 12; 13 20; 21 30; 31 45; 46 80]; % hanslmayr
  %   cfg_ana.frequencies = [3 7; 8 12; 8 10; 11 12; 13 20; 23 30; 32 47; 51 80];
  
  cfg_ana.frequencies = [ ...
  ana.freq.theta; ...
  ana.freq.alpha; ...
  ana.freq.alpha_lower; ...
  ana.freq.alpha_upper; ...
  ana.freq.beta_lower; ...
  ana.freq.beta_upper; ...
  ana.freq.gamma_lower; ...
  ana.freq.gamma_upper];

%   cfg_ana.frequencies = [8 10; 10 12];
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgF'];
end
%cfg_ana.latencies = [0 0.5; 0.5 1.0; 1.0 1.615];

if strcmp(cfg_ft.avgoverchan,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgC'];
end

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  for fr = 1:size(cfg_ana.frequencies,1)
    cfg_ft.frequency = cfg_ana.frequencies(fr,:);
    
    [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_pow);
  end
end

%% plot the cluster statistics

files.saveFigs = 1;
files.figPrintFormat = 'png';

cfg_ft = [];
%cfg_ft.alpha = .025;
cfg_ft.alpha = .05;
% cfg_ft.alpha = .1;
cfg_ft.avgoverfreq = cfg_ana.avgFrq;

cfg_plot = [];
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.frequencies = cfg_ana.frequencies;
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.dirStr = cfg_ana.dirStr;

if strcmp(cfg_ft.avgoverfreq,'no')
  % not averaging over frequencies - only works with ft_multiplotTFR
  files.saveFigs = 0;
  cfg_ft.interactive = 'yes';
  cfg_plot.mask = 'yes';
  %cfg_ft.maskstyle = 'saturation';
  cfg_ft.maskalpha = 0.3;
  cfg_plot.ftFxn = 'ft_multiplotTFR';
  % http://mailman.science.ru.nl/pipermail/fieldtrip/2009-July/002288.html
  % http://mailman.science.ru.nl/pipermail/fieldtrip/2010-November/003312.html
end

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  for fr = 1:size(cfg_plot.frequencies,1)
    cfg_ft.frequency = cfg_plot.frequencies(fr,:);
    
    mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs);
  end
end

%% line plots

files.saveFigs = 1;
files.figPrintFormat = 'png';

cfg = [];
cfg.parameter = 'powspctrm';

%cfg.times = [-0.2:0.05:0.9; -0.1:0.05:1.0]';
%cfg.times = [-0.2:0.1:0.9; -0.1:0.1:1.0]';
% cfg.times = [-0.2:0.2:0.8; 0:0.2:1.0]';
cfg.times = cfg_ana.latencies;

% cfg.freqs = [4 8; 8 12; 12 28; 28 50; 50 100];
%cfg.freqs = [4 8];
cfg.freqs = cfg_ana.frequencies;

% cfg.rois = {...
%   {'LAS'},{'FS'},{'RAS'},...
%   {'LPS'},{'PS'},{'RPS'},...
%   };

cfg.rois = {...
  {'LAS2'},{'FC'},{'RAS2'},...
  {'LT'},{'C'},{'RT'},...
  {'LPS'},{'PS'},{'RPS'},...
  {'LPI2'},{'PI'},{'RPI2'},...
  {'Oz'}
  };

% % frontal
% cfg.rois = {...
%   {'LAS2'},{'FC'},{'RAS2'},...
%   {'LT'},{'C'},{'RT'},...
%   {'Pz'}
%   };

% % posterior
% cfg.rois = {...
%   {'LPS'},{'PS'},{'RPS'},...
%   {'LPI2'},{'PI'},{'RPI2'},...
%   {'Oz'}
%   };

% cfg.conditions = ana.eventValues;
% cfg.conditions = cfg_ana.conditions;

% cfg.conditions = {...
%   {'word_onePres', 'word_RgH_rc_spac_p1', 'word_RgH_rc_mass_p1'} ...
%   {'word_onePres', 'word_RgH_fo_spac_p1', 'word_RgH_fo_mass_p1'} ...
%   {'word_onePres', 'word_RgH_rc_spac_p1', 'word_RgH_fo_spac_p1'} ...
%   {'word_onePres', 'word_RgH_rc_mass_p1', 'word_RgH_fo_mass_p1'} ...
%   {'word_onePres', 'word_RgH_rc_spac_p2', 'word_RgH_rc_mass_p2'} ...
%   {'word_onePres', 'word_RgH_fo_spac_p2', 'word_RgH_fo_mass_p2'} ...
%   {'word_onePres', 'word_RgH_rc_spac_p2', 'word_RgH_fo_spac_p2'} ...
%   {'word_onePres', 'word_RgH_rc_mass_p2', 'word_RgH_fo_mass_p2'} ...
%   };

cfg.conditions = {...
  {'img_onePres', 'img_RgH_rc_spac_p1', 'img_RgH_rc_mass_p1'} ...
  {'img_onePres', 'img_RgH_fo_spac_p1', 'img_RgH_fo_mass_p1'} ...
  {'img_onePres', 'img_RgH_rc_spac_p1', 'img_RgH_fo_spac_p1'} ...
  {'img_onePres', 'img_RgH_rc_mass_p1', 'img_RgH_fo_mass_p1'} ...
  {'img_onePres', 'img_RgH_rc_spac_p2', 'img_RgH_rc_mass_p2'} ...
  {'img_onePres', 'img_RgH_fo_spac_p2', 'img_RgH_fo_mass_p2'} ...
  {'img_onePres', 'img_RgH_rc_spac_p2', 'img_RgH_fo_spac_p2'} ...
  {'img_onePres', 'img_RgH_rc_mass_p2', 'img_RgH_fo_mass_p2'} ...
  {'word_onePres', 'word_RgH_rc_spac_p1', 'word_RgH_rc_mass_p1'} ...
  {'word_onePres', 'word_RgH_fo_spac_p1', 'word_RgH_fo_mass_p1'} ...
  {'word_onePres', 'word_RgH_rc_spac_p1', 'word_RgH_fo_spac_p1'} ...
  {'word_onePres', 'word_RgH_rc_mass_p1', 'word_RgH_fo_mass_p1'} ...
  {'word_onePres', 'word_RgH_rc_spac_p2', 'word_RgH_rc_mass_p2'} ...
  {'word_onePres', 'word_RgH_fo_spac_p2', 'word_RgH_fo_mass_p2'} ...
  {'word_onePres', 'word_RgH_rc_spac_p2', 'word_RgH_fo_spac_p2'} ...
  {'word_onePres', 'word_RgH_rc_mass_p2', 'word_RgH_fo_mass_p2'} ...
  };

%   {'img_onePres', 'img_RgH_rc_spac_p2'} ...
%   {'img_onePres', 'img_RgH_rc_mass_p2'} ...
%   {'img_onePres', 'img_RgH_fo_spac_p2'} ...
%   {'img_onePres', 'img_RgH_fo_mass_p2'} ...
%   {'word_RgH_rc_spac_p1', 'word_RgH_rc_mass_p1'} ...
%   {'word_RgH_fo_spac_p1', 'word_RgH_fo_mass_p1'} ...
%   {'img_RgH_fo_spac_p1', 'img_RgH_fo_mass_p1'} ...
%   {'img_RgH_rc_spac_p1', 'img_RgH_rc_mass_p1'} ...
%   {'img_RgH_rc_spac_p2', 'img_RgH_rc_mass_p2'} ...
%   {'img_RgH_fo_spac_p2', 'img_RgH_fo_mass_p2'} ...
%   {'img_RgH_rc_spac_p2', 'img_RgH_fo_spac_p2'} ...
%   {'img_RgH_rc_mass_p2', 'img_RgH_fo_mass_p2'} ...

cfg.plotTitle = true;
cfg.plotLegend = true;

cfg.plotErrorBars = false;
cfg.eb_transp = true;

cfg.plotClusSig = true;
% cfg.clusAlpha = 0.1;
cfg.clusAlpha = 0.05;
cfg.clusTimes = cfg_ana.latencies;
% cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
cfg.clusLimits = true;

cfg.linewidth = 2;
cfg.limitlinewidth = 0.5;

cfg.textFontSize = 10;

%cfg.ylim = [-0.6 0.6];
%cfg.ylim = [-0.5 0.2];
cfg.nCol = 3;

% whole power
cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-300_-100';
cfg.clusDirStr = '_avgT_avgF';
cfg.ylabel = 'Z-Trans Pow';
mm_ft_lineTFR(cfg,ana,exper,files,dirs,ga_pow);

% % induced power
% cfg.type = 'line_pow_induced';
% cfg.clusDirStr = '_pow_induced_-300_-100';
% cfg.ylabel = 'Log Pow';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_pow);

% % induced power
% cfg.type = 'line_pow_indu';
% cfg.clusDirStr = '_zpow_indu_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_pow);

% % evoked power
% cfg.type = 'line_pow_evok';
% cfg.clusDirStr = '_zp_evok_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_evoked);

% cfg.type = 'line_coh';
% cfg.clusDirStr = '_coh_-400_-200';
% cfg.ylabel = 'ITC';
% mm_ft_lineTFR(cfg,ana,files,dirs,ga_coh);

%% let me know that it's done
emailme = 1;
if emailme
  subject = sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    };
  send_gmail(subject,mail_message);
end

%% RM ANOVA

% stimType = 'word_';
stimType = 'img_';
memType = 'RgH_';

% spacings = {'mass', 'spac', 'onePres'};
% % oldnew = {'p1'};
% oldnew = {'p2'};
% memConds = {'all'};

% didn't test new words, so can't assess memory, but can use p1
spacings = {'mass', 'spac'};
% spacings = {'spac'};
% oldnew = {'p1', 'p2'};
% oldnew = {'p1'};
oldnew = {'p2'};
memConds = {'rc','fo'};
% memConds = {'rc'};

measure = 'powspctrm';

% latencies = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
% latencies = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap

latencies = [0.02:0.1:0.92; 0.1:0.1:1.0]'; % 100 no overlap
% latencies = [0.02:0.2:0.92; 0.2:0.2:1.0]'; % 200 no overlap
% latencies = [0.02:0.25:0.92; 0.25:0.25:1.0]'; % 250 no overlap
% latencies = [0.02:0.33:0.92; 0.33:0.33:1.0]'; % 330 no overlap
% latencies = [0.02:0.5:0.92; 0.5:0.5:1.0]'; % 500 no overlap

% % theta
freqs = ana.freq.theta;
% latencies = [0.6 1.0]; % word
% latencies = [0.1 0.4]; % img
% latencies = [0.02 0.5; 0.52 1.0];

% roi = {'LAI'}; % yes, neg **
% roi = {'LFP'};
% roi = {'FC'};
% roi = {'RFP'};
% roi = {'RAI'}; % something awry word mass p1 forgot, values are too high

% roi = {{'LAS2','FS'}}; % yes, pos ***
% roi = {'LAS2'}; % yes, pos ***
% roi = {'LAS'}; % yes, pos **
% roi = {'FS'}; % yes, pos **
% roi = {'C'}; % yes, pos
% roi = {'RAS'}; % yes, pos *
% roi = {'RAS2'}; % yes

% roi = {{'LT','FS'}}; % 
% roi = {'LT'}; % 
% roi = {'LPS2'}; % 
% roi = {'LPS'}; % 
roi = {'PS'}; % 
% roi = {'RPS'}; % 
% roi = {'RPS2'}; % 
% roi = {'RT'}; % 

% roi = {'LPI2'}; % yes
% roi = {'PI'}; % no
% roi = {'RPI2'}; % yes

% roi = {{'E23','Fz','E3'}}; % AF3 Fz AF4
% roi = {{'E67','Pz','E77'}}; % PO3 Pz P04

% % % alpha
% freqs = ana.freq.alpha;
% % freqs = ana.freq.alpha_lower;
% % freqs = ana.freq.alpha_upper;
% % latencies = [0.6 1.0]; % img
% % roi = {'LAS2'};
% 
% % latencies = [0.1 0.3]; % word, LT, early alpha effect Spac x Mem
% % latencies = [0.4 0.7]; % 
% % latencies = [0.5 0.9]; % 
% % latencies = [0.4 1.0]; % 
% roi = {'LT'}; % word **
% 
% % latencies = [0.3 0.7]; % word
% % latencies = [0.5 0.7]; % word
% % latencies = [0.8 1.0]; % word
% % roi = {'LAS2'}; %
% % roi = {'PS'}; % word **
% % roi = {'LPS'}; % 
% % roi = {'PI'}; % word ** 
% 
% % latencies = [0.4 0.6];
% % latencies = [0.4 0.7];
% % roi = {'FC'}; % no
% % roi = {'FS'}; % maybe
% % roi = {'FS2'}; % yes
% % roi = {'PS'}; %
% % roi = {'LPI2'}; %

latency = cell(1,size(latencies,1));
for i = 1:length(latency)
  latency{i} = sprintf('%dto%d',round(latencies(i,1)*1000),round(latencies(i,2)*1000));
end
latStr = sprintf(repmat('_%s',1,length(latency)),latency{:});
latStr = latStr(2:end);

freq_str = sprintf('%dfreq%dto%d',size(freqs,1),round(freqs(1,1)),round(freqs(end,end)));
freqIdx = (nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,freqs(1,1)):nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,freqs(1,2)));

factorNames = {'spacings', 'oldnew', 'memConds', 'roi', 'latency'};

nVariables = nan(size(factorNames));
keepTheseFactors = false(size(factorNames));
levelNames_teg = cell(size(factorNames)); % TEG
for c = 1:length(factorNames)
  % need to have a variable set to this exact factor name
  thisFac = eval(factorNames{c});
  nVariables(c) = length(thisFac);
  if nVariables(c) > 1
    keepTheseFactors(c) = true;
  end
  
  % go through the levels within this factor and save strings
  thisFac_str = cell(1,nVariables(c));
  for l = 1:length(thisFac)
    if iscell(thisFac{l})
      tf = thisFac{l};
      thisFac_str{l} = sprintf(repmat('_%s',1,length(tf)),tf{:});
      thisFac_str{l} = thisFac_str{l}(2:end);
    elseif ischar(thisFac{l})
      thisFac_str{l} = thisFac{l};
    end
  end
  levelNames_teg{c} = thisFac_str; % TEG
end

variableNames = cell(1,prod(nVariables));
levelNames = cell(prod(nVariables),length(factorNames));

ses=1;
nSub = sum(~exper.badSub);
anovaData = nan(nSub,prod(nVariables));
rmaov_data_teg = nan(nSub*prod(nVariables),length(factorNames) + 2);

fprintf('Collecting ANOVA data for %d subjects:\n\t',nSub);
fprintf('%s\n',stimType);
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(spacings)),spacings{:}),factorNames{1});
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(oldnew)),oldnew{:}),factorNames{2});
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(memConds)),memConds{:}),factorNames{3});
if iscell(roi{1})
  fprintf('%d ROIs (%s),',length(roi),factorNames{4});
elseif ischar(roi{1})
  fprintf('%s (%s),',sprintf(repmat(' %s',1,length(roi)),roi{:}),factorNames{4});
end
fprintf('%s (%s),',latStr,factorNames{5});
fprintf('\n\tFreq: %s...',freq_str);

lnDone = false;
vnDone = false;
subCount = 0;
rmCount = 0;
for sub = 1:length(exper.subjects)
  if ~exper.badSub(sub)
    subCount = subCount + 1;
  else
    continue
  end
  for ses = 1:length(exper.sesStr)
    lnCount = 0;
    vnCount = 0;
    
    for sp = 1:length(spacings)
      for on = 1:length(oldnew)
        for mc = 1:length(memConds)
          cond_str = [];
          cond_str_tmp = [];
          if strcmp(spacings{sp},'onePres')
            % single presentation or first presentation
            if strcmp(memConds{mc},'all')
              cond_str = sprintf('%s%s',stimType,spacings{sp});
            end
          elseif strcmp(spacings{sp},'mass') || strcmp(spacings{sp},'spac')
            cond_str = sprintf('%s%s%s_%s_%s',stimType,memType,memConds{mc},spacings{sp},oldnew{on});
            if strcmp(memConds{mc},'all')
              % manually input rc and fo
              cond_str_tmp = {sprintf('%s%s%s_%s_%s',stimType,memType,'rc',spacings{sp},oldnew{on}), sprintf('%s%s%s_%s_%s',stimType,memType,'fo',spacings{sp},oldnew{on})};
            end
          end
          
          for r = 1:length(roi)
            if iscell(roi{r})
              roi_str = sprintf(repmat('%s',1,length(roi{r})),roi{r}{:});
            elseif ischar(roi{r})
              roi_str = roi{r};
            end
            if ~isempty(cond_str_tmp)
              chanIdx = ismember(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,roi{r})})));
            else
              chanIdx = ismember(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,roi{r})})));
            end
            
            for lat = 1:length(latency)
              if ~isempty(cond_str_tmp)
                tbeg = nearest(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.time,latencies(lat,1));
                tend = nearest(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.time,latencies(lat,2));
              else
                tbeg = nearest(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,1));
                tend = nearest(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,2));
              end
              latIdx = tbeg:tend;
              
              if ~lnDone
                lnCount = lnCount + 1;
                levelNames{lnCount,1} = spacings{sp};
                levelNames{lnCount,2} = oldnew{on};
                levelNames{lnCount,3} = memConds{mc};
                levelNames{lnCount,4} = roi{r};
                levelNames{lnCount,5} = latency{lat};
              end
              
              vnCount = vnCount + 1;
              if ~vnDone
                variableNames{vnCount} = sprintf('Y%d',vnCount);
              end
              
              if ~isempty(cond_str_tmp)
                thisData = (mean(mean(mean(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.(measure)(chanIdx,freqIdx,latIdx),3),2),1) + ...
                  mean(mean(mean(data_pow.(exper.sesStr{ses}).(cond_str_tmp{2}).sub(sub).data.(measure)(chanIdx,freqIdx,latIdx),3),2),1)) / 2;
              else
                thisData = mean(mean(mean(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.(measure)(chanIdx,freqIdx,latIdx),3),2),1);
              end
              anovaData(subCount,vnCount) = thisData;

              rmCount = rmCount + 1;
              rmaov_data_teg(rmCount,:) = [anovaData(subCount,vnCount) sp on mc r lat sub];
            end
          end
        end
      end
    end
    lnDone = true;
    vnDone = true;
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

% TEG RM ANOVA

fprintf('=======================================\n');
fprintf('This ANOVA: %s, Freq: %s\n\n',stimType,freq_str);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: %s, Freq: %s\n',stimType,freq_str);
fprintf('=======================================\n');

%% plot spacing x memory x time interaction

theseSub = noNans & passTrlThresh;

ses=1;

RgH_rc_spac = squeeze(mean_similarity.(sprintf('%s_RgH_rc_spac',data_str))(theseSub,ses,latInd));
RgH_fo_spac = squeeze(mean_similarity.(sprintf('%s_RgH_fo_spac',data_str))(theseSub,ses,latInd));
RgH_rc_mass = squeeze(mean_similarity.(sprintf('%s_RgH_rc_mass',data_str))(theseSub,ses,latInd));
RgH_fo_mass = squeeze(mean_similarity.(sprintf('%s_RgH_fo_mass',data_str))(theseSub,ses,latInd));

plotMeanLines = true;
plotSub = false;
plotSubLines = false;

if plotMeanLines
  s_rc_mark = 'ko-';
  s_fo_mark = 'ko-.';
  m_rc_mark = 'rx-';
  m_fo_mark = 'rx-.';
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
    s_fo_mark_sub = 'ko-.';
    m_rc_mark_sub = 'rx-';
    m_fo_mark_sub = 'rx-.';
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
      plot(RgH_rc_spac(s,:),s_rc_mark_sub,'LineWidth',1);
      plot(RgH_fo_spac(s,:),s_fo_mark_sub,'LineWidth',1);
      % massed
      plot(RgH_rc_mass(s,:),m_rc_mark_sub,'LineWidth',1);
      plot(RgH_fo_mass(s,:),m_fo_mark_sub,'LineWidth',1);
    end
  else
    for t = 1:length(latInd)
      % spaced
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), RgH_rc_spac(:,t),s_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), RgH_fo_spac(:,t),s_fo_mark_sub,'LineWidth',1);
      % massed
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), RgH_rc_mass(:,t),m_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), RgH_fo_mass(:,t),m_fo_mark_sub,'LineWidth',1);
    end
  end
end
if plotMeanLines
  % spaced
  hsr = plot(mean(RgH_rc_spac,1),s_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  hsf = plot(mean(RgH_fo_spac,1),s_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  % massed
  hmr = plot(mean(RgH_rc_mass,1),m_rc_mark,'LineWidth',3,'MarkerSize',meanSizeM);
  hmf = plot(mean(RgH_fo_mass,1),m_fo_mark,'LineWidth',3,'MarkerSize',meanSizeM);
else
  for t = 1:length(latInd)
    % spaced
    hsr = plot(t,mean(RgH_rc_spac(:,t),1),s_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    hsf = plot(t,mean(RgH_fo_spac(:,t),1),s_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    % massed
    hmr = plot(t,mean(RgH_rc_mass(:,t),1),m_rc_mark,'LineWidth',3,'MarkerSize',meanSizeM);
    hmf = plot(t,mean(RgH_fo_mass(:,t),1),m_fo_mark,'LineWidth',3,'MarkerSize',meanSizeM);
  end
end

% horiz
plot([-length(latInd)-1, length(latInd)+1], [0 0],'k--','LineWidth',2);

hold off
axis square
xlim([0.75 length(latInd)+0.25]);
ylim([-0.35 0.35]);

set(gca,'XTick', 1:length(latInd));
set(gca,'XTickLabel',latencySec);
xlabel('Time (Sec)');

ylabel('Neural Similarity');

title(sprintf('Spacing \\times Memory \\times Time: %s',data_str));
legend([hsr, hsf, hmr, hmf],{'Spaced Rc','Spaced Fo','Massed Rc','Massed Fo'},'Location','North');

% ticFontSize = 20;
ticFontSize = 18;
publishfig(gcf,0,ticFontSize,[],[]);

% print(gcf,'-depsc2',sprintf('~/Desktop/similarity_spacXmemXtime_%s_%s_%s_%s_%s_%s.eps',data_str,origDataType,roi_str,latStr,eig_criterion,sim_method));
