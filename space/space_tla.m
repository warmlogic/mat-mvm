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
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
else
  error('Data directory not found.');
end

% procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
% procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');
procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');
% procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_continuousICA_nsClassic/tla');
% procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_continuousICA_ftAuto/tla');

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

files.figPrintFormat = 'png';
files.saveFigs = true;

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
  'img_RgM_spac_p1','img_RgM_spac_p2','img_RgM_mass_p1','img_RgM_mass_p2' ...
  } ...
  { ...
  'word_onePres' ...
  'word_RgH_rc_spac_p1','word_RgH_rc_spac_p2','word_RgH_rc_mass_p1','word_RgH_rc_mass_p2' ...
  'word_RgH_fo_spac_p1','word_RgH_fo_spac_p2','word_RgH_fo_mass_p1','word_RgH_fo_mass_p2' ...
  'word_RgM_spac_p1','word_RgM_spac_p2','word_RgM_mass_p1','word_RgM_mass_p2' ...
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
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
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
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
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
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_image'))) ...
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
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 1 & lag > 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0 & spaced == 0 & lag == 0 & presNum == 1',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
    sprintf('eventNumber == %d & targ == 1 & cr_recog_acc == 0  & spaced == 0 & lag == 0 & presNum == 2',find(ismember(exper.eventValues{sesNum},'multistudy_word'))) ...
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

%% load in the subject data

% % make sure ana.eventValues is set properly
% if ~iscell(ana.eventValues{1})
%   ana.eventValues = {ana.eventValues};
% end
% if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
%   ana.eventValues = {exper.eventValues};
% end

% keeptrials = true;
keeptrials = false;
% [data_tla,exper] = mm_ft_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');
[data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo',true);

% %% get rid of the bad channels
% 
% cfg = [];
% cfg.printRoi = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% [data_tla] = mm_rmBadChan(cfg,exper,ana,data_tla);

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

%% save

% saveDir = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
% saveDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
saveDir = dirs.saveDirProc;
% save(fullfile(saveDir,'space_word_img_data_ga_tla.mat'),'data_tla','ga_tla','exper','ana','dirs','files','-v7.3');
save(fullfile(saveDir,'space_word_img_data_tla.mat'),'data_tla','exper','ana','dirs','files','-v7.3');
% clear data_tla

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

% loadDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
% loadDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');
loadDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');
% loadDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_continuousICA_nsClassic/tla');
% loadDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_continuousICA_ftManual/tla');

% loadDir = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
% loadDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';

% load(fullfile(loadDir,'space_word_img_data_ga_tla.mat'));
load(fullfile(loadDir,'space_word_img_data_tla.mat'));

replaceDataType = {};

[dirs] = mm_checkDirs(dirs,replaceDataType,subDir);

%% decide who to kick out based on trial counts

% Subjects with bad behavior
% exper.badBehSub = {{}};
% exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE030','SPACE039'}};
exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE030'}};

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
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.ylim = [-10 10];
cfg_ft.layout = ft_prepare_layout([],ana);
sub = 1;
ses = 1;
for typ = 1:length(ana.eventValues{ses})
  for ev = 1:length(ana.eventValues{ses}{typ})
    figure
    ft_multiplotER(cfg_ft,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{ev}).sub(sub).data);
    title(strrep(ana.eventValues{ses}{typ}{ev},'_','-'));
  end
end

% cfg_ft = [];
% % cfg_ft.channel = {'E70'};
% cfg_ft.channel = {'E83'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% sub = 1;
% ses = 1;
% typ = 1;
% ev = 1;
% ev2 = 2;
% figure
% ft_singleplotER(cfg_ft,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{ev}).sub(sub).data,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{ev2}).sub(sub).data);

% cfg_ft = [];
% cfg_ft.channel = {'E20'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% %figure
% %ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data,data_tla.(exper.eventValues{2}).sub(1).ses(1).data,data_tla.(exper.eventValues{3}).sub(1).ses(1).data);
% %legend(strrep(exper.eventValues,'_','-'),'Location','SouthEast');
% figure
% cfg_ft.graphcolor = 'b';
% ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% legend(strrep(exper.eventValues{1},'_','-'),'Location','SouthEast');
% hold on
% plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
% plot([0 0],[-5 5],'k--'); % vertical
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cfg_ft = [];
% cfg_ft.channel = {'all'};
% cfg_ft.latency = [.5 .8];
% data_topo = ft_timelockanalysis(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% cfg_ft = [];
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.interactive = 'yes';
% cfg_ft.colormap = 'jet';
% %cfg_ft.colormap = 'hot';
% cfg_ft.colorbar = 'yes';
% cfg_ft.xlim = 'maxmin';
% cfg_ft.layout = ft_prepare_layout([],ana);
% figure
% ft_topoplotER(cfg_ft,data_topo);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Remove the dof field. I'm not sure what degrees of freedom refers to, but the tutorial says to.
% for evVal = 1:length(exper.eventValues)
%   for sub = 1:length(exper.subjects)
%     for ses = 1:length(exper.sesStr)
%       if isfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof')
%         data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data = rmfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof');
%       end
%     end
%   end
% end

%% lowpass filter and segment for ERPs

data_tla_backup = data_tla;

lpfilt = true;

if lpfilt
  sampleRate = exper.sampleRate;
  lpfreq = 40;
  lofiltord = 3;
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
          % process individual trials
          if isfield(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,'trial')
            for i = 1:size(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial,1)
              data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial(i,:,:) = ft_preproc_lowpassfilter( ...
                squeeze(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial(i,:,:)), ...
                sampleRate,lpfreq,lofiltord,lpfilttype);
            end
          end
          
          % process subject average
          if isfield(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,'avg')
            data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.avg = ft_preproc_lowpassfilter( ...
              data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.avg, ...
              sampleRate,lpfreq,lofiltord,lpfilttype);
          end
        end
        
        % select
        data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_selectdata(cfg_sel,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
        
      end
      
%       % grand average
%       if exist('ga_tla','var')
%         if lpfilt
%           if isfield(ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}),'avg')
%             ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).avg = ft_preproc_lowpassfilter( ...
%               ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).avg, ...
%               sampleRate,lpfreq,lofiltord,lpfilttype);
%           end
%         end
%         
%         % select
%         ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = ft_selectdata_new(cfg_sel,ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}));
%       end
    end
  end
end

%% get the grand average

% might need to remove avg field from data until this bug is resolved:
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=2471

ga_tla = struct;

for ses = 1:length(exper.sesStr)
  % set up strings to put in grand average function
  cfg_ana = [];
  cfg_ana.is_ga = 0;
  %cfg_ana.conditions = ana.eventValues{ses};
  cfg_ana.data_str = 'data_tla';
  % cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);
  %cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
  
  cfg_ft = [];
  cfg_ft.keepindividual = 'no';
  %cfg_ft.parameter = 'trial';
  for typ = 1:length(ana.eventValues{ses})
    cfg_ana.conditions = ana.eventValues{ses}{typ};
    cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
    
    for evVal = 1:length(ana.eventValues{ses}{typ})
      %tic
      fprintf('Running ft_timelockgrandaverage on %s...',ana.eventValues{ses}{typ}{evVal});
      %ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{evVal}){ses}));
      ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
      fprintf('Done.\n');
      %toc
    end
  end
end

% % turn keeptrial data into average for statistical functions because proper
% % processing of dimord is currently broken
% 
% data_tla_avg = struct;
% 
% cfg = [];
% cfg.keeptrials = 'no';
% 
% for ses = 1:length(exper.sesStr)
%   for typ = 1:length(ana.eventValues{ses})
%     for evVal = 1:length(ana.eventValues{ses}{typ})
%       for sub = 1:length(exper.subjects)
%         fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
%         if isfield(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,'avg')
%           data_tla_avg.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_timelockanalysis(cfg,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
%         end
%       end
%     end
%   end
% end

%% plot the conditions - simple

cfg_ft = [];
cfg_ft.xlim = [-0.2 1.0];
% cfg_ft.xlim = [-0.2 0.4];
% cfg_ft.xlim = [-1.0 2.0];
cfg_ft.parameter = 'avg';

cfg_ft.layout = ft_prepare_layout([],ana);

cfg_plot = [];

cfg_plot.is_ga = true;
% cfg_plot.is_ga = false;
cfg_plot.excludeBadSub = 1;

% %cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

% cfg_plot.rois = {{'LAS'},{'RAS'}};
% cfg_plot.ylims = [-4 2; -4 2];
% cfg_plot.legendlocs = {'SouthEast','NorthWest'};

% cfg_plot.rois = {{'FC'},{'FS'}};
% cfg_plot.ylims = [-4 2; -4 2];
% cfg_plot.legendlocs = {'SouthEast','NorthWest'};

% % cfg_plot.rois = {{'LPS'},{'RPS'}};
% cfg_plot.rois = {{'LPI3'},{'RPI3'}};
% cfg_plot.ylims = [-1 4; -1 4];
% cfg_plot.legendlocs = {'NorthWest','NorthWest'};

% cfg_plot.rois = {{'Fz'},{'Cz'},{'Pz'},{'Oz'}};
% cfg_plot.ylims = [-3 2; -2 3; -1 4; -1 4];
% cfg_plot.legendlocs = {'NorthEast','NorthEast','SouthEast','SouthEast'};

% LPC
% cfg_plot.rois = {{'E62','E72','E76','E77','E78','E84','E85'}}; % Centered on E77 *** using this for spacing analysis
% cfg_plot.rois = {{'C'}};
% cfg_plot.rois = {{'Pz'}}; % Pz
% cfg_plot.rois = {{'PS'}}; % Centered on Pz
% cfg_plot.rois = {{'PS2'}}; % Centered on E72
% cfg_plot.rois = {{'LPS'},{'RPS'}};
cfg_plot.rois = {{'LPS2'},{'RPS2'}};
% cfg_plot.rois = {{'LPS2','RPS2'}};
cfg_plot.ylims = [-1 5; -1 5];
% cfg_plot.legendlocs = {'SouthEast'};
cfg_plot.legendlocs = {'NorthWest','NorthWest'};

% % N400
% % cfg_plot.rois = {{'LAS'}}; % Center=20
% % cfg_plot.rois = {{'FS2'}}; % Center=6
% cfg_plot.rois = {{'C'}};
% % cfg_plot.rois = {{'Cz'}};
% cfg_plot.ylims = [-4 2];
% % cfg_plot.legendlocs = {'NorthEast'};
% cfg_plot.legendlocs = {'SouthWest'};

% % N2
% % cfg_plot.rois = {{'E58'},{'E96'}}; % T5, T6
% % cfg_plot.rois = {{'E50','E51','E57','E58','E59','E64','E65'}}; % T5 (L)
% % cfg_plot.rois = {{'E90','E91','E95','E96','E97','E100','E101'}}; % T6 (R)
% cfg_plot.rois = {{'E50','E51','E57','E58','E59','E64','E65','E90','E91','E95','E96','E97','E100','E101'}}; % T5+T6
% cfg_plot.ylims = [-3 2; -3 2];
% cfg_plot.legendlocs = {'NorthEast','NorthEast'};

% cfg_plot.rois = {{'Oz'},{'PI'}};
% cfg_plot.ylims = [-2 3; -2 3];
% cfg_plot.legendlocs = {'NorthEast','NorthEast'};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

%cfg_plot.condByTypeByROI = repmat({{{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}},size(cfg_plot.rois));

sesNum = 1;
% cfg_plot.condByROI = repmat(ana.eventValues(sesNum),size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'Face' 'House'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'word_RgH_rc_spac_p1', 'word_RgH_rc_spac_p2', 'word_RgH_fo_spac_p1', 'word_RgH_fo_spac_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_rc_spac_p1', 'img_RgH_rc_spac_p2', 'img_RgH_fo_spac_p1', 'img_RgH_fo_spac_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'word_RgH_rc_mass_p1', 'word_RgH_rc_mass_p2', 'word_RgH_fo_mass_p1', 'word_RgH_fo_mass_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_rc_mass_p1', 'img_RgH_rc_mass_p2', 'img_RgH_fo_mass_p1', 'img_RgH_fo_mass_p2'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'word_spac_p2', 'word_mass_p2', 'word_onePres'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_spac_p2', 'img_mass_p2', 'img_onePres'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'word_RgH_spac_p2', 'word_RgH_mass_p2', 'word_RgH_spac_p1', 'word_RgH_mass_p1', 'word_onePres'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_spac_p2', 'img_RgH_mass_p2', 'img_RgH_spac_p1', 'img_RgH_mass_p1', 'img_onePres'}},size(cfg_plot.rois));

% VanSEtal2007
cfg_plot.condByROI = repmat({{'word_RgH_rc_spac_p1', 'word_RgH_fo_spac_p1', 'word_RgH_rc_mass_p1', 'word_RgH_fo_mass_p1', 'word_RgH_rc_spac_p2', 'word_RgH_fo_spac_p2', 'word_RgH_rc_mass_p2', 'word_RgH_fo_mass_p2' ,'word_onePres'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'word_RgH_rc_spac_p2', 'word_RgH_fo_spac_p2', 'word_RgH_rc_mass_p2', 'word_RgH_fo_mass_p2' ,'word_onePres'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_rc_spac_p1', 'img_RgH_fo_spac_p1', 'img_RgH_rc_mass_p1', 'img_RgH_fo_mass_p1', 'img_RgH_rc_spac_p2', 'img_RgH_fo_spac_p2', 'img_RgH_rc_mass_p2', 'img_RgH_fo_mass_p2' ,'img_onePres'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_rc_spac_p2', 'img_RgH_fo_spac_p2', 'img_RgH_rc_mass_p2', 'img_RgH_fo_mass_p2' ,'img_onePres'}},size(cfg_plot.rois));

cfg_plot.axisxy = false;

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  if exist('linspecer','file')
    cfg_ft.graphcolor = linspecer(length(cfg_plot.conditions));
  end
  
  mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
%   mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,sesNum,data_tla);
  %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
end

%% find peaks
cfg = [];

% images
% cfg.conditions = cellflat(ana.eventValues{1}{1});
% words
% cfg.conditions = cellflat(ana.eventValues{1}{2});

% % words, all together
cfg.conditions = {'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};
% cfg.conditions = {'word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};

% P2 only
% cfg.conditions = {'word_onePres','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};
% cfg.conditions = {'word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};

% spaced
% cfg.conditions = {'word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2'};

% % massed
% cfg.conditions = {'word_RgH_rc_mass_p1','word_RgH_fo_mass_p1','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};

% % single presentation or first presentation
% cfg.conditions = {'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1'};
% cfg.conditions = {'word_onePres'};

% % NB: make sure pattern of conditions is seen at peak electrode(s)

% ==================================================
% LPC (positive late pareital component)
% ==================================================

% % step 1: find peak electrode in large area, collapsing across conditions
% cfg.datadim = 'elec';
% % cfg.roi = {'center101'};
% cfg.roi = {'posterior_noPeriph'};
% cfg.latency = [0.4 0.8]; % LPC
% cfg.order = 'descend'; % descend = positive peaks first
% % LPC: electrode cluster around E77 includes peak E84 E85 etc

% % step 2: find peak time at peak electrode(s)
% cfg.datadim = 'time';
% cfg.roi = {'E62','E72','E76','E77','E78','E84','E85'}; % Centered on E77 (596ms)
% % cfg.roi = {'RPS2'}; % Centered on E85 (500ms)
% % cfg.roi = {'RPS'}; % Centered on E86 (ms)
% % cfg.roi = {'LPS2','RPS2'}; % Bilateral, centered on E60+E85 (576ms)
% cfg.latency = [0.4 0.8]; % LPC
% cfg.order = 'descend'; % descend = positive peaks first

% % step 3: select a window for analysis around peak; cfg.outputSubjects=true
% % lpcPeak = 0.588;
% lpcPeak = 0.576; % bilateral
% cfg.datadim = 'time';
% cfg.roi = {'E62','E72','E76','E77','E78','E84','E85'}; % Centered on E77
% % cfg.roi = {'RPS2'}; % Centered on E85
% % cfg.roi = {'LPS2','RPS2'}; % Bilateral, centered on E60+E85
% % cfg.latency = [lpcPeak-0.05 lpcPeak+0.05]; % LPC - around GA peak (space+mass) +/- 50
% % cfg.latency = [lpcPeak-0.1 lpcPeak+0.1]; % LPC - around GA peak (space+mass) +/- 100
% cfg.latency = [lpcPeak-0.15 lpcPeak+0.15]; % LPC - around GA peak (space+mass) +/- 150
% cfg.avgovertime = true;
% cfg.order = 'descend'; % descend = positive peaks first

% ==================================================
% N400 (negative frontocentral component)
% ==================================================

% % step 1: find peak electrode in large area, collapsing across conditions
% cfg.datadim = 'elec';
% cfg.roi = {'center101'};
% cfg.latency = [0.3 0.5]; % N400
% cfg.order = 'ascend'; % ascend = negative peaks first

% % step 2: find peak time at peak electrode(s)
% cfg.datadim = 'time';
% % cfg.roi = {'FS2'}; % Centered on E6 (360ms)
% cfg.roi = {'C'}; % Centered on Cz (372ms)
% cfg.latency = [0.3 0.5]; % N400
% cfg.order = 'ascend'; % ascend = negative peaks first

% % step 3: select a window for analysis around peak
% n400Peak = 0.360; % FS2
% % n400Peak = 0.372; % C
% cfg.datadim = 'time';
% cfg.roi = {'FS2'}; % Centered on E6
% % cfg.roi = {'C'}; % Centered on Cz
% cfg.latency = [n400Peak-0.05 n400Peak+0.05]; % N400 - around GA peak (space+mass) +/- 50
% % cfg.latency = [n400Peak-0.1 n400Peak+0.1]; % N400 - around GA peak (space+mass) +/- 100
% cfg.avgovertime = true;
% cfg.order = 'ascend'; % ascend = negative peaks first

% ==================================================
% N2 (negative posterior attentional component)
% ==================================================

% % step 1: find peak electrode in large area, collapsing across conditions
% cfg.datadim = 'elec';
% % cfg.roi = {'center101'};
% % cfg.roi = {'posterior'};
% cfg.roi = {'posterior_noPeriph'};
% cfg.latency = [0.15 0.25]; % N2
% cfg.order = 'ascend'; % ascend = negative peaks first

% % step 2: find peak time at peak electrode(s)
% cfg.datadim = 'time';
% cfg.roi = {'E50','E51','E57','E58','E59','E64','E65'}; % Centered on E58/T5 (172ms)
% % cfg.roi = {'E50','E51','E57','E58','E59','E64','E65','E90','E91','E95','E96','E97','E100','E101'}; % Centered on T5+T6 (172ms)
% cfg.latency = [0.1 0.3]; % N2
% cfg.order = 'ascend'; % ascend = negative peaks first

% % step 3: select a window for analysis around peak
% n2Peak = 0.172; % 'E58'
% cfg.datadim = 'time';
% cfg.roi = {'E50', 'E51', 'E57', 'E58', 'E59', 'E64', 'E65'}; % Centered on 
% cfg.latency = [n2Peak-0.05 n2Peak+0.05]; % around GA peak (space+mass) +/- 50
% cfg.avgovertime = true;
% cfg.order = 'ascend'; % ascend = negative peaks first

cfg.is_ga = false;
% cfg.is_ga = true;
cfg.outputSubjects = false;
% cfg.outputSubjects = false;
cfg.sesNum = 1;

cfg.plotit = true;
% cfg.voltlim = [-1 5]; % LPC
cfg.voltlim = [-5 5]; % N400
% cfg.voltlim = [-3 3]; % N2

% % only for datadim='elec' and datadim='peak2peak'
% cfg.plottype = 'topo';
% % cfg.plottype = 'multi';

% % only for datadim='peak2peak'
% cfg.datadim = 'peak2peak';
% cfg.order = 'descend'; % descend = positive peaks first
% cfg.roi = {'posterior'};
% cfg.pospeak = [0.08 0.14];
% cfg.negpeak = [0.14 0.2];

% peakInfo = mm_findPeak(cfg,ana,exper,ga_tla);
peakInfo = mm_findPeak(cfg,ana,exper,data_tla);

%% subplots of each subject's ERPs

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'LAS'},{'LPS'}};
% cfg_plot.rois = {{'LPS'},{'RPS'}};
% cfg_plot.rois = {{'Cz'},{'Pz'}};
% cfg_plot.rois = {{'E70'},{'E83'}};
% cfg_plot.rois = {{'E72'}};
% cfg_plot.rois = {{'LPS'}};
cfg_plot.rois = {{'LPS2'}};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-0.2 1.0];
cfg_plot.ylim = [-10 10];
cfg_plot.parameter = 'avg';

% 18 and 21 seem bery noisy

cfg_plot.titleTrialCount = false;

% cfg_plot.rois = {{'E70'},{'E83'}};
% cfg_plot.ylims = [-10 10; -10 10];

% cfg_plot.rois = {{'E83'}};
% cfg_plot.xlim = [-0.2 1.0];
% cfg_plot.ylim = [-10 10];

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% % abbreviations for the condition types
% cfg_plot.typesByROI = {...
%   {'C2','C6'},...
%   {'C2','C6'}};

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

ses = 1;

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  cfg_plot.types = {'image','word'};
  
  mm_ft_subjplotER(cfg_plot,ana,exper,data_tla,ses);
end

% %% plot the conditions
% 
% cfg_ft = [];
% cfg_ft.xlim = [-0.2 1.5];
% cfg_ft.parameter = 'avg';
% 
% cfg_plot = [];
% 
% %cfg_plot.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -1 6; -1 6; -1 6];
% % vertical solid lines to plot
% cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
% cfg_plot.plotLegend = 0;
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
% cfg_plot.plotTitle = 0;
% 
% cfg_plot.is_ga = 1;
% cfg_plot.excludeBadSub = 1;
% 
% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% 
% % cfg_plot.condByTypeByROI = {...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   };
% % cfg_plot.typesByROI = repmat({{'C2','C6'}},size(cfg_plot.condByTypeByROI));
% 
% %cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));
% 
% for r = 1:length(cfg_plot.rois)
%   cfg_plot.roi = cfg_plot.rois{r};
%   %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
%   cfg_plot.conditions = cfg_plot.condByROI{r};
%   cfg_ft.ylim = cfg_plot.ylims(r,:);
%   cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
%   if cfg_plot.plotLegend
%     cfg_plot.legendloc = cfg_plot.legendlocs{r};
%   end
%   
%   mm_ft_plotERP(cfg_ft,cfg_plot,ana,exper,files,dirs,ga_tla);
% end

% %% standardize data to have mean=0 and std=1 across all trials and elecs within a condition
% 
% parameter = 'trial';
% ses = 1;
% 
% % collapse across time and channel and across both conditions; trial is the
% % only remaining dimension
% %
% % NB: DO SEPARATELY FOR THE TRAINING AND TEST SETS
% 
% standardize = false;
% 
% if standardize
%   fprintf('Standardizing the data...');
%   
%   fn = fieldnames(data_tla);
%   for f = 1:length(fn)
%     for sub = 1:length(exper.subjects)
% 
%   
%   % timeInd = data1.time >= cfg_stat.latency(1) & data1.time <= cfg_stat.latency(2) + 0.0001;
%   % data_cat = cat(1,data1.(parameter)(:,:,timeInd),data2.(parameter)(:,:,timeInd));
%   %data_cat = cat(1,data1.(parameter),data2.(parameter));
%   dim = size(data_tla.(fn{f}).sub(sub).ses(ses).data.(parameter));
%   dat = reshape(data_tla.(fn{f}).sub(sub).ses(ses).data.(parameter), dim(1), prod(dim(2:end)));
%   
%   % m = dml.standardizer;
%   % m = m.train(dat);
%   % dat_z2 = reshape(m.test(dat),dim);
%   
%   mu = nanmean(dat);
%   sigma = nanstd(dat);
%   
%   dat_z = dat;
%   idx = ~isnan(mu);
%   dat_z(:,idx) = bsxfun(@minus,dat(:,idx),mu(idx));
%   idx = ~isnan(sigma) & ~(sigma==0);
%   dat_z(:,idx) = bsxfun(@rdivide,dat_z(:,idx),sigma(idx));
%   dat_z_ravel = reshape(dat_z,dim);
%   
%   % put it back
%   data_tla.(fn{f}).sub(sub).ses(ses).data.(parameter) = dat_z_ravel;
%   %data1.(parameter) = dat_z_ravel(1:size(data1.(parameter),1),:,:);
%   %data2.(parameter) = dat_z_ravel(size(data2.(parameter),1)+1:size(data1.(parameter),1)+size(data2.(parameter),1),:,:);
%   
%     end
%   end
%   
%   fprintf('Done.\n');
% else
%   fprintf('Not standardizing data!\n');
% end

%% make some GA plots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
cfg_ft.xlim = [-0.2 1.0]; % time
% cfg_ft.xlim = [-0.2 0.5]; % time
%cfg_ft.xlim = [-0.2 2.0]; % time

cfg_ft.parameter = 'avg';

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

cfg_plot.ftFxn = 'ft_singleplotER';

% % N1 @ T5/T6
% % cfg_plot.rois = {{'E50','E51','E57','E58','E59','E64','E65'},{'E90','E91','E95','E96','E97','E100','E101'}};
% cfg_plot.rois = {{'E50','E51','E57','E58','E59','E64','E65'}};
% cfg_plot.ylims = [-3 4; -3 4];
% cfg_plot.legendlocs = {'SouthEast','SouthEast'};
% n1Peak = 0.172; % E58
% cfg_plot.x_bounds = [n1Peak-0.05 n1Peak+0.05]; % time

% % N400 @ C
% % cfg_plot.rois = {{'FS'}};
% % cfg_plot.rois = {{'LAS'}};
% cfg_plot.rois = {{'C'}};
% cfg_plot.ylims = [-3 4];
% cfg_plot.legendlocs = {'NorthEast'};
% n400Peak = 0.372; % C
% cfg_plot.x_bounds = [n400Peak-0.05 n400Peak+0.05]; % time

% % LPC @ E77
cfg_plot.rois = {{'E62','E72','E76','E77','E78','E84','E85'}};
cfg_plot.ylims = [-1 6];
cfg_plot.legendlocs = {'NorthWest'};
lpcPeak = 0.596; % centered on E77
cfg_plot.x_bounds = [lpcPeak-0.1 lpcPeak+0.1]; % time

cfg_plot.condByROI = repmat({{{'word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2'}}},size(cfg_plot.rois));
cfg_plot.rename_condByROI = repmat({{{'Space P1 Recalled','Space P1 Forgot','Space P2 Recalled','Space P2 Forgot'}}},size(cfg_plot.rois));
cfg_ft.graphcolor = 'kgrm';
cfg_ft.linestyle = {'-','--','-','--'};
% cfg_plot.condByROI = repmat({{{'word_RgH_rc_mass_p1','word_RgH_fo_mass_p1','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'Mass P1 Recalled','Mass P1 Forgot','Mass P2 Recalled','Mass P2 Forgot'}}},size(cfg_plot.rois));
% cfg_ft.graphcolor = 'kgbc';
% cfg_ft.linestyle = {'-','--','-','--'};


% cfg_plot.rois = {{'LPS2'},{'RPS2'}};
% cfg_plot.ylims = [-1 6; -1 6];
% cfg_plot.legendlocs = {'NorthWest','NorthWest'};

% % P4/E92
% cfg_plot.rois = {{'E85','E85','E91','E92','E93','E97','E98'}};
% cfg_plot.ylims = [-1 6];
% cfg_plot.legendlocs = {'NorthWest'};

% cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5; -2 5];
% cfg_plot.rois = {{'PI'}};

% cfg_plot.rois = {{'LPI2'}};
% cfg_plot.ylims = [-3 3];

% cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
% cfg_plot.x_bounds = [0.13 0.19];
% cfg_plot.x_bounds = [0.14 0.2];
cfg_plot.plotLegend = 1;

cfg_plot.xlabel = 'Time (s)';
cfg_plot.ylabel = 'Voltage (\muV)';
% cfg_plot.xlabel = '';
% cfg_plot.ylabel = '';

% cfg_plot.ftFxn = 'ft_topoplotER';
% % cfg_plot.ylims = [-2 2];
% % %cfg_plot.ylims = 'maxmin';
% cfg_ft.marker = 'off';
% % cfg_ft.marker = 'labels';
% % cfg_ft.markerfontsize = 9;
% 
% % cfg_ft.xlim = [0.13 0.19];
% cfg_ft.xlim = [0:0.1:1.0];

% %cfg_ft.comment = 'no';
% % cfg_plot.rois = {'all'};
% % cfg_ft.xlim = [0 1.2]; % time
% % cfg_plot.subplot = 1;
% cfg_plot.rois = {{'LAS'}};
% cfg_ft.xlim = [0.3 0.5]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [0.5 0.8]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [1.0 1.5]; % time

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% %cfg_plot.rois = {{'FS'},{'LAS','RAS'},{'LPS','RPS'}};
% %cfg_plot.rois = {{'FS'},{'PS'}};
% %cfg_plot.rois = {'E71'};
% cfg_plot.rois = {'all'};
% cfg_plot.ylims = [-5 5]; % voltage in multiplot
% %cfg_plot.ylims = repmat('maxmin',size(cfg_plot.rois,2),1); % voltage in multiplot

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
% cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{'FSC','FSI','N'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'FSC','FSI','CR'}}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));
%cfg_plot.rename_condByROI = repmat({{{'SC','SI','CR'}}},size(cfg_plot.rois));


% cfg_plot.conditions = {{'word_RgH_rc_spac_p2','word_onePres'},{'word_RgH_rc_mass_p2','word_onePres'},{'word_RgH_fo_spac_p2','word_onePres'},{'word_RgH_fo_mass_p2','word_onePres'}};
% cfg_plot.conditions = {{'word_RgH_rc_spac_p2','word_RgH_rc_mass_p2'},{'word_RgH_fo_spac_p2','word_RgH_fo_mass_p2'}};

% cfg_plot.condByROI = repmat({{{'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'OnePres','Space P1 Recalled','Space P1 Forgot','Mass P1 Recalled','Mass P1 Forgot'}}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{'word_onePres','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'OnePres','Space P2 Recalled','Space P2 Forgot','Mass P2 Recalled','Mass P2 Forgot'}}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'OnePres','Space P1 Recalled','Space P1 Forgot','Mass P1 Recalled','Mass P1 Forgot','Space P2 Recalled','Space P2 Forgot','Mass P2 Recalled','Mass P2 Forgot'}}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{{'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'OnePres','Space P1 Recalled','Space P1 Forgot','Space P2 Recalled','Space P2 Forgot'}}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{'word_onePres','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'OnePres','Mass P1 Recalled','Mass P1 Forgot','Mass P2 Recalled','Mass P2 Forgot'}}},size(cfg_plot.rois));


% cfg_plot.condByROI = repmat({{{'img_onePres','img_RgH_rc_spac_p1','img_RgH_fo_spac_p1','img_RgH_rc_spac_p2','img_RgH_fo_spac_p2'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'OnePres','Space P1 Recalled','Space P1 Forgot','Space P2 Recalled','Space P2 Forgot'}}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{{'img_onePres','img_RgH_rc_mass_p1','img_RgH_fo_mass_p1','img_RgH_rc_mass_p2','img_RgH_fo_mass_p2'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'OnePres','Mass P1 Recalled','Mass P1 Forgot','Mass P2 Recalled','Mass P2 Forgot'}}},size(cfg_plot.rois));


% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% cfg_plot.condByTypeByROI = repmat({{{'HSC2','HSI2','CR2'},{'HSC6','HSI6','CR6'}}},size(cfg_plot.rois));
% % cfg_plot.condByTypeByROI = {...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   };
% cfg_plot.typesByROI = repmat({{'C2','C6'}},size(cfg_plot.condByTypeByROI));

sesNum = 1;

% cfg_plot.ticFontSize = 20;
% cfg_plot.labelFontSize = 24;
cfg_plot.legendFontSize = 16;

files.saveFigs = true;
files.figPrintFormat = 'png';

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.types = cfg_plot.typesByROI{r};
  cfg_plot.rename_conditions = cfg_plot.rename_condByROI{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  
  if strcmp(cfg_plot.ftFxn,'ft_singleplotER')
    if isfield(cfg_plot,'x_bounds') && ~isempty(cfg_plot.x_bounds)
      cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
    else
      cfg_plot.x_bound = [];
    end
    if cfg_plot.plotLegend
      cfg_plot.legendloc = cfg_plot.legendlocs{r};
    end
    
%     if exist('linspecer','file')
%       cfg_ft.graphcolor = linspecer(length(cfg_plot.conditions{1}));
%     end
  end
  
  mm_ft_plotER(cfg_ft,cfg_plot,exper,ana,files,dirs,ga_tla,sesNum);
end

%% plot the contrasts

files.saveFigs = false;
files.figPrintFormat = 'png';

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_ft = [];
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
cfg_ft.colormap = 'jet';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'yes';
% cfg_ft.colorbar = 'no';

%cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';
%cfg_plot.conditions = {{'HSC2','CR2'},{'HSI2','CR2'},{'HSC2','HSI2'},{'HSC6','CR6'},{'HSI6','CR6'},{'HSC6','HSI6'}}; % {'H2','CR2'}, {'H6','CR6'},
%cfg_plot.conditions = {{'RH','RCR'},{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}};
% cfg_plot.conditions = {{'SC','CR'},{'SC','SI'},{'SI','CR'}};
%cfg_plot.conditions = {{'RHSC','RHSI'}};
%cfg_plot.conditions = {{'FSC','RSSI'}};
% cfg_plot.conditions = {{'Face','House'}};

% % cfg_plot.conditions = {{'word_RgH_rc_spac_p2','word_onePres'},{'word_RgH_rc_mass_p2','word_onePres'},{'word_RgH_fo_spac_p2','word_onePres'},{'word_RgH_fo_mass_p2','word_onePres'}};

% spacing - word
cfg_plot.conditions = {{'word_RgH_rc_mass_p2','word_RgH_rc_spac_p2'},{'word_RgH_fo_mass_p2','word_RgH_fo_spac_p2'}};
cfg_plot.cond_rename = {{'Mass P2 Recalled','Space P2 Recalled'},{'Mass P2 Forgot','Space P2 Forgot'}};

% spacing - img
% cfg_plot.conditions = {{'img_RgH_rc_mass_p2','img_RgH_rc_spac_p2'},{'img_RgH_fo_mass_p2','img_RgH_fo_spac_p2'}};
% cfg_plot.cond_rename = {{'Mass P2 Recalled','Space P2 Recalled'},{'Mass P2 Forgot','Space P2 Forgot'}};

% SME - word
cfg_plot.conditions = {{'word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'},{'word_RgH_rc_spac_p2','word_RgH_fo_spac_p2'}};
cfg_plot.cond_rename = {{'Mass P2 Recalled','Mass P2 Forgot'},{'Space P2 Recalled','Space P2 Forgot'}};

% % SME - img
% cfg_plot.conditions = {{'img_RgH_rc_mass_p2','img_RgH_fo_mass_p2'},{'img_RgH_rc_spac_p2','img_RgH_fo_spac_p2'}};
% cfg_plot.cond_rename = {{'Mass P2 Recalled','Mass P2 Forgot'},{'Space P2 Recalled','Space P2 Forgot'}};


cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-2 2]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
% cfg_ft.comment = 'no';
cfg_ft.comment = 'xlim';
cfg_ft.commentpos = 'middletop';

cfg_ft.shading = 'interp';

% cfg_plot.roi = {'PI'};
% cfg_plot.roi = {'E73'};
% cfg_plot.roi = {'E70','E71','E66','E83','E76','E84'};
% cfg_plot.roi = {'E70','E71','E66','E75','E83','E76','E84'};
% cfg_ft.xlim = [0.01 0.8]; % time


% cfg_plot.roi = {'E62','E72','E76','E77','E78','E84','E85'}; % E77, LPC
% lpcPeak = 0.596; % centered on E77
% cfg_ft.xlim = [lpcPeak-0.1 lpcPeak+0.1]; % time

% cfg_plot.roi = {'C'}; % N400
% n400Peak = 0.372; % C
% cfg_ft.xlim = [n400Peak-0.05 n400Peak+0.05]; % time

% cfg_plot.roi = {'E50','E51','E57','E58','E59','E64','E65'}; % E58, N1
% n1Peak = 0.172; % E58
% cfg_ft.xlim = [n1Peak-0.05 n1Peak+0.05]; % time


% cfg_plot.roi = {'LPS','RPS'};
cfg_plot.roi = {'LPS2','RPS2'};
% cfg_ft.xlim = [0.5 0.8]; % time

%cfg_plot.roi = {'RAS'};
%cfg_ft.xlim = [1.1 1.9]; % time
% cfg_ft.xlim = [0 1.5]; % time
% cfg_plot.roi = {'all'};

cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time
cfg_ft.xlim = (0:0.05:1.0); % time

% cfg_ft.xlim = [0.13 0.19];
% cfg_ft.xlim = [0.14 0.19];

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.xlim = [-0.2 1.5]; % time
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.ylim = [-1 1]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_ft.xlim = [-0.2 1.5]; % time
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.ylim = [-1 1]; % volt

sesNum = 1;
mm_ft_contrastER(cfg_ft,cfg_plot,exper,ana,files,dirs,ga_tla,sesNum);

%% line plots

files.saveFigs = true;
files.figPrintFormat = 'png';

cfg_ft = [];
cfg_ft.parameter = 'avg';

cfg_plot = [];

cfg_plot.conditions = {'word_RgH_rc_mass_p2','word_RgH_fo_mass_p2','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2'};
cfg_plot.plot_order = cfg_plot.conditions;
cfg_plot.rename_conditions = {'Mass P2 Recalled','Mass P2 Forgot','Space P2 Recalled','Space P2 Forgot'};


cfg_plot.linespec = {'bo','cx','ro','mx'};
% cfg_plot.markcolor = {'w','k','w','k','w','k','w','k'};
cfg_plot.markcolor = {'none','none','none','none'};
cfg_plot.marksize = 20;

% cfg_plot.roi = {'E50','E51','E57','E58','E59','E64','E65'}; % E58
% n1Peak = 0.172; % E58
% cfg_plot.latency = [n1Peak-0.05 n1Peak+0.05]; % time
% cfg_plot.ylim = [-2 2];
% cfg_plot.legendloc = 'NorthEast';

% cfg_plot.roi = {'C'};
% n400Peak = 0.372; % C
% cfg_plot.latency = [n400Peak-0.05 n400Peak+0.05]; % time
% cfg_plot.ylim = [-2 2];
% cfg_plot.legendloc = 'SouthWest';

cfg_plot.roi = {'E62','E72','E76','E77','E78','E84','E85'}; % E77
lpcPeak = 0.596; % centered on E77
cfg_plot.latency = [lpcPeak-0.1 lpcPeak+0.1]; % time
cfg_plot.ylim = [1.5 5.5];
cfg_plot.legendloc = 'SouthWest';

% cfg_plot.roi = {'LPS2'};
% 
% % cfg_plot.latency = [0.4 0.8];
% cfg_plot.latency = [0.5 0.8];
% % cfg_plot.latency = [0.4 0.6];
% % cfg_plot.latency = [0.6 0.8];
% 
% cfg_plot.ylim = [-1 4];

cfg_plot.plotLegend = true;
% cfg_plot.legendloc = 'SouthEast';
cfg_plot.xlabel = '';

sesNum = 1;
mm_ft_lineplotER_multiSes(cfg_ft,cfg_plot,ana,exper,files,dirs,sesNum,data_tla);


%% RM ANOVA

stimType = 'word_';
% stimType = 'img_';
memType = 'RgH_';

% % spacings = {'mass', 'spac', 'onePres'};
% spacings = {'mass', 'spac'};
% % oldnew = {'p1'};
% % oldnew = {'p2'};
% oldnew = {'p1', 'p2'};
% memConds = {'all'};

% didn't test new words, so can't assess memory, but can use p1
spacings = {'mass', 'spac'};
% spacings = {'spac'};
oldnew = {'p1', 'p2'};
% oldnew = {'p1'};
% oldnew = {'p2'};
memConds = {'rc','fo'};
% memConds = {'rc'};

measure = 'avg';

% % roi = {'LPS'};
% roi = {'LPS2'};
% latencies = [0.4 0.6; 0.604 0.804];

% roi = {'PS'};
% roi = {{'LPS2'},{'RPS2'}};
% roi = {{'LPS2','RPS2'}};
% latencies = [0.5 0.8];
% latencies = [0.5 0.648; 0.652 0.8];
% latencies = [0.604 0.804];
% latencies = [0.4 0.65; 0.55 0.8];

% roi = {{'E62','E72','E76','E77','E78','E84','E85'}}; % E77
% lpcPeak = 0.596; % centered on E77
% latencies = [lpcPeak-0.1 lpcPeak+0.1]; % time

roi = {{'C'}};
n400Peak = 0.372; % C
latencies = [n400Peak-0.05 n400Peak+0.05]; % time

% roi = {{'E50','E51','E57','E58','E59','E64','E65'}}; % E58
% n1Peak = 0.172; % E58
% latencies = [n1Peak-0.05 n1Peak+0.05]; % time

latency = cell(1,size(latencies,1));
for i = 1:length(latency)
  latency{i} = sprintf('%dto%d',latencies(i,1)*1000,latencies(i,2)*1000);
end
latStr = sprintf(repmat('_%s',1,length(latency)),latency{:});
latStr = latStr(2:end);

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
              %cond_str = sprintf('%s%s_%s',stimType,spacings{sp});
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
              chanIdx = ismember(data_tla.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,roi{r})})));
            else
              chanIdx = ismember(data_tla.(exper.sesStr{ses}).(cond_str).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,roi{r})})));
            end
            
            for lat = 1:length(latency)
              if ~isempty(cond_str_tmp)
                tbeg = nearest(data_tla.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.time,latencies(lat,1));
                tend = nearest(data_tla.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.time,latencies(lat,2));
              else
                tbeg = nearest(data_tla.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,1));
                tend = nearest(data_tla.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,2));
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
                thisData = (mean(mean(mean(data_tla.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.(measure)(chanIdx,latIdx),3),2),1) + ...
                  mean(mean(mean(data_tla.(exper.sesStr{ses}).(cond_str_tmp{2}).sub(sub).data.(measure)(chanIdx,latIdx),3),2),1)) ./ 2;
              else
                thisData = mean(mean(mean(data_tla.(exper.sesStr{ses}).(cond_str).sub(sub).data.(measure)(chanIdx,latIdx),3),2),1);
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
fprintf('This ANOVA: %s\n\n',stimType);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: %s\n',stimType);
fprintf('=======================================\n');

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
% and the times that correspond to each set of ROIs

%cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'},{'LPS','RPS'},{'LPS','RPS'},{'LAS','RAS'},{'LPS','RPS'},{'LPS','RPS'},{'LPS','RPS'}};
%cfg_ana.latencies = [0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];

% cfg_ana.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_ana.latencies = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];

cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_ana.rois = {{'FC'}};
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% cfg_ana.rois = {{'LPI2'}};
% % cfg_ana.rois = {{'E70','E71','E66','E75','E83','E76','E84'}};
% cfg_ana.latencies = [0.14 0.2];

% cfg_ana.rois = {{'PI'}};
% % cfg_ana.rois = {{'E70','E71','E66','E75','E83','E76','E84'}};
% cfg_ana.latencies = [0.13 0.19];
% % cfg_ana.latencies = [0.14 0.19];

% cfg_ana.rois = {{'E70','E71','E66'},{'E83','E76','E84'}};
% % cfg_ana.rois = {{'E70','E66','E65','E69'},{'E83','E84','E89','E90'}};
% % cfg_ana.rois = {{'E70'},{'E83'}};
% cfg_ana.latencies = [0.14 0.19; 0.14 0.19];

% % LF O/N
% cfg_ana.rois = {{'RAS'},{'RAS'},{'RAI'},{'RAI'}};
% cfg_ana.latencies = [1.1 1.4; 1.4 1.9; 1.1 1.4; 1.4 1.9];

% % LPN
% cfg_ana.rois = {{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_ana.latencies = [1.2 1.8; 1.2 1.8; 1.2 1.8];

cfg_ana.conditions = {{'word_RgH_rc_spac_p2','word_onePres'},{'word_RgH_rc_mass_p2','word_onePres'},{'word_RgH_fo_spac_p2','word_onePres'},{'word_RgH_fo_mass_p2','word_onePres'}};
% cfg_ana.conditions = {{'word_RgH_rc_spac_p2','word_RgH_rc_mass_p2'},{'word_RgH_fo_spac_p2','word_RgH_fo_mass_p2'}};

%cfg_ana.conditions = {{'all'}};
%cfg_ana.conditions = {{'all_within_types'}};
%cfg_ana.conditions = {{'all_across_types'}};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.parameter = 'avg';
cfg_ft.correctm = 'fdr';

% line plot parameters
cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
%cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 1 4; 1 4; 1 4; -4 -1; 1 4; 1 4; 1 4];
%cfg_plot.ylims = [-4 -1; 1 4; 1 4];
% cfg_plot.ylims = [-5 -2; -5 -2; -5 -2; 1.5 4.5; 1.5 4.5;];
cfg_plot.ylims = [-5 -2; 1.5 4.5];
%cfg_plot.plot_order = {'CR2','H2','HSC2','HSI2','CR6','H6','HSC6','HSI6'};
%cfg_plot.plot_order = {'RHSC','RHSI','RCR'};
% cfg_plot.plot_order = {'SC','SI','CR'};
%cfg_plot.rename_conditions = {'SC','SI','CR'};
% cfg_plot.xlabel = 'Condition';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

sesNum = 1;

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  if cfg_plot.individ_plots || cfg_plot.line_plots
    cfg_plot.ylim = cfg_plot.ylims(r,:);
  end
  
  % use data_tla_avg because FieldTrip doesn't deal with dimord properly
  % when single-trials exist
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla_avg,sesNum);
  
%   mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,ga_tla,sesNum);
end

%% output some values

cfg = [];

cfg.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg.latencies = [0.3 0.5; 0.5 0.8];
% cfg.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg.latencies = [0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];
cfg.condByROI = repmat({{'SC','SI','CR'}},size(cfg.rois));

% cfg.rois = {{'FC'}};
% cfg.latencies = [0.3 0.5];
% cfg.condByROI = repmat({{'FSC','FSI','N'}},size(cfg.rois));

cfg.parameter = 'avg';

cfg.excludeBadSub = false;

%cfg.direction = 'columns';
cfg.direction = 'rows';

mm_printDataToText(cfg,exper,ana,dirs,data_tla);

%% 3-way ANOVA: Hemisphere x Block Type x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% IV2: abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% IV3: outermost cell holds one cell for each ROI; each ROI cell holds one
% cell for each event type; each event type cell holds strings for its
% conditions
cfg_ana.condByTypeByROI = {...
  %{{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
  {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}},...
  {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

% For each ROI, what's common among the conditions in each type
cfg_ana.condCommonByROI = {...
  %{'CR','H','HSC','HSI'},...
  {'CR','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Block Type','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByTypeByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov33ER(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Hemisphere x Condition - separate colors

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% % IV2: define the conditions tested for each set of ROIs
% cfg_ana.condByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% 
% % For each ROI, what's common among the conditions in each type
% cfg_ana.condCommonByROI = {...
%   {'CR','H','HSC','HSI'},...
%   {'CR','HSC','HSI'}};
% 
% % abbreviations for the condition types
% cfg_ana.typesByROI = {...
%   {'C2','C6'},...
%   {'C2','C6'}};

% IV2: define the conditions tested for each set of ROIs
%cfg_ana.condByROI = {{'RCR','RH'},{'RCR','RHSC','RHSI'}};
cfg_ana.condByROI = {{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
% for the one-way ANOVAs
cfg_ana.condCommonByROI = {...
  %{'CR','H'},...
  {'CR','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  %cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% IV2: define the conditions tested for each set of ROIs
%cfg_ana.condByROI = {{'RH','RCR'},{'RCR','RHSC','RHSI'}};
%cfg_ana.condByROI = {{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
cfg_ana.condByROI = {exper.eventValues,exper.eventValues};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% cfg_ana.condCommonByROI = {...
%   {'CR','HSC','HSI'},...
%   {'CR','HSC','HSI'}};
cfg_ana.condCommonByROI = {...
  exper.eventValues,...
  exper.eventValues};

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% 1-way ANOVA: Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 0;

% define which regions to average across for the test
cfg_ana.rois = {{'FC'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5];

cfg_ana.condByROI = repmat({{'FSC', 'FSI', 'N'}},size(cfg_ana.rois));
% Define the IVs (type: event, roi, latency)
cfg_ana.IV1.name = 'FSC/FSI/CR';
cfg_ana.IV1.cond = {'FSC', 'FSI', 'N'};
cfg_ana.IV1.type = 'event';

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  
  mm_ft_rmaov1ER_spec(cfg_ana,exper,ana,data_tla);
end

%% cluster statistics

cfg_ft = [];
cfg_ft.avgovertime = 'no';
% cfg_ft.avgovertime = 'yes';
% cfg_ft.avgoverchan = 'no';
cfg_ft.avgoverchan = 'yes';

cfg_ft.parameter = 'avg';

cfg_ana = [];
% cfg_ana.roi = 'all';
% cfg_ana.roi = {'center109'};
cfg_ana.roi = {'LPS'};
%cfg_ana.latencies = [0 1.0; 1.0 2.0];
%cfg_ana.latencies = [0.2 0.6];
% cfg_ana.latencies = [0.3 0.5];
% cfg_ana.latencies = [0.5 0.8];
cfg_ana.latencies = [0 1.0];

cfg_ft.numrandomization = 500;
% cfg_ft.clusteralpha = 0.05;
cfg_ft.clusteralpha = 0.1;
cfg_ft.alpha = 0.05;

% extra directory info
cfg_ana.dirStr = '';
if strcmp(cfg_ft.avgovertime,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end
if strcmp(cfg_ft.avgoverchan,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgC'];
end

% cfg_ana.conditions = {...
%   {'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},...
%   {'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'},...
%   {'CR2','CR6'},{'H2','H6'},{'HSC2','HSC6'},{'HSI2','HSI6'}};
% cfg_ana.conditions = {'all_within_types'};
%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};

cfg_ana.conditions = {{'word_RgH_rc_spac_p2','word_onePres'},{'word_RgH_fo_spac_p2','word_onePres'},{'word_RgH_rc_mass_p2','word_onePres'},{'word_RgH_fo_mass_p2','word_onePres'}};
cfg_ana.conditions = {{'word_RgH_rc_spac_p2','word_onePres'}};

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  
  stat_clus = mm_ft_clusterstatER(cfg_ft,cfg_ana,exper,ana,dirs,data_tla);
end

%% plot the cluster statistics

%files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = 0.1;

cfg_plot = [];
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.dirStr = cfg_ana.dirStr;

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  
  mm_ft_clusterplotER(cfg_ft,cfg_plot,ana,files,dirs);
end

%% let me know that it's done
emailme = 1;
if emailme
  subject = sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    };
  send_gmail(subject,mail_message);
end

%% correlations

cfg_ana = [];

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.dpTypesByROI = {...
  {'Item','Source'},...
  {'Item','Source'}};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds two cells, one for each
% d' type; each d' cell contains strings for its conditions
cfg_ana.condByROI = {...
  {{{'CR2','H2'},{'HSC2','HSI2'}}, {{'CR6','H6'},{'HSC6','HSI6'}}}...
  {{{'CR2','H2'},{'HSC2','HSI2'}}, {{'CR6','H6'},{'HSC6','HSI6'}}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% C2 d' values
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_item =  abs([]);
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_source =  abs([]);
% C6 d' values
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_item =  abs([]);
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_source =  abs([]);

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data_tla);
end
