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
procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');

subjects = {
  %'SPACE001'; low trial counts
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
  %'SPACE019'; low trial counts
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
  %'SPACE030'; low trial counts
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
% cfg.baseline_type = 'absolute';
% cfg.baseline_type = 'relchange';
% cfg.baseline_type = 'relative';
% cfg.baseline_type = 'db';

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
if strcmp(cfg.rmevoked,'yes') && ~exist('data_evoked','var')
  % load('/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0_art_zeroVar/tla_-1000_2000_avg/data_evoked.mat');
  
  % local testing
  %load('/Users/matt/data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar_badChanManual_badChanEP/tla_-1000_2000_avg/data_evoked.mat');
end

if isfield(cfg,'equatetrials') && strcmp(cfg.equatetrials,'yes')
  eq_str = '_eq';
else
  eq_str = '';
end
if isfield(cfg,'keeptrials') && strcmp(cfg.keeptrials,'yes')
  kt_str = '_trials';
else
  kt_str = '_avg';
end
if isfield(cfg,'rmevoked') && strcmp(cfg.rmevoked,'yes')
  indu_str = '_induced';
else
  indu_str = '_whole';
end
saveFile = fullfile(dirs.saveDirProc,sprintf('data_%s%s%s%s.mat',cfg.output,eq_str,kt_str,indu_str));

if exist(saveFile,'file')
  fprintf('Loading saved file: %s\n',saveFile);
  load(saveFile);
else
  fprintf('Running mm_ft_loadData_multiSes\n');
  if exist('data_evoked','var')
    [data_pow,exper] = mm_ft_loadData_multiSes(cfg,exper,dirs,ana,data_evoked);
  else
    [data_pow,exper] = mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
  end
  if cfg.saveFile
    fprintf('Saving %s...\n',saveFile);
    save(saveFile,sprintf('data_%s',cfg.output),'exper','cfg');
  end
end
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

%% decide who to kick out based on trial counts

% Subjects with bad behavior
% exper.badBehSub = {{}};
exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE039'}};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

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

% fieldtrip's baseline method
cfg_bl = [];
cfg_bl.baseline = [-0.3 -0.1];
cfg_bl.baselinetype = 'absolute';
for ses = 1:length(exper.sesStr)
  for typ = 1:length(ana.eventValues{ses})
    for evVal = 1:length(ana.eventValues{ses}{typ})
      for sub = 1:length(exper.subjects)
        data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_freqbaseline(cfg_bl,data_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
      end
    end
  end
end

% % my zscore method
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
%         blm = nanmean(subSesEvData.(cfg.param)(:,:,:,blt),4);
%         
%         if strcmp(cfg.baseline_type,'zscore')
%           fprintf('Z-transforming data relative to mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline_time(1),cfg.baseline_time(2));
%           % std across time, then avg across events (lower freqs often get smaller std)
%           blstd = nanmean(nanstd(subSesEvData.(cell2mat(fn)).(cfg.param)(:,:,:,blt),0,4),1);
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

%% save

% saveDir = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
% saveDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
saveDir = dirs.saveDirProc;
% save(fullfile(saveDir,'space_word_data_ga_pow.mat'),'data_pow','ga_pow','exper','ana','dirs','files','-v7.3');
save(fullfile(saveDir,'space_word_img_data_ga_pow.mat'),'data_pow','ga_pow','exper','ana','dirs','files','-v7.3');
% clear data_pow

%% load

loadDir = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
% loadDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow';
% load(fullfile(loadDir,'space_word_data_ga_pow.mat'));
load(fullfile(loadDir,'space_word_img_data_ga_pow.mat'));

[dirs] = mm_checkDirs(dirs);

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
cfg_ft.zlim = [-2 2];
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
cfg_plot.rois = {{'RPI2'}};
%cfg_plot.rois = {'E124'};
%cfg_plot.rois = {'E25'};
%cfg_plot.rois = {'RAS'};
%cfg_plot.rois = {'LPS','RPS'};
%cfg_plot.rois = {'LPS'};
cfg_plot.excludeBadSub = 1;
cfg_plot.numCols = 5;

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
ses=1;
% cfg_plot.condByROI = repmat(ana.eventValues{ses},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{{{'word_RgH_rc_spac_p2','word_RgH_rc_mass_p2'}}}},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{{{'word_RgH_fo_spac_p2'}}}},size(cfg_plot.rois));

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

%% line plots

files.saveFigs = true;
files.figPrintFormat = 'png';

cfg = [];
cfg.parameter = 'powspctrm';

%cfg.times = [-0.2:0.05:0.9; -0.1:0.05:1.0]';
%cfg.times = [-0.2:0.1:0.9; -0.1:0.1:1.0]';
% cfg.times = [-0.2:0.2:0.8; 0:0.2:1.0]';

%cfg.times = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
cfg.times = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap

% cfg.freqs = [3 8; 8 12; 12 28; 28 50; 50 80];
% cfg.freqs = [3 8; 8 10; 10 12];

% cfg.freqs = [3 7; 8 12; 8 10; 11 12; 13 20; 23 30; 32 47; 51 80];
cfg.freqs = [3 7; 8 12; 13 20; 23 30; 32 47; 51 80];

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
cfg.conditions = {{'word_RgH_rc_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_spac_p2','word_RgH_fo_mass_p2'}};

cfg.plotTitle = true;
cfg.plotLegend = true;

cfg.plotClusSig = false;
cfg.clusAlpha = 0.1;
%cfg.clusTimes = cfg.times;
% cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
%cfg.clusTimes = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
cfg.clusTimes = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
cfg.clusLimits = false;

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

%% plot the conditions - simple

cfg_ft = [];
cfg_ft.xlim = [-0.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

% %cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

cfg_plot.rois = {{'LAS'},{'RAS'}};
cfg_plot.ylims = [-4 2; -4 2];
cfg_plot.rois = {{'FC'},{'FS'}};
cfg_plot.ylims = [-4 2; -4 2];
% cfg_plot.rois = {{'posterior'}};
cfg_plot.rois = {{'LPS'},{'RPS'}};
cfg_plot.ylims = [-2 3; -2 3];
cfg_plot.legendlocs = {'SouthEast','NorthWest'};

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% cfg_ft.xlim = [-0.2 1.0];
% cfg_plot.rois = {{'E70'},{'E83'}};
% cfg_plot.ylims = [-10 10; -10 10];
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

cfg_plot.condByROI = repmat({{'word_RgH_rc_spac_p1', 'word_RgH_rc_spac_p2', 'word_RgH_fo_spac_p1', 'word_RgH_fo_spac_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_rc_spac_p1', 'img_RgH_rc_spac_p2', 'img_RgH_fo_spac_p1', 'img_RgH_fo_spac_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'word_RgH_rc_mass_p1', 'word_RgH_rc_mass_p2', 'word_RgH_fo_mass_p1', 'word_RgH_fo_mass_p2'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_rc_mass_p1', 'img_RgH_rc_mass_p2', 'img_RgH_fo_mass_p1', 'img_RgH_fo_mass_p2'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'word_spac_p2', 'word_mass_p2', 'word_onePres'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_spac_p2', 'img_mass_p2', 'img_onePres'}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'word_RgH_spac_p2', 'word_RgH_mass_p2', 'word_RgH_spac_p1', 'word_RgH_mass_p1', 'word_onePres'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'img_RgH_spac_p2', 'img_RgH_mass_p2', 'img_RgH_spac_p1', 'img_RgH_mass_p1', 'img_onePres'}},size(cfg_plot.rois));


for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,sesNum,ga_pow);
  %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
end

%% subplots of each subject's ERPs

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.rois = {{'LAS'},{'LPS'}};
% cfg_plot.rois = {{'E70'},{'E83'}};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-0.2 1.0];
cfg_plot.ylim = [-10 10];
cfg_plot.parameter = 'avg';

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

%%

% a = ones(2,2,4);
a = ones(2,4);
for i = 1:size(a,1)
  for j = 1:size(a,2)
    a(i,j,:) = j * a(i,j,:);
  end
end

dim = size(a);
b = reshape(a, dim(1), prod(dim(2:end)));

%% make some GA plots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [-0.2 1.0]; % time
cfg_ft.xlim = [-0.2 1.5]; % time
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
cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5; -2 5];
cfg_plot.rois = {{'FC'}};
cfg_plot.ylims = [-7.5 2];
cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.plotLegend = 0;
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};

% cfg_plot.xlabel = 'Time (s)';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

% cfg_plot.ftFxn = 'ft_topoplotER';
% cfg_plot.ylims = [-6 6];
% %cfg_plot.ylims = 'maxmin';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
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
cfg_plot.condByROI = repmat({{{'FSC','FSI','N'}}},size(cfg_plot.rois));
cfg_plot.rename_condByROI = repmat({{{'FSC','FSI','CR'}}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));
%cfg_plot.rename_condByROI = repmat({{{'SC','SI','CR'}}},size(cfg_plot.rois));

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

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.types = cfg_plot.typesByROI{r};
  cfg_plot.rename_conditions = cfg_plot.rename_condByROI{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  
  if strcmp(cfg_plot.ftFxn,'ft_singleplotER')
    cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
    if cfg_plot.plotLegend
      cfg_plot.legendloc = cfg_plot.legendlocs{r};
    end
  end
  
  mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,ga_pow);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_ft = [];
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
cfg_ft.colormap = 'jet';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'no';

%cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';
%cfg_plot.conditions = {{'HSC2','CR2'},{'HSI2','CR2'},{'HSC2','HSI2'},{'HSC6','CR6'},{'HSI6','CR6'},{'HSC6','HSI6'}}; % {'H2','CR2'}, {'H6','CR6'},
%cfg_plot.conditions = {{'RH','RCR'},{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}};
% cfg_plot.conditions = {{'SC','CR'},{'SC','SI'},{'SI','CR'}};
%cfg_plot.conditions = {{'RHSC','RHSI'}};
%cfg_plot.conditions = {{'FSC','RSSI'}};
cfg_plot.conditions = {{'Face','House'}};


cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-5.5 5.5]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
% cfg_ft.comment = 'no';
cfg_ft.comment = 'xlim';
cfg_ft.commentpos = 'middletop';

cfg_plot.roi = {'E73'};
%cfg_plot.roi = {'LAS'};
% cfg_ft.xlim = [0.01 0.8]; % time

% cfg_plot.roi = {'LPS','RPS'};
% cfg_ft.xlim = [0.5 0.8]; % time

%cfg_plot.roi = {'RAS'};
%cfg_ft.xlim = [1.1 1.9]; % time
% cfg_ft.xlim = [0 1.5]; % time
% cfg_plot.roi = {'all'};

% cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time
cfg_ft.xlim = (0:0.05:1.0); % time

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

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_pow);

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
  cfg_ana.latencies = [0 0.48; 0.5 1.0];
elseif strcmp(cfg_ft.avgovertime,'yes')
  %cfg_ana.latencies = [-0.2:0.1:0.9; -0.1:0.1:1.0]'; % some overlap
  %cfg_ana.latencies = [-0.2:0.2:0.8; 0:0.2:1.0]'; % some overlap
  
  %cfg_ana.latencies = [-0.18:0.12:0.92; -0.1:0.12:1.0]'; % 100 no overlap
  %cfg_ana.latencies = [-0.18:0.2:0.82; -0.02:0.2:1.0]'; % 200 no overlap
  
  % new analyses
  %cfg_ana.latencies = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
  cfg_ana.latencies = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
  
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end

if strcmp(cfg_ft.avgoverfreq,'no')
  cfg_ana.frequencies = [4 100];
elseif strcmp(cfg_ft.avgoverfreq,'yes')
  % % cfg_ana.frequencies = [4 8; 8.1 14; 14.1 21; 21.1 28; 28.1 42; 42.1 64; 64.1 80];
  % % cfg_ana.frequencies = [4 8; 8 12; 12 28; 28 50; 50 100];
%   cfg_ana.frequencies = [3 8; 9 12; 9 10; 11 12; 13 28; 29 50; 51 80];
%   cfg_ana.frequencies = [3 7; 8 12; 8 10; 11 12; 13 20; 21 30; 31 45; 46 80]; % hanslmayr
  
  cfg_ana.frequencies = [3 7; 8 12; 8 10; 11 12; 13 20; 23 30; 32 47; 51 80];
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

cfg.plotClusSig = true;
% cfg.clusAlpha = 0.1;
cfg.clusAlpha = 0.05;
cfg.clusTimes = cfg.times;
% cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
cfg.clusLimits = true;

%cfg.ylim = [-0.6 0.6];
%cfg.ylim = [-0.5 0.2];
cfg.nCol = 3;

% whole power
cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-300_-100';
cfg.clusDirStr = '_avgT_avgF';
cfg.ylabel = 'Log10 Z-Trans Pow';
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

stimType = 'word_';
% stimType = 'img_';
memType = 'RgH_';

% spacings = {'mass', 'spac', 'onePres'};
% oldnew = {'p1'};
% % oldnew = {'p2'};
% memConds = {'all'};

% didn't test new words, so can't assess memory, but can use p1
spacings = {'mass', 'spac'};
% spacings = {'spac'};
oldnew = {'p1', 'p2'};
% oldnew = {'p1'};
% oldnew = {'p2'};
memConds = {'rc','fo'};

measure = 'powspctrm';

% % theta
freqs = [3 7];
% roi = {'LAS'};
% roi = {'PS'};
% roi = {'FC'};
% roi = {'FS'};
roi = {{'E23','Fz','E3'}}; % AF3 Fz AF4
% roi = {{'E67','Pz','E77'}}; % PO3 Pz P04
% roi = {'LPS'};
% roi = {'PS'};
% roi = {'RPS'}; % theta word P2 SME interaction 620-980
% roi = {'RPI2'};
% latencies = [0.3 1.0];
latencies = [0.62 0.98];
% latencies = [0.62 0.78];

% % % alpha
% % freqs = [8 10];
% % freqs = [11 12];
% freqs = [8 12];
% % roi = {'FC'}; % no
% % roi = {'FS'}; % maybe
% roi = {'FS2'}; % yes
% latencies = [0.4 0.6];

% freqs = [3 7; 8 12; 13 20; 21 30; 31 45; 46 80];

% latencies = [0.02 0.46; 0.5 0.98];
% latencies = [0.02 0.46];
% latencies = [0.5 0.98];

% freqs = [3 7];
% freqs = [8 12];
% freqs = [13 20];
% freqs = [21 30];
% freqs = [31 45];
% freqs = [46 80];

latency = cell(1,size(latencies,1));
for i = 1:length(latency)
  latency{i} = sprintf('%dto%d',latencies(i,1)*1000,latencies(i,2)*1000);
end
latStr = sprintf(repmat('_%s',1,length(latency)),latency{:});
latStr = latStr(2:end);

freq_str = sprintf('%dfreq%dto%d',size(freqs,1),freqs(1,1),freqs(end,end));
freqIdx = (nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,freqs(1,1)):nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,freqs(1,2)));

factorNames = {'spacings', 'oldnew', 'memConds', 'roi', 'latency'};

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
nSub = sum(~exper.badSub);
anovaData = nan(nSub,prod(nVariables));
rmaov_data_teg = nan(nSub*prod(nVariables),length(factorNames) + 2);

fprintf('Collecting ANOVA data for %d subjects:\n\t',nSub);
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
          if strcmp(spacings{sp},'onePres')
            % single presentation or first presentation
            if strcmp(memConds{mc},'all');
              cond_str = sprintf('%s%s_%s',stimType,spacings{sp});
            end
          elseif strcmp(spacings{sp},'mass') || strcmp(spacings{sp},'spac')
            cond_str = sprintf('%s%s%s_%s_%s',stimType,memType,memConds{mc},spacings{sp},oldnew{on});
          end
          
          for r = 1:length(roi)
            if iscell(roi{r})
              roi_str = sprintf(repmat('%s',1,length(roi{r})),roi{r}{:});
            elseif ischar(roi{r})
              roi_str = roi{r};
            end
            chanIdx = ismember(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,roi{r})})));
            
            for lat = 1:length(latency)
              latIdx = (nearest(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,1)):nearest(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,2)));
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
              
              anovaData(subCount,vnCount) = mean(mean(mean(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.(measure)(chanIdx,freqIdx,latIdx),3),2),1);
              
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
fprintf('This ANOVA: Freq: %s\n\n',freq_str);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: Freq: %s\n',freq_str);
fprintf('=======================================\n');

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
