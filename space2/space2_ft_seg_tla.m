% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SPACE2';

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
% exper.eegFileExt = 'raw';
exper.eegFileExt = 'mff';

% types of events to find in the NS file; these must be the same as the
% events in the NS files; or space_trialfun.m must be set up to find the
% corrct events
exper.eventValues = {{'multistudy_image', 'multistudy_word', 'cued_recall_stim'}};
% exper.eventValues = {{'multistudy_image', 'multistudy_word', 'cued_recall_stim'}, {'multistudy_image', 'multistudy_word', 'cued_recall_stim'}};
% exper.eventValues = {{'expo_stim'}};
% exper.eventValues = {{'expo_stim'},{'expo_stim'}};

% pre- and post-stimulus times to read, in seconds (pre is negative).
% Construct as a cell with one Nx2 matrix per session where N is
% length(exper.eventValues{ses}) Order must correspond to the event order
% in exper.eventValues.
exper.prepost = {[-1.0 2.0; -1.0 2.0; -1.0 2.0]};
% exper.prepost = {[-1.0 2.0; -1.0 2.0; -1.0 2.0], [-1.0 2.0; -1.0 2.0; -1.0 2.0]};
% exper.prepost = {[-1.0 2.0]};
% exper.prepost = {[-0.2 1.0], [-0.2 1.0]};

exper.subjects = {
  'SPACE2001';
%   'SPACE2002';
%   'SPACE2003';
%   'SPACE2004';
%   'SPACE2005';
%   'SPACE2006';
%   'SPACE2007';
%   'SPACE2008';
%   'SPACE2009';
%   'SPACE2010';
%   'SPACE2011';
%   'SPACE2012';
%   'SPACE2013';
%   'SPACE2014';
%   'SPACE2015';
%   'SPACE2016';
%   'SPACE2017';
%   'SPACE2018';
%   'SPACE2019';
%   'SPACE2020';
  };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, other directories (ns_evt, ns_bci, etc).
% They are not necessarily the session directory names where the FieldTrip
% data is saved for each subject because of the option to combine sessions.
% See 'help create_ft_struct' for more information.
% exper.sessions = {{'session_1'}};
exper.sessions = {{'session_1', 'session_2'}};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.subDir = '';
% dirs.dataDir = fullfile(exper.name,'EEG','Sessions','face_house_ratings','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
dirs.behDir = fullfile(exper.name,'Behavioral','Sessions',dirs.subDir);
% dirs.dataDir = fullfile(exper.name,'EEG','Sessions','ftpp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
dirs.dataDir = fullfile(exper.name,'EEG','Sessions','ftpp',dirs.subDir);
% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'data');

% pick the right dirs.dataroot
if isfield(dirs,'serverDir') && exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
  %runLocally = 1;
elseif isfield(dirs,'serverLocalDir') && exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
  %runLocally = 1;
elseif isfield(dirs,'dreamDir') && exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
  %runLocally = 0;
elseif isfield(dirs,'localDir') && exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

% Use the FT chan locs file
files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
exper.refChan = 'Cz';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';

ana.offsetMS = 36;

ana.continuous = 'yes';
% ana.trialFxn = 'space_trialfun';
ana.trialFxn = 'space2_trialfun_mff';
ana.allowTrialOverlap = true;
ana.renumberSamplesContiguous = true;
% files used when adding metadata to segmented trials
ana.useMetadata = true;
ana.metadata.types = {'eventStruct','nsEvt'};
ana.useExpInfo = true;
% ana.evtToleranceMS = 8; % 2 samples @ 250 Hz
ana.usePhotodiodeDIN = true;
ana.photodiodeDIN_toleranceMS = 20;
ana.photodiodeDIN_str = 'DIN ';
if ana.useExpInfo
  % possible sessions and phases
  %ana.sessionNames = {'oneDay'};
  ana.sessionNames = {'day1','day2'};
  
  % phases occur within a session; for dealing with events.mat
  ana.phaseNames = {{'multistudy', 'cued_recall_only'}};
  %ana.phaseNames = {{'multistudy', 'cued_recall_only'}, {'multistudy', 'cued_recall_only'}};
  %ana.phaseNames = {{'expo'}};
  %ana.phaseNames = {{'expo'},{'expo'}};
  %ana.phaseNames = {{'multistudy'}};
  %ana.phaseNames = {{'distract_math'}};
  %ana.phaseNames = {{'cued_recall'}};
  
  % types of event info to store in trialinfo field; must correspond to
  % values listed in exper.eventValues
  % ana.trl_order.expo_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'expo_response', 'rt', 'cr_recog_acc', 'cr_recall_resp', 'cr_recall_spellCorr'};
  ana.trl_order.multistudy_image = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'catNum', 'targ', 'spaced', 'lag', 'presNum', 'pairOrd', 'pairNum', 'cr_recall_resp', 'cr_recall_spellCorr'};
  ana.trl_order.multistudy_word = ana.trl_order.multistudy_image;
  %ana.trl_order.distract_math_stim = {'eventNumber', 'sesType', 'phaseType', 'response', 'acc', 'rt'};
  %ana.trl_order.cued_recall_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'pairNum', 'recog_resp', 'recog_acc', 'recog_rt', 'new_resp', 'new_acc', 'new_rt', 'recall_resp', 'recall_spellCorr', 'recall_rt'};
  ana.trl_order.cued_recall_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'pairNum', 'recall_resp', 'recall_spellCorr', 'recall_rt'};
end

% preprocess continuous data in these ways
ana.cfg_cont.lpfilter = 'yes';
ana.cfg_cont.lpfreq = 100;
ana.cfg_cont.hpfilter = 'yes';
ana.cfg_cont.hpfreq = 0.1;
ana.cfg_cont.hpfilttype = 'but';
ana.cfg_cont.hpfiltord = 4;
ana.cfg_cont.bsfilter = 'yes';
ana.cfg_cont.bsfreq = [59 61];

% ana.artifact.continuousRepair = true;
% ana.artifact.continuousReject = true;
% ana.artifact.continuousICA = true;

% this will turn off ana.cfg_cont and ana.artifact.continuousXXX so no
% processing is done; will only save FieldTrip events
ana.returnAfterSavingFtEvents = false;

% artifact settings
ana.artifact.reject = 'complete';
ana.artifact.preArtBaseline = 'yes'; % yes=entire trial

% % ana.artifact.type = {'nsClassic','ftAuto'};
% ana.artifact.type = {'nsClassic'};
% 
% % set up for nsClassic
% ana.artifact.checkArtSec = [0 1.0];
% % ana.artifact.blink_threshold = 70;
% % ana.artifact.fast_threshold = 100;
% % ana.artifact.diff_threshold = 50;
% ana.artifact.blink_threshold = 40;
% ana.artifact.fast_threshold = 70;
% ana.artifact.diff_threshold = 40;
% ana.artifact.rejectTrial_nBadChan = 10;
% ana.artifact.repairChan_percentBadTrials = 20;
% ana.artifact.allowBadNeighborChan = false;

ana.artifact.type = {'ftAuto'};
% set up for ftAuto following nsClassic
% negative trlpadding: don't check that time (on both sides) for artifacts.
% IMPORTANT: Not used for threshold artifacts. only use if segmenting a lot
% of extra time around trial epochs. Otherwise set to zero.
ana.artifact.trlpadding = -1.0;
% ana.artifact.trlpadding = 0;
ana.artifact.artpadding = 0.1;
ana.artifact.fltpadding = 0;

ana.artifact.thresh = true;
ana.artifact.threshmin = -100;
ana.artifact.threshmax = 100;
ana.artifact.threshrange = 150;
ana.artifact.basic_art = true;
ana.artifact.basic_art_z = 30;
ana.artifact.jump_art = true;
ana.artifact.jump_art_z = 50;
% eog_art is only used with ftAuto
ana.artifact.eog_art = false;
% ana.artifact.eog_art_z = 3.5;

% %%%%%
% % START OF FTMANUAL+FTICA ARTIFACT CHECKING
% %%%%%
% 
% ana.artifact.type = {'ftManual', 'ftICA'};
% % ana.artifact.type = {'ftAuto', 'ftICA'};
% % ana.artifact.type = {'ftAuto'};
% ana.artifact.resumeManArtFT = false;
% ana.artifact.resumeICACompFT = false;
% % negative trlpadding: don't check that time (on both sides) for artifacts.
% % IMPORTANT: Not used for threshold artifacts. only use if segmenting a lot
% % of extra time around trial epochs. Otherwise set to zero.
% ana.artifact.trlpadding = -0.7;
% % ana.artifact.trlpadding = 0;
% ana.artifact.artpadding = 0.1;
% ana.artifact.fltpadding = 0;
% 
% % set up for ftManual/ftAuto
% ana.artifact.thresh = true;
% % ana.artifact.threshmin = -150;
% % ana.artifact.threshmax = 150;
% % ana.artifact.threshrange = 250;
% ana.artifact.threshmin = -200;
% ana.artifact.threshmax = 200;
% ana.artifact.threshrange = 350;
% ana.artifact.basic_art = true;
% ana.artifact.basic_art_z = 60;
% ana.artifact.jump_art = true;
% ana.artifact.jump_art_z = 70;
% % % eog_art is only used with ftAuto
% % ana.artifact.eog_art = false;
% % ana.artifact.eog_art_z = 3.5;
% 
% % set up for ftICA
% ana.artifact.thresh_postICA = true;
% ana.artifact.threshmin_postICA = -100;
% ana.artifact.threshmax_postICA = 100;
% ana.artifact.threshrange_postICA = 150;
% ana.artifact.basic_art_postICA = true;
% ana.artifact.basic_art_z_postICA = 30;
% ana.artifact.jump_art_postICA = true;
% ana.artifact.jump_art_z_postICA = 50;
% 
% %%%%%
% % END OF FTMANUAL+FTICA ARTIFACT CHECKING
% %%%%%


% %%%%%
% % USE THIS: START OF POST-CONTINUOUS ARTIFACT CHECKING
% %%%%%
% 
% % set up for ftManual/ftAuto
% % ana.artifact.type = {'ftManual'};
% ana.artifact.type = {'ftAuto'};
% % ana.artifact.resumeManArtFT = false;
% ana.artifact.resumeManArtContinuous = false;
% ana.artifact.resumeICACompContinuous = false;
% % negative trlpadding: don't check that time (on both sides) for artifacts.
% % IMPORTANT: Not used for threshold artifacts. only use if segmenting a lot
% % of extra time around trial epochs. Otherwise set to zero.
% ana.artifact.trlpadding = -1.0;
% % ana.artifact.trlpadding = 0;
% ana.artifact.artpadding = 0.1;
% ana.artifact.fltpadding = 0;
% 
% % set up for ftAuto after continuous ICA rejection
% ana.artifact.checkAllChan = true;
% ana.artifact.thresh = true;
% ana.artifact.threshmin = -100;
% ana.artifact.threshmax = 100;
% ana.artifact.threshrange = 150;
% ana.artifact.basic_art = true;
% ana.artifact.basic_art_z = 30;
% ana.artifact.jump_art = true;
% ana.artifact.jump_art_z = 30;
% % eog_art is only used with ftAuto
% ana.artifact.eog_art = false;
% % ana.artifact.eog_art_z = 3.5;
% 
% %%%%%
% % END OF POST-CONTINUOUS ARTIFACT CHECKING
% %%%%%

% process the data
ana.ftFxn = 'ft_timelockanalysis';
% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'tla';
ana.overwrite.raw = 1;
ana.overwrite.proc = 1;

% any preprocessing? (run after processing artifacts)
cfg_pp = [];
% average rereference
cfg_pp.reref = 'yes';
cfg_pp.refchannel = 'all';
cfg_pp.implicitref = exper.refChan;
% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];
% single precision to save space
%cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.keeptrials = 'yes';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs_multiSes(exper,ana,cfg_proc,dirs,files,'tla',true);
% [dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct_multiSes(ana,cfg_pp,exper,dirs,files);
% [exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);

process_ft_data_multiSes(ana,cfg_proc,exper,dirs,files,cfg_pp);
% process_ft_data(ana,cfg_proc,exper,dirs);

% %% get the bad channel information
% 
% cfg = [];
% cfg.badChanManual = false;
% cfg.badChanEP = true;
% [exper] = mm_getBadChan(cfg,exper,dirs);

% % save the analysis details
% 
% backup_orig_AD = true;
% % whether to sort by subject number
% sortBySubj = true;
% % whether to overwite existing subjects in the struct
% replaceOrig = true;
% 
% % concatenate the additional ones onto the existing ones
% saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
% if ~exist(saveFile,'file')
%   fprintf('Saving analysis details: %s...',saveFile);
%   save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
%   fprintf('Done.\n');
% else
%   additional_AD_file = fullfile(dirs.saveDirProc,sprintf('analysisDetails%s.mat',sprintf(repmat('_%s',1,length(exper.subjects)),exper.subjects{:})));
%   fprintf('Temporarily saving new analysis details: %s...',additional_AD_file);
%   save(additional_AD_file,'exper','ana','dirs','files','cfg_proc','cfg_pp');
%   
%   [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails(saveFile,additional_AD_file,backup_orig_AD,sortBySubj,replaceOrig);
% end
% 
% %% let me know that it's done
% emailme = 0;
% if emailme
%   subject = sprintf('Done with%s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
%   mail_message = {...
%     sprintf('Done with%s %s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
%     sprintf('%s',saveFile),...
%     };
%   send_gmail(subject,mail_message);
% end
