% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SPACE';

exper.sampleRate = 250;

% % pre- and post-stimulus times to read, in seconds (pre is negative)
% exper.prepost = [-1.0 1.0];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files; or space_trialfun.m must be set up to find the
% corrct events
% % [exper.eventValues, evInd] = sort({'expo_face', 'expo_house'});
% [exper.eventValues, evInd] = sort({'expo_stim'});
% [exper.eventValues, evInd] = sort({'multistudy_image', 'multistudy_word'});
% [exper.eventValues, evInd] = sort({'distract_math_stim'});
% [exper.eventValues, evInd] = sort({'cued_recall_stim'});
[exper.eventValues, evInd] = sort({'expo_stim', 'multistudy_image', 'multistudy_word', 'cued_recall_stim'});

% pre- and post-stimulus times to read, in seconds (pre is negative);
% because they get sorted, must correspond to the order listed in
% exper.eventValues
exper.prepost = [...
  -1.0 2.0; ...
  -1.0 2.0; ...
  -1.0 2.0; ...
  -1.0 2.0];
% exper.prepost = [-1.0 2.0];
% exper.prepost = [-0.2 1.0];
exper.prepost = exper.prepost(evInd,:);

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 0;
exper.eventValuesExtra.equateExtrasSeparately = 0;

exper.subjects = {
  'SPACE001';
  'SPACE002';
  'SPACE003';
  'SPACE004';
  'SPACE005';
  'SPACE006';
  'SPACE007';
  %'SPACE008';
  'SPACE009';
  'SPACE010';
  'SPACE011';
  'SPACE012';
  'SPACE013';
  'SPACE014';
  'SPACE015';
  };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_1'};

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
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';

ana.continuous = 'yes';
ana.trialFxn = 'space_trialfun';
% files used when adding metadata to segmented trials
ana.useMetadata = true;
ana.metadata.types = {'eventStruct','nsEvt'};
ana.useExpInfo = true;
ana.usePhotodiodeDIN = true;
ana.photodiodeDIN_thresholdMS = 50;
ana.photodiodeDIN_str = 'DIN ';
if ana.useExpInfo
  % possible sessions and phases
  ana.sessionNames = {'oneDay'};
  
  % phases occur within a session; for dealing with events.mat
  ana.phaseNames = {{'expo', 'multistudy', 'cued_recall'}};
  %ana.phaseNames = {{'expo'}};
  %ana.phaseNames = {{'multistudy'}};
  %ana.phaseNames = {{'distract_math'}};
  %ana.phaseNames = {{'cued_recall'}};
  
  % types of event info to store in trialinfo field; must correspond to
  % values listed in exper.eventValues
  ana.trl_order.expo_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'expo_response', 'rt', 'cr_recog_acc', 'cr_recall_resp', 'cr_recall_spellCorr'};
  ana.trl_order.multistudy_image = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'catNum', 'targ', 'spaced', 'lag', 'presNum', 'pairOrd', 'pairNum', 'cr_recog_acc', 'cr_recall_resp', 'cr_recall_spellCorr'};
  ana.trl_order.multistudy_word = ana.trl_order.multistudy_image;
  %ana.trl_order.distract_math_stim = {'eventNumber', 'sesType', 'phaseType', 'response', 'acc', 'rt'};
  ana.trl_order.cued_recall_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'pairNum', 'recog_resp', 'recog_acc', 'recog_rt', 'new_resp', 'new_acc', 'new_rt', 'recall_resp', 'recall_spellCorr', 'recall_rt'};
end

% preprocess continuous data in these ways
ana.cfg_cont.lpfilter = 'yes';
ana.cfg_cont.lpfreq = 100;
ana.cfg_cont.hpfilter = 'yes';
ana.cfg_cont.hpfreq = 0.1;
ana.cfg_cont.hpfilttype = 'but';
ana.cfg_cont.hpfiltord = 4;
ana.cfg_cont.bsfilter = 'yes';
ana.cfg_cont.bsfreq = 59:61;

% artifact settings
ana.artifact.type = {'ftManual', 'ftICA'};
ana.artifact.resumeManArtFT = false;
ana.artifact.resumeICACompFT = false;
% negative trlpadding: don't check that time (on both sides) for artifacts
ana.artifact.trlpadding = -0.5;
% ana.artifact.trlpadding = 0;
ana.artifact.artpadding = 0.1;
ana.artifact.fltpadding = 0;
ana.artifact.threshmin = -150;
ana.artifact.threshmax = 150;
ana.artifact.basic_art_z = 30;
ana.artifact.muscle_art_z = 40;
ana.artifact.jump_art_z = 50;
ana.artifact.threshmin_postICA = -150;
ana.artifact.threshmax_postICA = 150;
ana.artifact.basic_art_z_postICA = 30;
ana.artifact.muscle_art_z_postICA = 40;
ana.artifact.jump_art_z_postICA = 50;
ana.overwrite.raw = 1;

% process the data
ana.ftFxn = 'ft_timelockanalysis';
% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'tla';
ana.overwrite.proc = 1;

% any preprocessing? (run after processing artifacts)
cfg_pp = [];
% average rereference
cfg_pp.reref = 'yes';
cfg_pp.refchannel = 'all';
cfg_pp.implicitref = 'Cz';
% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];
% single precision to save space
%cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.keeptrials = 'yes';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);
process_ft_data(ana,cfg_proc,exper,dirs);

% %% get the bad channel information
% 
% cfg = [];
% cfg.badChanManual = false;
% cfg.badChanEP = true;
% [exper] = mm_getBadChan(cfg,exper,dirs);

% save the analysis details

backup_orig_AD = true;
% whether to sort by subject number
sortBySubj = true;
% whether to overwite existing subjects in the struct
replaceOrig = true;

% concatenate the additional ones onto the existing ones
saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
if ~exist(saveFile,'file')
  fprintf('Saving analysis details: %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  fprintf('Done.\n');
else
  additional_AD_file = fullfile(dirs.saveDirProc,sprintf('analysisDetails%s.mat',sprintf(repmat('_%s',1,length(exper.subjects)),exper.subjects{:})));
  fprintf('Temporarily saving new analysis details: %s...',additional_AD_file);
  save(additional_AD_file,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  
  [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails(saveFile,additional_AD_file,backup_orig_AD,sortBySubj,replaceOrig);
end

%% let me know that it's done
emailme = 0;
if emailme
  subject = sprintf('Done with%s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with%s %s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    sprintf('%s',saveFile),...
    };
  send_gmail(subject,mail_message);
end
