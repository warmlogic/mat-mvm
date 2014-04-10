% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'EBIRD';

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files; or space_trialfun.m must be set up to find the
% corrct events
exper.eventValues = {{'match_stim'}};
% exper.eventValues = {{'match_stim'}, {'nametrain_stim', 'name_stim'}};

% pre- and post-stimulus times to read, in seconds (pre is negative).
% Construct as a cell with one Nx2 matrix per session where N is
% length(exper.eventValues{ses}) Order must correspond to the event order
% in exper.eventValues.
exper.prepost = {[-0.2 1.0]};
% exper.prepost = {[-0.2 1.0], [-0.2 1.0; -0.2 1.0]};

exper.subjects = {
%   'EBIRD049';
%   'EBIRD002';
%   'EBIRD003';
%   'EBIRD004';
%   'EBIRD005';
%   'EBIRD006';
%   'EBIRD007';
%   'EBIRD008';
%   'EBIRD009';
%   'EBIRD010';
%   'EBIRD011';
%   'EBIRD012';
%   'EBIRD013';
%   'EBIRD014';
%   'EBIRD015';
%   'EBIRD016';
%   'EBIRD017';
%   'EBIRD018';
%   'EBIRD019';
%   'EBIRD020';
  'EBIRD021';
  };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {{'session_1'}};
% exper.sessions = {{'session_2'}};
% exper.sessions = {{'session_1'}, {'session_8'}, {'session_9'}};
% exper.sessions = {...
%   {'session_1'}, ...
%   {'session_2'}, ...
%   {'session_3'}, ...
%   {'session_4'}, ...
%   {'session_5'}, ...
%   {'session_6'}, ...
%   {'session_7'}, ...
%   {'session_8'}, ...
%   {'session_9'}};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.subDir = '';
dirs.behDir = fullfile(exper.name,'Behavioral','Sessions',dirs.subDir);
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
ana.trialFxn = 'ebird_trialfun';
ana.allowTrialOverlap = false;
ana.renumberSamplesContiguous = false;
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
  %ana.sessionNames = {'pretest','train1','train2','train3','train4','train5','train6','posttest','posttest_delay'};
  ana.sessionNames = {'pretest'};
%   ana.sessionNames = {'pretest','posttest','posttest_delay'};
%   ana.sessionNames = {'pretest', 'train1'};
%   ana.sessionNames = {'train1'};
  %ana.sessionNames = {'train2'};
  
  % phases occur within a session; for dealing with events.mat
%   ana.phaseNames = {...
%     {'match'}, {'nametrain', 'name', 'name'}, {'name', 'name', 'name', 'name'}, ...
%     {'name', 'name', 'name', 'name'}, {'name', 'name', 'name', 'name'}, {'name', 'name', 'name', 'name'}, ...
%     {'name', 'name', 'name', 'name'}, {'match'}, {'match'}};
%   ana.phaseNames = {{'match'},{'nametrain', 'name', 'name'}};
  ana.phaseNames = {{'match'}};
%   ana.phaseNames = {{'match'}, {'match'}, {'match'}};
%   ana.phaseNames = {{'nametrain', 'name', 'name'}};
  %ana.phaseNames = {{'name', 'name', 'name', 'name'}};
  
  % types of event info to store in trialinfo field; must correspond to
  % values listed in exper.eventValues
  ana.trl_order.match_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'stimNum', 'imgCond', 'isSubord', 'trained', 'sameSpecies', 'response', 'rt', 'acc'};
  ana.trl_order.nametrain_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'block', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'imgCond', 'isSubord', 'response', 'rt', 'acc'};
  ana.trl_order.name_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'imgCond', 'isSubord', 'response', 'rt', 'acc'};
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

% artifact settings
% ana.artifact.type = {'ftManual', 'ftICA'};
ana.artifact.type = {'ftAuto'};
ana.artifact.reject = 'complete';
ana.artifact.resumeManArtFT = false;
ana.artifact.resumeICACompFT = false;
% % negative trlpadding: don't check that time (on both sides) for artifacts
% IMPORTANT: Not used for threshold artifacts. only use if segmenting a lot
% of extra time around trial epochs. Otherwise set to zero.
% ana.artifact.trlpadding = -0.5;
ana.artifact.trlpadding = 0;
ana.artifact.artpadding = 0.1;
ana.artifact.fltpadding = 0;
% ana.artifact.threshmin = -150;
% ana.artifact.threshmax = 150;
% ana.artifact.threshrange = 250;
ana.artifact.threshmin = -200;
ana.artifact.threshmax = 200;
ana.artifact.threshrange = 350;
ana.artifact.basic_art_z = 60;
% ana.artifact.muscle_art_z = 70;
ana.artifact.jump_art_z = 70;
ana.artifact.threshmin_postICA = -100;
ana.artifact.threshmax_postICA = 100;
ana.artifact.threshrange_postICA = 150;
ana.artifact.basic_art_z_postICA = 30;
% ana.artifact.muscle_art_z_postICA = 50;
ana.artifact.jump_art_z_postICA = 50;

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
cfg_pp.implicitref = 'Cz';
% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];
% single precision to save space
%cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.keeptrials = 'yes';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs_multiSes(exper,ana,cfg_proc,dirs,files,'tla',false);
% [dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct_multiSes(ana,cfg_pp,exper,dirs,files);
% [exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);

process_ft_data_multiSes(ana,cfg_proc,exper,dirs,files,cfg_pp);

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
%   [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails_multiSes(saveFile,additional_AD_file,backup_orig_AD,sortBySubj,replaceOrig);
% end

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
