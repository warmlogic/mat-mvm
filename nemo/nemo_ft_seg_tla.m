% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'FAV';

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files; or space_trialfun.m must be set up to find the
% corrct events
exper.eventValues = {{'ao_standard_corr','ao_target_corr', 'study_prime', 'study_targ_AA_rel',...
    'study_targ_AA_unrel','study_targ_CA_rel','study_targ_CA_unrel','test_targ_A_new_corr',...
    'test_targ_AA_unrel_corr','test_targ_CA_unrel_corr'}};

% pre- and post-stimulus times to read, in seconds (pre is negative);
% because they get sorted, must correspond to the order listed in
% exper.eventValues
exper.prepost = {[...
  -0.2 0.65; ...
  -0.2 0.65; ...
  -0.2 1.0; ...
  -0.2 1.0; ...
  -0.2 1.0; ...
  -0.2 1.0; ...
  -0.2 1.0; ...
  -0.2 1.0; ...
  -0.2 1.0; ...
  -0.2 1.0]};

exper.subjects = {
%  'TC_NEMO-001';
%  'TC_NEMO-002';
%  'TC_NEMO-003';
%  'TC_NEMO-004';
%  'TC_NEMO-005';
%  'TC_NEMO-006';
%  'TC_NEMO-007';
%  'TC_NEMO-008';
%  'TC_NEMO-009';
%  'TC_NEMO-010';
  'TC_NEMO-011';
  'TC_NEMO-012';
%  'TC_NEMO-013';
%  'TC_NEMO-014';
%  'TC_NEMO-015';
%  'TC_NEMO-016';
%  'TC_NEMO-017';
%  'TC_NEMO-018';
%  'TC_NEMO-019';
%  'TC_NEMO-020';
%  'TC_NEMO-021';
%  'TC_NEMO-022';
%  'TC_NEMO-023';
%  'TC_NEMO-024';
%  'TC_NEMO-025';
%  'TC_NEMO-026';
%  'TC_NEMO-027';
%  'TC_NEMO-028';
%  'TC_NEMO-029';
%  'TC_NEMO-030';
%  'TC_NEMO-031';
%  'TC_NEMO-032';
%  'TC_NEMO-033';
%  'TC_NEMO-034';
%  'TC_NEMO-035';
%  'TC_NEMO-036';
%  'TC_NEMO-037';
%  'TC_NEMO-038';
%  'TC_NEMO-039';
%  'TC_NEMO-040';
%  'TC_NEMO-041';
%  'TC_NEMO-042';
%  'TC_NEMO-043';
%  'TC_NEMO-044';
%  'TC_NEMO-045';
%  'TC_NEMO-046';
%  'TC_NEMO-047';
%  'TC_NEMO-048';
%  'TC_NEMO-049';
%  'TC_NEMO-050';
%  'TC_NEMO-051';
%  'TC_NEMO-052';
%  'TC_NEMO-053';
%  'TC_NEMO-054';
%  'TC_NEMO-055';
  };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {{'session_1'}};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.subDir = '';
% dirs.dataDir = fullfile(exper.name,'EEG','Sessions','face_house_ratings','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
dirs.behDir = fullfile(exper.name,'NEMO','Behavioral','Sessions',dirs.subDir);
% dirs.dataDir = fullfile(exper.name,'EEG','Sessions','ftpp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
dirs.dataDir = fullfile(exper.name,'NEMO','EEG','Sessions','ftpp',dirs.subDir);
% Possible locations of the data files (dataroot)
%dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
%dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
%dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'Desktop','markupfiles');


% pick the right dirs.dataroot
if exist('/Volumes/CurranImacLab3','dir')
    dirs.dataroot = dirs.localDir;
else
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
end

% Use the FT chan locs file
files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';

ana.continuous = 'yes';
ana.trialFxn = 'nemo_trialfun';
ana.allowTrialOverlap = true;
ana.renumberSamplesContiguous = true;
% files used when adding metadata to segmented trials
ana.useMetadata = true;
ana.metadata.types = {'nsEvt'};
ana.useExpInfo = true;
ana.evtToleranceMS = 8; % 2 samples @ 250 Hz
ana.usePhotodiodeDIN = false;
ana.photodiodeDIN_toleranceMS = 40;
ana.photodiodeDIN_str = 'DIN ';
if ana.useExpInfo
  % possible sessions and phases
   ana.sessionNames = {'NEMO'};
  
  % phases occur within a session; for dealing with events.mat
  ana.phaseNames = {{'TC_NEMO_fN400study','TC_NEMO_AO','TC_NEMO_fN400test'}};
  %ana.phaseNames = {{'expo'}};
  %ana.phaseNames = {{'multistudy'}};
  %ana.phaseNames = {{'distract_math'}};
  %ana.phaseNames = {{'cued_recall'}};
  
  % types of event info to store in trialinfo field; must correspond to
  % values listed in exper.eventValues
   ana.trl_order.ao_standard_corr = {'eventNumber', 'phaseType', 'trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number'};
 
  ana.trl_order.ao_target_corr = {'eventNumber', 'phaseType', 'trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number'}; 
  
  ana.trl_order.study_prime = {'eventNumber', 'phaseType', 'trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
  
  ana.trl_order.study_targ_AA_rel = {'eventNumber', 'phaseType', 'trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
 
  ana.trl_order.study_targ_AA_unrel = {'eventNumber', 'phaseType', 'trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
  
  ana.trl_order.study_targ_CA_rel = {'eventNumber', 'phaseType', 'trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
  
  ana.trl_order.study_targ_CA_unrel = {'eventNumber', 'phaseType', 'trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
  
  ana.trl_order.test_targ_A_new_corr = {'eventNumber', 'phaseType','trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
  
  ana.trl_order.test_targ_AA_unrel_corr = {'eventNumber', 'phaseType','trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
  
  ana.trl_order.test_targ_CA_unrel_corr = {'eventNumber', 'phaseType','trial','cell_label',...
      'times_studied','resp_value','accuracy','reaction_time','random_number',...
      'FSG_target','BSG_target','CNC_target','KFFRQ_target','NLET_target','NPHON_target',...
      'NSYLL_target','orthoN_target','phonoN_target','AOA_target','BFRQ_target',...
      'FAM_target','IMG_target','CMEAN_target','PMEAN_target','TLFRQ_target',...
      'wordid_target','CNC_prime','KFFRQ_prime','NLET_prime','NPHON_prime',...
      'NSYLL_prime','orthoN_prime','phonoN_prime','AOA_prime','BFRQ_prime','FAM_prime',...
      'IMG_prime','CMEAN_prime','PMEAN_prime','TLFRQ_prime','wordid_prime','word_pair_id'};
  
  
  
%   ana.trl_order.multistudy_word = ana.trl_order.multistudy_image;
  %ana.trl_order.distract_math_stim = {'eventNumber', 'sesType', 'phaseType', 'response', 'acc', 'rt'};
%  ana.trl_order.test_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'pairNum', 'recog_resp', 'recog_acc', 'recog_rt', 'new_resp', 'new_acc', 'new_rt', 'recall_resp', 'recall_spellCorr', 'recall_rt'};
end

% preprocess continuous data in these ways
ana.cfg_cont.lpfilter = 'yes';
ana.cfg_cont.lpfreq = 40; %originally 100, NS tool set to 40Hz lowpass
ana.cfg_cont.hpfilter = 'yes';
ana.cfg_cont.hpfreq = 0.1;
ana.cfg_cont.hpfilttype = 'but';
ana.cfg_cont.hpfiltord = 4;
ana.cfg_cont.bsfilter = 'no'; %originally 'yes', but no bandpass filter used in NS processing
ana.cfg_cont.bsfreq = [59 61];

% artifact settings
ana.artifact.type = {'ftAuto'}; %, {'ftManual', 'ftICA'};
ana.artifact.reject = 'complete';
ana.artifact.resumeManArtFT = false;
ana.artifact.resumeICACompFT = false;
% negative trlpadding: don't check that time (on both sides) for artifacts
% ana.artifact.trlpadding = -0.5;
ana.artifact.trlpadding = 0;
ana.artifact.artpadding = 0.1;
ana.artifact.fltpadding = 0;
ana.artifact.threshmin = -70; %originally -150; - NS art det classic uses 70
ana.artifact.threshmax = 70; %originally 150; - NS art det classic uses 70
ana.artifact.threshrange = 250;
ana.artifact.basic_art_z = 40;
% ana.artifact.muscle_art_z = 60;
ana.artifact.jump_art_z = 60;
ana.artifact.threshmin_postICA = -100;
ana.artifact.threshmax_postICA = 100;
ana.artifact.threshrange_postICA = 150;
ana.artifact.basic_art_z_postICA = 30;
% ana.artifact.muscle_art_z_postICA = 50;
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
cfg_pp.demean = 'no'; %originally 'yes', but no baseline corr used in NS processing
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

% % %% get the bad channel information
% % 
% % cfg = [];
% % cfg.badChanManual = false;
% % cfg.badChanEP = true;
% % [exper] = mm_getBadChan(cfg,exper,dirs);
% 
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
% emailme = 1;
% if emailme
%   subject = sprintf('Done with%s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
%   mail_message = {...
%     sprintf('Done with%s %s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
%     sprintf('%s',saveFile),...
%     };
%   send_gmail(subject,mail_message);
% end
