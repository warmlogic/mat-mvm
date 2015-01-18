function space2_ftprocess_tla_wrapper(whichStages)
% space_ftprocess_tla_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub space_ftprocess_tla_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There are two stages:
%  stage1 = call wrapper that calls create_ft_struct_multiSes (which calls seg2ft),
%  saves the raw files
%  stage2 = calls ft_timelockanalysis and saves the processed files
%
% Input:
%  whichStages: the stage number(s) to run (default = 1:2)
%
% Output:
%  timelock data

% check/handle arguments
error(nargchk(0,1,nargin))
% STAGES = 1:2;
STAGES = 1;
if nargin == 1
  STAGES = whichStages;
end

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% exper.eventValues = {{'expo_stim'}};

% pre- and post-stimulus times to read, in seconds (pre is negative).
% Construct as a cell with one Nx2 matrix per session where N is
% length(exper.eventValues{ses}) Order must correspond to the event order
% in exper.eventValues.
exper.prepost = {[-1.0 2.0; -1.0 2.0; -1.0 2.0]};
% exper.prepost = {[-1.0 2.0]};

exper.subjects = {
%   'SPACE2001';
%   'SPACE2002';
%   %'SPACE2003'; % DNF session 2
%   'SPACE2004';
%   'SPACE2005';
%   'SPACE2006';
%   'SPACE2007'; % terrible performance
%   'SPACE2008';
%   %'SPACE2009'; % DNF session 2
%   'SPACE2010';
%   'SPACE2011';
%   'SPACE2012';
%   %'SPACE2013'; % didn't record EEG, stopped session 1 in middle
%   'SPACE2014';
%   'SPACE2015';
%   'SPACE2016';
%   'SPACE2017'; % terrible performance
%   'SPACE2018';
%   'SPACE2019';
%   % 'SPACE2020'; % DNF session 2
%   'SPACE2021';
%   'SPACE2022';
%   %'SPACE2023'; % DNF session 2
%   %'SPACE2024'; % DNF session 2
%   'SPACE2025';
%   'SPACE2026';
%   %'SPACE2027'; % DNF session 2
%   'SPACE2028';
%   'SPACE2029';
%   'SPACE2030';
%   'SPACE2031';
%   'SPACE2032';
%   'SPACE2033';
%   'SPACE2034';
%   'SPACE2035';
%   'SPACE2036';
%   'SPACE2037';
  'SPACE2038';
  %'SPACE2039'; % DNF session 2
%   'SPACE2040';
  };


% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, other directories (ns_evt, ns_bci, etc).
% They are not necessarily the session directory names where the FieldTrip
% data is saved for each subject because of the option to combine sessions.
% See 'help create_ft_struct' for more information.
% exper.sessions = {{'session_1'}};
exper.sessions = {{'session_1', 'session_2'}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POSSIBLY MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  runLocally = 1;
elseif isfield(dirs,'serverLocalDir') && exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
  runLocally = 1;
elseif isfield(dirs,'dreamDir') && exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
  runLocally = 0;
elseif isfield(dirs,'localDir') && exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
  runLocally = 1;
else
  error('Data directory not found.');
end

% Use the FT chan locs file
files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
exper.refChan = 'Cz';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

%% set up analysis parameters

% raw data
ana.segFxn = 'seg2ft';

% 17ms avg from photocell test, + 8ms from EGI's A/D correction (NS 4.4)
ana.offsetMS = 25;

ana.continuous = 'yes';
ana.trialFxn = 'space2_trialfun_mff';
ana.allowTrialOverlap = true;
ana.renumberSamplesContiguous = true;
% files used when adding metadata to segmented trials
ana.useMetadata = true;
ana.metadata.types = {'eventStruct','nsEvt'};
ana.useExpInfo = true;
% ana.evtToleranceMS = 8; % 2 samples @ 250 Hz
ana.usePhotodiodeDIN = false;
ana.photodiodeDIN_toleranceMS = 20;
ana.photodiodeDIN_str = 'DIN ';
if ana.useExpInfo
  % possible sessions and phases
  ana.sessionNames = {'day1','day2'};
  
  % phases occur within a session; for dealing with events.mat
  ana.phaseNames = {{'multistudy', 'cued_recall_only'}};
  
  % types of event info to store in trialinfo field; must correspond to
  % values listed in exper.eventValues
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

% this will turn off ana.cfg_cont and ana.artifact.continuousXXX so no
% processing is done; will only save FieldTrip events
ana.returnAfterSavingFtEvents = false;

% artifact settings
ana.artifact.reject = 'complete';
ana.artifact.preArtBaseline = 'yes';

ana.artifact.type = {'nsClassic','ftAuto'};
% ana.artifact.type = {'nsClassic'};

% set up for nsClassic
ana.artifact.checkArtSec = [0 1.0];
ana.artifact.blink_threshold = 70;
ana.artifact.fast_threshold = 100;
ana.artifact.diff_threshold = 50;
ana.artifact.rejectTrial_nBadChan = 10;
ana.artifact.repairChan_percentBadTrials = 20;
ana.artifact.allowBadNeighborChan = false;

% set up for ftAuto following nsClassic
% negative trlpadding: don't check that time (on both sides) for artifacts.
% IMPORTANT: Not used for threshold artifacts. only use if segmenting a lot
% of extra time around trial epochs. Otherwise set to zero.
ana.artifact.trlpadding = -1.0;
% ana.artifact.trlpadding = 0;
ana.artifact.artpadding = 0.1;
ana.artifact.fltpadding = 0;

ana.artifact.thresh = true;
ana.artifact.threshmin = -150;
ana.artifact.threshmax = 150;
ana.artifact.threshrange = 250;
ana.artifact.basic_art = true;
ana.artifact.basic_art_z = 30;
ana.artifact.jump_art = true;
ana.artifact.jump_art_z = 50;
% eog_art is only used with ftAuto
ana.artifact.eog_art = false;
% ana.artifact.eog_art_z = 3.5;

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
% % eog_art is only used with ftAuto
% ana.artifact.eog_art = false;
% ana.artifact.eog_art_z = 3.5;
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

%% set up processing parameters

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
% [dirs,files] = mm_ft_setSaveDirs_multiSes(exper,ana,cfg_proc,dirs,files,'tla',false);
% [dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

%% set up for running stages and specifics for Dream

% name(s) of the functions for different stages of processing
stageFun = {@stage1};
% timeOut  = {2}; % in HOURS
timeOut  = {18}; % in HOURS
% stageFun = {@stage1,@stage2};
% timeOut  = {2,2}; % in HOURS

if runLocally == 0
  % need to export DISPLAY to an offscreen buffer for MATLAB DCS graphics
  if exist('parcluster','file')
    sched = parcluster();
  else
    sched = findResource();
  end
  if strcmp(sched.Type, 'generic')
    setenv('DISPLAY', 'dream:99');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%capture diary and time statistics
thisRun = [exper.name,'_overview_',datestr(now,'ddmmmyyyy-HHMMSS')];
diary(fullfile(dirs.saveDirProc,[thisRun '.log']));
tStart = tic;
fprintf('START TIME: %s\n',datestr(now,13));
for i = STAGES
  tS = tic;
  fprintf('STAGE%d START TIME: %s\n',i, datestr(now,13));
  
  % execute the processing stages
  if i == 1
    stageFun{i}(ana,cfg_pp,cfg_proc,exper,dirs,files,runLocally,timeOut{i});
  %elseif i == 2
  %  stageFun{i}(ana,cfg_pp,cfg_proc,exper,dirs,files,runLocally,timeOut{i});
  end
  
  fprintf('STAGE%d END TIME: %s\n',i, datestr(now,13));
  fprintf('%.3f -- elapsed time STAGE%d (seconds)\n', toc(tS), i);
end
time = toc(tStart);
fprintf('%.3f -- elapsed time OVERALL (seconds)\n', time);
fprintf('END TIME: %s\n',datestr(now,13));
diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stage1(ana,cfg_pp,cfg_proc,exper,dirs,files,runLocally,timeOut)
% stage1: process the input files with FieldTrip based on the analysis
% parameters

%% Process the data
if runLocally == 0
  %% Dream: create one task for each subject (i.e., submit one to each node)
  
  % start a new job
  job = newJob(dirs);
  
  % save the original subjects array so we can set exper to have single
  % subjects, one for each task created
  allSubjects = exper.subjects;
  
  for i = 1:length(allSubjects)
    fprintf('Processing %s...\n',allSubjects{i});
    
    % Dream: create one task for each subject
    exper.subjects = allSubjects(i);
    
    % inArg = {ana,cfg_pp,exper,dirs,files};
    % save the exper struct (output 1) so we can use it later
    %createTask(job,@create_ft_struct_multiSes,1,inArg);
    
    inArg = {ana,cfg_pp,exper,dirs,files,cfg_proc};
    createTask(job,@create_ft_struct_multiSes_cluster,1,inArg);
    
  end
  
  runJob(job,timeOut,fullfile(dirs.saveDirProc,[exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
  
  % % exper = getAllOutputArguments(job);
  %tasks = get(job,'Tasks');
  %exper = tasks(i).OutputArguments{1};
  %process_ft_data_multiSes(ana,cfg_proc,exper,dirs,files,cfg_pp);
  
  % get the trial counts together across subjects, sessions, and events
  %[exper] = mm_ft_concatTrialCounts_cluster(job,exper,allSubjects);

%   % save the analysis details; overwrite if it already exists
%   saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails.mat'));
%   %if ~exist(saveFile,'file')
%   fprintf('Saving %s...',saveFile);
%   save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
%   fprintf('Done.\n');
%   %else
%   %  error('Not saving! %s already exists.\n',saveFile);
%   %end
  
  % final step: destroy the job because this doesn't happen in runJob
  destroy(job);
  
else
  %% run the function locally
  
  % create a log of the command window output
  thisRun = [exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS')];
  % turn the diary on
  diary(fullfile(dirs.saveDirProc,[thisRun,'.log']));
  
  % use the peer toolbox
  %ana.usePeer = 1;
  ana.usePeer = 0;
  
  % Local: run all the subjects
  [exper] = create_ft_struct_multiSes(ana,cfg_pp,exper,dirs,files);
  
  process_ft_data_multiSes(ana,cfg_proc,exper,dirs,files,cfg_pp);
  
%   % save the analysis details; overwrite if it already exists
%   saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails.mat'));
%   %if ~exist(saveFile,'file')
%   fprintf('Saving %s...',saveFile);
%   save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
%   fprintf('Done.\n');
%   %else
%   %  error('Not saving! %s already exists.\n',saveFile);
%   %end
  
  % turn the diary off
  diary off
end

% function stage2(ana,cfg_pp,cfg_proc,exper,dirs,files,runLocally,timeOut)
% % stage2: process the input files with FieldTrip based on the analysis
% % parameters
% 
% %% Process the data
% if runLocally == 0
%   %% Dream: create one task for each subject (i.e., submit one to each node)
%   
%   % start a new job
%   job = newJob(dirs);
%   
%   %adFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
%   %[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
%   
%   %sdFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr,'subjectDetails.mat');
%   %load(sdFile,'exper','ana','dirs','files','cfg_pp');
%   
%   % save the original subjects array so we can set exper to have single
%   % subjects, one for each task created
%   allSubjects = exper.subjects;
%   
%   for i = 1:length(allSubjects)
%     fprintf('Processing %s...\n',allSubjects{i});
%     
%     % Dream: create one task for each subject
%     exper.subjects = allSubjects(i);
%     
%     %inArg = {ana,cfg_proc,exper,dirs};
%     inArg = {ana,cfg_proc,exper,dirs,files,cfg_pp};
%     
%     % save the exper struct (output 1) so we can use it later
%     createTask(job,@process_ft_data_multiSes,0,inArg);
%   end
%   
%   runJob(job,timeOut,fullfile(dirs.saveDirProc,[exper.name,'_stage2_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
%   
%   % final step: destroy the job because this doesn't happen in runJob
%   destroy(job);
%   
% else
%   %% run the function locally
%   
%   % create a log of the command window output
%   thisRun = [exper.name,'_stage2_',datestr(now,'ddmmmyyyy-HHMMSS')];
%   % turn the diary on
%   diary(fullfile(dirs.saveDirProc,[thisRun,'.log']));
%   
%   % use the peer toolbox
%   %ana.usePeer = 1;
%   ana.usePeer = 0;
%   
%   % Local: run all the subjects
%   %process_ft_data(ana,cfg_proc,exper,dirs);
%   process_ft_data_multiSes(ana,cfg_proc,exper,dirs,files,cfg_pp);
%   
%   % turn the diary off
%   diary off
% end

function job = newJob(dirs)
% newJob Creates a new PCT job and sets job's dependencies
%
%   dirs -- data structure with necessary fields like data locations

% Set up scheduler, job
if exist('parcluster','file')
  sched = parcluster();
else
  sched = findResource();
end
job = createJob(sched);
if verLessThan('matlab','8.3')
  % define the directories to add to worker sessions' matlab path
  homeDir = getenv('HOME');
  myMatlabDir = fullfile(homeDir,'Documents','MATLAB');
  p = path();
  set(job, 'PathDependencies', {homeDir, myMatlabDir, pwd(), p, dirs.dataroot});
end

function runJob( job, timeOut, logFile )
% runJob Submits and waits on job to finish or timeout
%  runJob will submit the supplied job to the scheduler and will
% wait for the job to finish or until the timeout has been reached. 
% If the job finishes, then the command window outputs of all tasks
% are appended to the log file and the job is destroyed.
%   If the timeout is reached, an error is reported but the job is not
% destroyed.
%
%   job -- the job object to submit
%   timeOut -- the timeout value in hours
%   logFile -- full file name of the log file to append output to
%
% Example:
%       runJob( job, 5, 'thisrun.log');

% check/handle arguments
error(nargchk(1,3,nargin))
TIMEOUT=3600*5; % default to 5 hours 
if nargin > 1
  TIMEOUT=timeOut*3600;
end
LOGFILE=[job.Name '.log'];
if nargin > 2
  LOGFILE = logFile;
end

% Capture command window output from all tasks
alltasks = get(job, 'Tasks');
if verLessThan('matlab','8.3')
  set(alltasks, 'CaptureCommandWindowOutput', true);
end

% Submit Job/Tasks and wait for completion (or timeout)
submit(job)
finished = waitForState(job, 'finished', TIMEOUT);
if finished
  errors = logOutput(alltasks, LOGFILE);
  if errors
    error([mfilename ':logOutput'],'%s had %d errors',job.Name, errors)
  %elseif ~errors
  %  destroy(job);
  end
else
  error([mfilename ':JobTimeout'],'%s: Timed out waiting for job...NAME: %s',...
    datestr(now, 13), job.Name, job.ID, job.StartTime)
end

function numErrors=logOutput( tasks, logFile )
% logOutput - concatenates tasks' output into a logfile
%   tasks -- the tasks to capture output from 
%   logFile -- the file to log the output to
%   numErrors -- number of tasks which failed

% check for argument(s)
error(nargchk(2,2,nargin))

numErrors=0;
try
  fid=fopen(logFile, 'a+');
  for i=1:length(tasks)
    fprintf(fid,'\n***** START TASK %d *****\n',i);
    fprintf(fid,'%s\n', tasks(i).CommandWindowOutput);
    if ~isempty(tasks(i).Error.stack)
      numErrors = numErrors +1;
      % write to log file
      fprintf( fid, 'ERROR: %s\n', tasks(i).Error.message );
      fprintf( fid, '%s\n', tasks(i).Error.getReport );
      % write to standard error
      fprintf( 2, 'ERROR: %s\n', tasks(i).Error.message );
      fprintf( 2, '%s\n', tasks(i).Error.getReport );
    end
    fprintf(fid,'\n***** END TASK %d *****\n',i);
  end
  fclose(fid);
catch ME
  disp(ME)
  warning([mfilename ':FailOpenLogFile'],...
    'Unable to write log file with task output...');
end
