function ebird_ftprocess_tla_wrapper_2stage(whichStages)
% ebird_ftprocess_tla_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub ebird_ftprocess_tla_wrapper.m
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
STAGES = 1:2;
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

exper.name = 'EBIRD';

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files; or space_trialfun.m must be set up to find the
% corrct events
% exper.eventValues = {{'match_stim'}};
% exper.eventValues = {{'match_stim'}, {'match_stim'}, {'match_stim'}};
% exper.eventValues = {{'match_stim'}, {'nametrain_stim', 'name_stim'}};
exper.eventValues = {{'match_stim'}, {'nametrain_stim', 'name_stim'}, ...
  {'name_stim'}, {'name_stim'}, {'name_stim'}, {'name_stim'}, {'name_stim'}, ...
  {'match_stim'}, {'match_stim'}};

% pre- and post-stimulus times to read, in seconds (pre is negative).
% Construct as a cell with one Nx2 matrix per session where N is
% length(exper.eventValues{ses}) Order must correspond to the event order
% in exper.eventValues.
% exper.prepost = {[-0.2 1.0]};
% exper.prepost = {[-0.2 1.0], [-0.2 1.0], [-0.2 1.0]};
% exper.prepost = {[-0.2 1.0], [-0.2 1.0; -0.2 1.0]};
exper.prepost = {[-0.2 1.0], [-0.2 1.0; -0.2 1.0], [-0.2 1.0], [-0.2 1.0], [-0.2 1.0], [-0.2 1.0], [-0.2 1.0], [-0.2 1.0], [-0.2 1.0]};

exper.subjects = {
  %'EBIRD049'; % Pilot. (due to short ses1 match, missing ses2 name)
  %'EBIRD002'; % Pilot. (due to short ses1 match, missing ses2 name)
  %'EBIRD003'; % Pilot. (due to missing ses7 name) - NB: LAST PILOT TO BE REPLACED
  %'EBIRD004'; % DNF. Dropout. Last session: 8.
  'EBIRD005';
  %'EBIRD006'; % DNF. Dropout. Last session: 2.
  'EBIRD007';
  'EBIRD008';
  'EBIRD009';
  'EBIRD010';
  'EBIRD011';
  'EBIRD012';
  %'EBIRD013'; % DNF. Dropout. Last session: 5. Lost session 6 in HD crash.
  %'EBIRD014'; % DNF. Rejected. Last session: 1.
  %'EBIRD015'; % DNF. Lost in HD crash.
  %'EBIRD016'; % DNF. Lost in HD crash.
  %'EBIRD017'; % DNF. Lost in HD crash.
  'EBIRD018';
  'EBIRD019';
  'EBIRD020';
  'EBIRD021';
  %'EBIRD022'; % DNF. Dropout. Last session: 8.
  %'EBIRD023'; % DNF. Dropout. Last session: 1.
  'EBIRD024';
  'EBIRD025';
  'EBIRD027';
  'EBIRD029';
  'EBIRD032';
  'EBIRD034';
  'EBIRD042';
  };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
% exper.sessions = {{'session_1'}};
% exper.sessions = {{'session_2'}};
% exper.sessions = {{'session_1'}, {'session_8'}, {'session_9'}};
exper.sessions = {...
  {'session_1'}, ...
  {'session_2'}, ...
  {'session_3'}, ...
  {'session_4'}, ...
  {'session_5'}, ...
  {'session_6'}, ...
  {'session_7'}, ...
  {'session_8'}, ...
  {'session_9'}};

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

ana.continuous = 'yes';
ana.trialFxn = 'ebird_trialfun';
ana.allowTrialOverlap = false;
ana.renumberSamplesContiguous = false;
% files used when adding metadata to segmented trials
ana.useMetadata = true;
ana.metadata.types = {'eventStruct','nsEvt'};
ana.useExpInfo = true;
ana.usePhotodiodeDIN = true;
ana.photodiodeDIN_toleranceMS = 20;
ana.photodiodeDIN_str = 'DIN ';
if ana.useExpInfo
  % possible sessions and phases
  ana.sessionNames = {'pretest','train1','train2','train3','train4','train5','train6','posttest','posttest_delay'};
%   ana.sessionNames = {'pretest'};
%   ana.sessionNames = {'pretest','posttest','posttest_delay'};
%   ana.sessionNames = {'pretest', 'train1'};
%   ana.sessionNames = {'train1'};
  %ana.sessionNames = {'train2'};
  
  % phases occur within a session; for dealing with events.mat
  ana.phaseNames = {...
    {'match'}, {'nametrain', 'name', 'name'}, {'name', 'name', 'name', 'name'}, ...
    {'name', 'name', 'name', 'name'}, {'name', 'name', 'name', 'name'}, {'name', 'name', 'name', 'name'}, ...
    {'name', 'name', 'name', 'name'}, {'match'}, {'match'}};
%   ana.phaseNames = {{'match'},{'nametrain', 'name', 'name'}};
%   ana.phaseNames = {{'match'}};
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
ana.artifact.reject = 'complete';
ana.artifact.preArtBaseline = 'yes';

ana.artifact.type = {'nsClassic','ftAuto'};

% set up for nsClassic
ana.artifact.checkArtSec = [-Inf Inf];
ana.artifact.fast_threshold = 100;
ana.artifact.diff_threshold = 50;
ana.artifact.rejectTrial_nBadChan = 10;
ana.artifact.repairChan_percentBadTrials = 20;
ana.artifact.allowBadNeighborChan = false;

% set up for ftAuto following nsClassic
% % negative trlpadding: don't check that time (on both sides) for artifacts
% IMPORTANT: Not used for threshold artifacts. only use if segmenting a lot
% of extra time around trial epochs. Otherwise set to zero.
% ana.artifact.trlpadding = -0.5;
ana.artifact.trlpadding = 0;
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
[dirs,files] = mm_ft_setSaveDirs_multiSes(exper,ana,cfg_proc,dirs,files,'tla',false);
% [dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

%% set up for running stages and specifics for Dream

% name(s) of the functions for different stages of processing
%stageFun = {@stage1};
%timeOut  = {2}; % in HOURS
stageFun = {@stage1,@stage2};
timeOut  = {2,2}; % in HOURS

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
    [exper] = stageFun{i}(ana,cfg_pp,exper,dirs,files,runLocally,timeOut{i});
  elseif i == 2
    stageFun{i}(ana,cfg_pp,cfg_proc,exper,dirs,files,runLocally,timeOut{i});
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

function [exper] = stage1(ana,cfg_pp,exper,dirs,files,runLocally,timeOut)
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
    
    inArg = {ana,cfg_pp,exper,dirs,files};
    
    % save the exper struct (output 1) so we can use it later
    createTask(job,@create_ft_struct_multiSes,1,inArg);
  end
  
  runJob(job,timeOut,fullfile(dirs.saveDirProc,[exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
  
  %exper = getAllOutputArguments(job);
  
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

function stage2(ana,cfg_pp,cfg_proc,exper,dirs,files,runLocally,timeOut)
% stage2: process the input files with FieldTrip based on the analysis
% parameters

%% Process the data
if runLocally == 0
  %% Dream: create one task for each subject (i.e., submit one to each node)
  
  % start a new job
  job = newJob(dirs);
  
  %adFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
  %[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
  
  %sdFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr,'subjectDetails.mat');
  %load(sdFile,'exper','ana','dirs','files','cfg_pp');
  
  % save the original subjects array so we can set exper to have single
  % subjects, one for each task created
  allSubjects = exper.subjects;
  
  for i = 1:length(allSubjects)
    fprintf('Processing %s...\n',allSubjects{i});
    
    % Dream: create one task for each subject
    exper.subjects = allSubjects(i);
    
    %inArg = {ana,cfg_proc,exper,dirs};
    inArg = {ana,cfg_proc,exper,dirs,files,cfg_pp};
    
    % save the exper struct (output 1) so we can use it later
    createTask(job,@process_ft_data_multiSes,0,inArg);
  end
  
  runJob(job,timeOut,fullfile(dirs.saveDirProc,[exper.name,'_stage2_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
  
  % final step: destroy the job because this doesn't happen in runJob
  destroy(job);
  
else
  %% run the function locally
  
  % create a log of the command window output
  thisRun = [exper.name,'_stage2_',datestr(now,'ddmmmyyyy-HHMMSS')];
  % turn the diary on
  diary(fullfile(dirs.saveDirProc,[thisRun,'.log']));
  
  % use the peer toolbox
  %ana.usePeer = 1;
  ana.usePeer = 0;
  
  % Local: run all the subjects
  %process_ft_data(ana,cfg_proc,exper,dirs);
  process_ft_data_multiSes(ana,cfg_proc,exper,dirs,files,cfg_pp);
  
  % turn the diary off
  diary off
end

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
% define the directories to add to worker sessions' matlab path
homeDir = getenv('HOME');
myMatlabDir = fullfile(homeDir,'Documents','MATLAB');
p = path();
set(job, 'PathDependencies', {homeDir, myMatlabDir, pwd(), p, dirs.dataroot});

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
set(alltasks, 'CaptureCommandWindowOutput', true);

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
