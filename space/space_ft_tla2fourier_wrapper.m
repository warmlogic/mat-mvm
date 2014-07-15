function space_ft_tla2fourier_wrapper(whichStages)
% space_ft_tla2fourier_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub space_ft_tla2fourier_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There is one stage:
%  stage1 = call wrapper that calls mm_tla2fourier, saves fourier and pow
%
% Input:
%  whichStages: the stage number(s) to run (default = 1)
%
% Output:
%  time-frequency data

% check/handle arguments
error(nargchk(0,1,nargin))
STAGES = 1;
if nargin == 1
  STAGES = whichStages;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Experiment-specific setup

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
  runLocally = 1;
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
  runLocally = 1;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
  runLocally = 0;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
  runLocally = 1;
else
  error('Data directory not found.');
end

% procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';
procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');

subjects = {
%   'SPACE001';
%   'SPACE002';
%   'SPACE003';
%   'SPACE004';
%   'SPACE005';
%   'SPACE006';
%   'SPACE007';
%   %'SPACE008'; % didn't perform task correctly, didn't perform well
%   'SPACE009';
%   'SPACE010';
%   'SPACE011';
%   'SPACE012';
%   'SPACE013';
%   'SPACE014';
%   'SPACE015';
%   'SPACE016';
%   'SPACE017'; % old assessment: really noisy EEG, half of ICA components rejected
%   'SPACE018';
%   'SPACE019';
%   'SPACE020';
%   'SPACE021';
%   'SPACE022';
%   'SPACE027';
%   'SPACE029';
%   'SPACE037';
%   %'SPACE039'; % noisy EEG; original EEG analyses stopped here
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

% allowRecallSynonyms = true;

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POSSIBLY MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set up analysis parameters

cfg_ana = [];
cfg_ana.orig_ftype = 'tla';
cfg_ana.orig_param = 'trial';

% must match with cfg_ft.output, below
cfg_ana.out_param = 'fourierspctrm';

% upon loading data, remove data.cfg.previous field to save memory
cfg_ana.rmPreviousCfg = true;

cfg_ana.output2alt = true;
cfg_ana.alt_ftype = 'pow';
cfg_ana.alt_param = 'powspctrm';
cfg_ana.useLockFiles = false;
cfg_ana.splitTrials = true;
if cfg_ana.splitTrials
  cfg_ana.splitSize = 100;
  % number of trials to lump in with the last trial split. Must be >= 1.
  cfg_ana.splitRemainderLump = 10;
  
  % whether to see if files are too big to combine (RAM limit issue)
  cfg_ana.checkSplitFileSizeForSaving = false;
  if cfg_ana.checkSplitFileSizeForSaving
    % approximate limit in MB that the combined files can reach; otherwise
    % will keep the split files
    cfg_ana.combineSavingLimitMB = 3000;
  end
end

cfg_ft = [];
cfg_ft.pad = 'maxperlen';
% cfg_ft.output = 'pow';
% % cfg_ft.output = 'powandcsd';
% cfg_ft.keeptrials = 'yes';
% cfg_ft.keeptapers = 'no';

% cfg_ft.output = 'pow';
% cfg_ft.keeptrials = 'no';
cfg_ft.output = 'fourier';
cfg_ft.keeptrials = 'yes';
% cfg_ft.keeptapers = 'yes';

% wavelet
cfg_ft.method = 'wavelet';
cfg_ft.width = 5;
%cfg_ft.toi = -0.8:0.04:3.0;
cfg_ft.toi = -0.5:0.04:1.0;
% cfg_ft.toi = 0:0.04:1.0;
% % evenly spaced frequencies, but not as many as foilim makes
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% % cfg_ft.foi = 3:freqstep:9;
% cfg_ft.foi = 3:freqstep:60;
% cfg_ft.foi = 4:1:100;
%cfg_ft.foi = 4:1:30;
cfg_ft.foi = 3:1:80;
% cfg_ft.foilim = [3 9];
% cfg_ft.foilim = [3 13];
% cfg_ft.foilim = [3 80];

cfg_alt = [];
cfg_alt.variance = 'no';
cfg_alt.jackknife = 'no';
cfg_alt.keeptrials = 'yes';

cfg_ana.saveroot = strrep(dirs.saveDirProc,cfg_ana.orig_ftype,cfg_ft.output);
if ~exist(cfg_ana.saveroot,'dir')
  [s] = mkdir(cfg_ana.saveroot);
  if ~s
    error('Could not make %s',cfg_ana.saveroot);
  end
end

cfg_ana.saveroot_alt = strrep(dirs.saveDirProc,cfg_ana.orig_ftype,cfg_ana.alt_ftype);
if ~exist(cfg_ana.saveroot_alt,'dir')
  [s] = mkdir(cfg_ana.saveroot_alt);
  if ~s
    error('Could not make %s',cfg_ana.saveroot_alt);
  end
end

%% set up for running stages and specifics for Dream

% name(s) of the functions for different stages of processing
stageFun = {@stage1};
timeOut  = {4}; % in HOURS
% stageFun = {@stage1,@stage2};
% timeOut  = {2,2}; % in HOURS

if runLocally == 0
  % need to export DISPLAY to an offscreen buffer for MATLAB DCS graphics
  sched = findResource();
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
    stageFun{i}(exper,dirs,cfg_ana,cfg_ft,cfg_alt,runLocally,timeOut{i});
  %elseif i == 2
  %  stageFun{i}(ana,cfg_proc,exper,dirs,runLocally,timeOut{i});
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

function stage1(exper,dirs,cfg_ana,cfg_ft,cfg_alt,runLocally,timeOut)
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
    
    inArg = {exper,dirs,cfg_ana,cfg_ft,cfg_alt};
    
    % save the exper struct (output 1) so we can use it later
    nOutArg = 0;
    createTask(job,@mm_tla2fourier,nOutArg,inArg);
  end
  
  runJob(job,timeOut,fullfile(dirs.saveDirProc,[exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
  
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
  %ana.usePeer = 0;
  
  % Local: run all the subjects
  mm_tla2fourier(exper,dirs,cfg_ana,cfg_ft,cfg_alt);
  
  % turn the diary off
  diary off
end

function job = newJob(dirs)
% newJob Creates a new PCT job and sets job's dependencies
%
%   dirs -- data structure with necessary fields like data locations

% Set up scheduler, job
sched = findResource();
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
