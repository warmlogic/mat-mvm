function space_pow_clusterStat_wrapper(whichStages)
% space_pow_clusterStat_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub2011 space_pow_clusterStat_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There is one stage:
%  stage1 =
%
% Input:
%  whichStages: the stage number(s) to run (default = 1)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expName = 'SPACE';

subDir = '';
dataDir = fullfile(expName,'EEG','Sessions','ftpp',subDir);
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

% loads in pre-saved subject average file from this directory

% saveDirProc = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow');
saveDirProc = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/pow');

% loadFile = 'space_word_img_data_pow.mat';
loadFile = 'space_word_img_data_pow_wordBL.mat';

% subjects = {
%   %'SPACE001'; % low trial counts
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
%   %'SPACE017'; % really noisy EEG, half of ICA components rejected
%   'SPACE018';
%   %'SPACE019'; % low trial counts
%   'SPACE020';
%   'SPACE021';
%   'SPACE022';
%   'SPACE027';
%   'SPACE029';
%   'SPACE037';
%   'SPACE039'; % noisy EEG; original EEG analyses stopped here
%   'SPACE023';
%   'SPACE024';
%   'SPACE025';
%   'SPACE026';
%   'SPACE028';
%   %'SPACE030'; % low trial counts
%   'SPACE032';
%   'SPACE034';
%   'SPACE047';
%   'SPACE049';
%   'SPACE036';
%   };
% 
% % only one cell, with all session names
% sesNames = {'session_1'};
% 
% % replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
% replaceDataroot = true;
% 
% replaceDatatype = {'tla','pow'};
% 
% allowRecallSynonyms = true;
% 
% % accurateClassifSelect = true;
% accurateClassifSelect = false;

%% set up frequencies

freq = mm_freqSet('ndtools');

%% analysis details

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
cfg_ana.avgFrq = cfg_ft.avgoverfreq;

allConditions = {...
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
  cfg_ana.latencies = [0 0.48; 0.5 1.0];
elseif strcmp(cfg_ft.avgovertime,'yes')
  cfg_ana.latencies = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
%   cfg_ana.latencies = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
  
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end

if strcmp(cfg_ft.avgoverfreq,'no')
  % cfg_ana.frequencies = [3 80];
  cfg_ana.frequencies = freq.fullRange;

elseif strcmp(cfg_ft.avgoverfreq,'yes')
  % cfg_ana.frequencies = [3 7; 8 12; 8 10; 11 12; 13 20; 23 30; 32 47; 51 80];
  
  cfg_ana.frequencies = [ ...
    freq.theta; ...
    freq.alpha; ...
    freq.alpha_lower; ...
    freq.alpha_upper; ...
    freq.beta_lower; ...
    freq.beta_upper; ...
    freq.gamma_lower; ...
    freq.gamma_upper];
  
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgF'];
end

if strcmp(cfg_ft.avgoverchan,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgC'];
end

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
thisRun = [expName,'_overview_',datestr(now,'ddmmmyyyy-HHMMSS')];
diary(fullfile(saveDirProc,[thisRun '.log']));
tStart = tic;
fprintf('START TIME: %s\n',datestr(now,13));
for i = STAGES
  tS = tic;
  fprintf('STAGE%d START TIME: %s\n',i, datestr(now,13));
  
  % execute the processing stages
  if i == 1
    stageFun{i}(expName,saveDirProc,loadFile,cfg_ft,cfg_ana,allConditions,runLocally,timeOut{i});
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

function stage1(expName,saveDirProc,loadFile,cfg_ft,cfg_ana,allConditions,runLocally,timeOut)
% stage1: process the input files with FieldTrip based on the analysis
% parameters

%% Process the data
if runLocally == 0
  %% Dream: create one task for each subject (i.e., submit one to each node)
  
  % start a new job
  job = newJob([]);
  
  for i = 1:length(allConditions)
    % create one job per comparison
    cfg_ana.conditions = allConditions(i);
    fprintf('Processing cluster statistics for %s vs %s...\n',cfg_ana.conditions{1}{1},cfg_ana.conditions{1}{2});
    
    inArg = {saveDirProc,loadFile,cfg_ft,cfg_ana};
    createTask(job,@space_pow_clusterStat_cluster,0,inArg);
  end
  
  runJob(job,timeOut,fullfile(saveDirProc,[expName,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
  
  % final step: destroy the job because this doesn't happen in runJob
  destroy(job);
  
else
  %% run the function locally
  
  %error('This setup does not make sense with the current data loading setup.');
  
  % create a log of the command window output
  thisRun = [expName,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS')];
  % turn the diary on
  diary(fullfile(saveDirProc,[thisRun,'.log']));
  
  % use the peer toolbox
  %ana.usePeer = 1;
  %ana.usePeer = 0;
  
  cfg_ana.conditions = allConditions;
  space_pow_clusterStat_cluster(saveDirProc,loadFile,cfg_ft,cfg_ana);
  
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
  if ~isempty(dirs)
    set(job, 'PathDependencies', {homeDir, myMatlabDir, pwd(), p, dirs.dataroot});
  else
    set(job, 'PathDependencies', {homeDir, myMatlabDir, pwd(), p});
  end
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
