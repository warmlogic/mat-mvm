function space_rsa_pca_tla_hilbert_classif_wrapper(whichStages)
% space_rsa_pca_tla_classif_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub space_rsa_pca_tla_classif_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There is one stage:
%  stage1 =
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expName = 'SPACE';

% saveDirProc = getenv('HOME');
% saveDirProc = fullfile(filesep,'data','projects','curranlab',expName);
saveDirProc = fullfile(filesep,'data','projects','curranlab',expName,'EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla');

runLocally = 0;

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
  %'SPACE017'; % old assessment: really noisy EEG, half of ICA components rejected
  'SPACE018';
  %'SPACE019';
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
  %'SPACE030'; % low trial counts
  'SPACE032';
  'SPACE034';
  'SPACE047';
  'SPACE049';
  'SPACE036';
  };

% only one cell, with all session names
sesNames = {'session_1'};

%% analysis details

allROIs = {{'LPI2','LPS','LT','RPI2','RPS','RT'},{'center109'},{'LPS','RPS'},{'LT','RT'},{'LPI2','RPI2'}};

allLats = {[0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
  0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
  0 0.3; 0.3 0.6; 0.6 0.9; ...
  0 0.5; 0.5 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0]};

% allFreqs = {[4 8; 8 12; 12 30; 30 50]};
% allFreqs = {[4 8] [8 12] [12 30] [30 50]};
allFreqs = {[4 8; 8 12; 12 30; 30 50] [4 8] [8 12] [12 30] [30 50]};

allSimMethod = {'cosine'};

% sim_method = 'cosine';
% sim_method = 'correlation';
% sim_method = 'spearman';
% sim_method = 'euclidean';

% allEigCrit = {'kaiser'};
allEigCrit = {'CV85','kaiser'};

% % keep components with eigenvalue >= 1
% eig_criterion = 'kaiser';

% % compute the percent explained variance expected from each component if
% % all events are uncorrelated with each other; keep it if above this level.
% % So, each component would explain 100/n, where n is the number of
% % events/components.
% eig_criterion = 'analytic';

% keep components that cumulatively explain at least 85% of the variance
% eig_criterion = 'CV85';

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
    stageFun{i}(expName,saveDirProc,subjects,sesNames,allROIs,allLats,allFreqs,allSimMethod,allEigCrit,runLocally,timeOut{i});
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

function stage1(expName,saveDirProc,subjects,sesNames,allROIs,allLats,allFreqs,allSimMethod,allEigCrit,runLocally,timeOut)
% stage1: process the input files with FieldTrip based on the analysis
% parameters

%% Process the data
if runLocally == 0
  %% Dream: create one task for each subject (i.e., submit one to each node)
  
  % start a new job
  job = newJob([]);
  
  % save the original subjects array so we can set exper to have single
  % subjects, one for each task created
  %allSubjects = subjects;
  
  for i = 1:length(subjects)
    % Dream: create one task for each subject
    thisSub = subjects{i};
    thisSes = sesNames{1};
    
    for r = 1:length(allROIs)
      thisROI = allROIs{r};
      if iscell(thisROI)
        roi_str = sprintf(repmat('%s',1,length(thisROI)),thisROI{:});
      elseif ischar(thisROI)
        roi_str = thisROI;
      end
      
      for l = 1:length(allLats)
        latencies = allLats{l};
        
        for f = 1:length(allFreqs)
          freqs = allFreqs{f};
          freq_str = sprintf('%dfreq%dto%d',size(freqs,1),freqs(1,1),freqs(end,end));
          
          for s = 1:length(allSimMethod)
            sim_method = allSimMethod{s};
            
            for e = 1:length(allEigCrit)
              eig_criterion = allEigCrit{e};
              
              fprintf('Processing %s %s, roi: %s, %d latencies, %s, sim: %s, eig: %s...\n',thisSub,thisSes,roi_str,size(latencies,1),freq_str,sim_method,eig_criterion);
              
              % inArg = {ana,cfg_pp,exper,dirs,files};
              % save the exper struct (output 1) so we can use it later
              %createTask(job,@create_ft_struct_multiSes,1,inArg);
              
              inArg = {thisSub,thisSes,thisROI,latencies,freqs,sim_method,eig_criterion};
              createTask(job,@space_rsa_pca_tla_hilbert_classif_cluster,0,inArg);
            end
          end
        end
      end
    end
  end
  
  runJob(job,timeOut,fullfile(saveDirProc,[expName,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
  
  % final step: destroy the job because this doesn't happen in runJob
  destroy(job);
  
else
  %% run the function locally
  
  error('This setup does not make sense with the current data loading setup.');
  
  % create a log of the command window output
  thisRun = [expName,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS')];
  % turn the diary on
  diary(fullfile(saveDirProc,[thisRun,'.log']));
  
  % use the peer toolbox
  %ana.usePeer = 1;
  %ana.usePeer = 0;
  
  for r = 1:length(allROIs)
    thisROI = allROIs{r};
    for l = 1:length(allLats)
      latencies = allLats{l};
      for f = 1:length(allFreqs)
        freqs = allFreqs{f};
        for s = 1:length(allSimMethod)
          sim_method = allSimMethod{s};
          for e = 1:length(allEigCrit)
            eig_criterion = allEigCrit{e};
            % Local: run all the subjects
            space_rsa_pca_tla_hilbert_classif_cluster(subjects,sesNames,thisROI,latencies,freqs,sim_method,eig_criterion);
          end
        end
      end
    end
  end
  
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
