function space2_rsa_pca_tla_classif_wrapper(whichStages)
% space2_rsa_pca_tla_classif_wrapper(whichStages)
%
% This is the first file to run for doing PCA-based RSA on dream.
%
% To run on dream, at the command line type: distmsub2011 space2_rsa_pca_tla_classif_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There is one stage:
%  stage1 = run space2_rsa_pca_tla_classif_cluster.m
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

expName = 'SPACE2';

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

saveDirProc = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');
% saveDirProc = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_continuousICA_ftAuto/tla');

subjects = {
  'SPACE2001';
  %'SPACE2002'; % really noisy EEG, finished both sessions, incl in beh
  %'SPACE2003'; % DNF session 2, exclude
  'SPACE2004'; % low trial counts, exclude via threshold
  'SPACE2005';
  'SPACE2006';
  %'SPACE2007'; % bad performance, low trial counts, exclude
  'SPACE2008';
  %'SPACE2009'; % DNF session 2, exclude
  'SPACE2010';
  'SPACE2011';
  'SPACE2012';
  %'SPACE2013'; % didn't record EEG, stopped session 1 in middle, exclude
  'SPACE2014';
  'SPACE2015';
  'SPACE2016'; % low trial counts, exclude via threshold
  'SPACE2017'; % not great performance, still including
  'SPACE2018';
  'SPACE2019';
  %'SPACE2020'; % DNF session 2, exclude
  'SPACE2021';
  'SPACE2022';
  %'SPACE2023'; % DNF session 2, exclude
  %'SPACE2024'; % DNF session 2, exclude
  %'SPACE2025'; % bad performance, low trial counts, exclude
  'SPACE2026';
  %'SPACE2027'; % really noisy EEG, DNF session 2, exclude
  'SPACE2028';
  'SPACE2029'; % low trial counts, exclude via threshold
  'SPACE2030';
  'SPACE2031';
  'SPACE2032'; % low trial counts, exclude via threshold
  'SPACE2033';
  'SPACE2034'; % low trial counts, exclude via threshold
  'SPACE2035'; % not great performance, still including
  'SPACE2036';
  'SPACE2037';
  'SPACE2038';
  %'SPACE2039'; % DNF session 2, exclude
  'SPACE2040';
  };

% only one cell, with all session names
sesNames = {'session_1_session_2'};

% Correct (1) and incorrect (0) are for typos and total intrusions,
% respectively. Coupled (2) means intrinsically linked words and must have
% the same word stem (staple/stapler, bank/banker, dance/dancer,
% serve/server)  Synonym (3) means words that strictly have the same
% meaning and can have the same word stem (sofa/couch, doctor/physician,
% home/house, pasta/noodle, woman/lady, cash/money). Homonym (4) means
% words that sound exactly the same (board/bored, brake/break). Related (5)
% means closely associated words but cannot have the same word stem
% (whiskey/rum, map/compass, sailor/boat, broccoli/carrot).

% gradeCorrect = '[1, 2, 3, 4]';
gradeCorrect = '[1, 2, 4]';

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

replaceDatatype = {};

% accurateClassifSelect = true;
accurateClassifSelect = false;

%% analysis details

date_string = datestr(now,1);

lpfilt = false;
subselect_eeg = true;

dataTypes = {'img_rc_mass','img_fo_mass', 'img_rc_spac2','img_fo_spac2', 'img_rc_spac12','img_fo_spac12', 'img_rc_spac32','img_fo_spac32'};

% dataTypes = {'word_rc_mass','word_fo_mass', 'word_rc_spac2','word_fo_spac2', 'word_rc_spac12','word_fo_spac12', 'word_rc_spac32','word_fo_spac32'};

% dataTypes = {'img_rc_mass','img_fo_mass', 'img_rc_spac2','img_fo_spac2', 'img_rc_spac12','img_fo_spac12', 'img_rc_spac32','img_fo_spac32', 'word_rc_mass','word_fo_mass', 'word_rc_spac2','word_fo_spac2', 'word_rc_spac12','word_fo_spac12', 'word_rc_spac32','word_fo_spac32'};

% allROIs = {{'LPI2','LPS','LT','RPI2','RPS','RT'},{'center109'},{'LPS','RPS'},{'LT','RT'},{'LPI2','RPI2'},{'LAS2','FS','RAS2'},{'LFP','FC','RFP'}};
allROIs = {{'FS'},{'C'},{'PS'},{'PI'},{'LT'},{'RT'},{'LPS'},{'RPS'}};

% allLats = {[0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
%   0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
%   0 0.3; 0.3 0.6; 0.6 0.9; ...
%   0 0.5; 0.5 1.0; ...
%   0.3 0.8; ...
%   0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
%   0 0.8; 0.1 0.9; 0.2 1.0;
%   0 1.0]};

allLats = {[0.0 0.2; 0.22 0.4; 0.42 0.6; 0.62 0.8; 0.82 1.0; ...
  0.1 0.3; 0.32 0.5; 0.52 0.7; 0.72 0.9; ...
  0 0.3; 0.32 0.6; 0.62 0.9; ...
  0 0.5; 0.52 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0]};

cfg_sel = [];
cfg_sel.avgoverchan = 'no';
% cfg_sel.avgovertime = 'yes';
cfg_sel.avgovertime = 'no';

allSimMethod = {'cosine'};

% sim_method = 'cosine';
% sim_method = 'correlation';
% sim_method = 'spearman';
% sim_method = 'euclidean';

allEigCrit = {'kaiser'};
% allEigCrit = {'CV85','kaiser'};

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
    stageFun{i}(expName,saveDirProc,replaceDataroot,replaceDatatype,gradeCorrect,accurateClassifSelect,date_string,dataTypes,cfg_sel,subjects,sesNames,allROIs,allLats,allSimMethod,allEigCrit,lpfilt,subselect_eeg,runLocally,timeOut{i});
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

function stage1(expName,saveDirProc,replaceDataroot,replaceDatatype,gradeCorrect,accurateClassifSelect,date_string,dataTypes,cfg_sel,subjects,sesNames,allROIs,allLats,allSimMethod,allEigCrit,lpfilt,subselect_eeg,runLocally,timeOut)
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
        
        for s = 1:length(allSimMethod)
          sim_method = allSimMethod{s};
          
          for e = 1:length(allEigCrit)
            eig_criterion = allEigCrit{e};
            
            fprintf('Processing %s %s, roi: %s, %d latencies, sim: %s, eig: %s...\n',thisSub,thisSes,roi_str,size(latencies,1),sim_method,eig_criterion);
            
            inArg = {saveDirProc,replaceDataroot,replaceDatatype,gradeCorrect,accurateClassifSelect,date_string,dataTypes,cfg_sel,thisSub,thisSes,thisROI,latencies,sim_method,eig_criterion,lpfilt,subselect_eeg};
            createTask(job,@space2_rsa_pca_tla_classif_cluster,0,inArg);
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
  
  error('This setup does not make sense with the current data loading setup because it will try to reload the data on each iteration.');
  
  % create a log of the command window output
  thisRun = [expName,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS')];
  % turn the diary on
  diary(fullfile(saveDirProc,[thisRun,'.log']));
  
  % use the peer toolbox
  %ana.usePeer = 1;
  %ana.usePeer = 0;
  
  % Local: run all the subjects
  for r = 1:length(allROIs)
    thisROI = allROIs{r};
    for l = 1:length(allLats)
      latencies = allLats{l};
      
      for s = 1:length(allSimMethod)
        sim_method = allSimMethod{s};
        for e = 1:length(allEigCrit)
          eig_criterion = allEigCrit{e};
          space2_rsa_pca_tla_classif_cluster(saveDirProc,replaceDataroot,replaceDatatype,gradeCorrect,accurateClassifSelect,date_string,dataTypes,cfg_sel,subjects,sesNames,thisROI,latencies,sim_method,eig_criterion);
        end
      end
    end
  end
  
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
