function blink_correction_wrapper(whichStages)
% blink_correction_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub blink_correction_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There is only one stage:
%  stage1 = ICA blink correction using EP Toolkit
%
% Input:
%  whichStages: the stage number(s) to run (default = 1)
%
% Output:
%  Currently outputs to ep_mat format for doing analyses in EP Toolkit
%
%  If you use another format:
%
%    egis files need to have type/creator changed to EGIS/NETs
%
%    raw files, need to have type/creator changed to eGLY/NETs
%

% check/handle arguments
error(nargchk(0,1,nargin))
STAGES = 1;
if nargin == 1
  STAGES = whichStages;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exper.name = 'SOSI';
% %exper.name = 'SOCO';
% %exper.name = 'COSI';
% dirs.eegDir = fullfile('eeg','eppp');
% % sampling rate of the EEG data
% exper.sampleRate = 250;
% % segment length created with NetStation (in seconds)
% exper.prepost = [-1 2];
% % baseline period (in milliseconds)
% exper.baseline_ms = [-200 0];

exper.name = 'COSI2';
dirs.eegDir = fullfile('eeg','eppp');
% sampling rate of the EEG data
exper.sampleRate = 500;
% segment length created with NetStation (in seconds)
exper.prepost = [-1 2];
% baseline period (in milliseconds)
exper.baseline_ms = [-200 0];

% exper.name = 'FRCE';
% dirs.eegDir = fullfile('EEG','eppp');
% % sampling rate of the EEG data
% exper.sampleRate = 500;
% % segment length created with NetStation (in seconds)
% exper.prepost = [-1.25 2.25];
% % baseline period (in milliseconds)
% exper.baseline_ms = [-200 0];

% Type of input file
exper.inputFileExt = 'egis';
%exper.inputFileExt = 'raw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POSSIBLY MODIFY THIS STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% name(s) of the functions for different stages of processing
stageFun = {@stage1};
timeOut  = {2}; % in HOURS

% Location of the data files (dataroot)s
if strcmp(exper.inputFileExt,'egis') || strcmp(exper.inputFileExt,'raw')
  dirs.dataDir = fullfile(dirs.eegDir,sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),sprintf('2_ns_%s',exper.inputFileExt));
else
  error('Filetype %s is unknown.',exper.inputFileExt);
end

dirs.homeDir = getenv('HOME');

% 2 potential Curran server dataroots
dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data',exper.name);
dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data',exper.name);
dirs.localDir = fullfile(dirs.homeDir,'data',exper.name);
% Dream dataroot
dirs.dreamDir = fullfile(filesep,'data','projects','curranlab',exper.name);

% pick the right dirs.dataroot
if exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
  runLocally = 1;
elseif exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
  runLocally = 1;
elseif exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
  runLocally = 1;
elseif exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
  runLocally = 0;
else
  error('Data directory not found.');
end

if runLocally == 0
  %% set up for running on Dream
  
  % need to export DISPLAY to an offscreen buffer for MATLAB DCS graphics
  sched = findResource();
  if strcmp(sched.Type, 'generic')
    setenv('DISPLAY', 'dream:99');
  end

  % add the EEG analysis toolboxes to the path if not in startup.m
  %add_eeg_path;
  
  % number of time points to read in (approx 100,000 per 1 GB of RAM)
  exper.memChunkSize = 3 * 100000;
  
else
  %% set up for running locally
  
  % number of time points to read in (approx 100,000 per 1 GB of RAM)
  
  % RAM for the Mac Pro server
  %exper.memChunkSize = 24 * 100000;
  % RAM for the middle Mac Pro
  exper.memChunkSize = 18 * 100000;
  % RAM for my laptop
  %exper.memChunkSize = 8 * 100000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%capture diary and time statistics
thisRun = [exper.name,'_overview_',datestr(now,'ddmmmyyyy-HHMMSS')];
%thisRun = [exper.name,'_overview_',datestr(now,7) datestr(now,3) datestr(now,10)];
diary(fullfile(dirs.dataroot,dirs.dataDir,[thisRun '.log']));
tStart = tic;
fprintf('START TIME: %s\n',datestr(now,13));
for i = STAGES
  tS = tic;
  fprintf('STAGE%d START TIME: %s\n',i, datestr(now,13));
  
  % execute the processing stage
  stageFun{i}(exper, dirs, runLocally, timeOut{i});
  
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

function stage1(exper, dirs, runLocally, timeOut)
% stage1: run ICA blink correction using EP Toolkit

if runLocally == 0
  %% start a new job
  job = newJob(dirs);
end

%% set up analysis parameters

if ~isempty(exper.baseline_ms)
  % calculate the baseline in samples
  exper.baseline_samp = abs(exper.prepost(1) * exper.sampleRate) + [((exper.baseline_ms(1) / (1000/exper.sampleRate)) + 1) (exper.baseline_ms(2) / (1000/exper.sampleRate))];
elseif isempty(exper.baseline_ms)
  exper.baseline_samp = [];
end
% struct for ep_readData, which gets called by ep_artifactCorrection
exper.textPrefs_struct = struct('firstRow',1,'firstCol',1,'lastCol',0,'orientation',1);

% initialize
exper.inputFile = '';
exper.blinkFile = '';

if strcmp(exper.inputFileExt,'egis')
  % EGIS format
  exper.inputFormat = 'egi_egis';
  %exper.outputFormat = 'egi_egis';
  exper.outputFormat = 'ep_mat';
elseif strcmp(exper.inputFileExt,'raw')
  % RAW format
  exper.inputFormat = 'egi_sbin';
  %exper.outputFormat = 'egi_sbin';
  exper.outputFormat = 'ep_mat';
end

exper.type = 'single_trial';
exper.montage = 'Hydrocel-128-1';
exper.ced = 'GSN-Hydrocel-129.ced';
exper.refChan = 129;

%% set up input arguments
inArg = {'files',exper.inputFile,...
  'inputFormat',exper.inputFormat,...
  'outputFormat',exper.outputFormat,...
  'baseline',exper.baseline_samp,...
  'timePoints',[],... % default: []. don't drop any timepoints
  'template','bothTemplate',...
  'sacctemplate','none',...
  'blinkFile',exper.blinkFile,...
  'saturation',[-Inf Inf],... % default: [-Inf Inf]
  'window',80,... % default: 80 ms
  'minmax',100,... % bad chan min max (default: 100 uV)
  'badnum',10,... % percent bad chan exceeded to declare bad trial (default: 10)
  'neighbors',6,... % default: 6
  'maxneighbor',30,... % default: 30 uV
  'badchan',0.4,... % minimum predictability from neighbors (default: 0.4)
  'blink',0.9,... % threshold correlation with blink template, 0 to 1 (default: .9)
  'trialminmax',100,... % bad trial voltage diff; default 100 uV
  'detrend',0,... % not recommended (default: 0)
  'badtrials',20,... % default 20%
  'channelMode','replace',... % bad channel correction: 'replace' interpolates bad channels; 'none' to do nothing
  'trialMode','none',... % movement correction: 'fix' to fix bad trial data; 'none' to do nothing
  'noadjacent',1,... % 1 to not allow adjacent bad channels (trial or subject declared bad) (default: 1)
  'chunkSize',exper.memChunkSize,...
  'minTrialsPerCell',15,...
  'movefacs',20,...
  'noFigure',0,...
  'reference',exper.refChan,... % channel 129 is Cz, implicitly present
  'textPrefs',exper.textPrefs_struct,...
  'type',exper.type,...
  'montage',exper.montage,...
  'ced',exper.ced,...
  };

%% find the files to be processed
fileList = dir(fullfile(dirs.dataroot,dirs.dataDir,[exper.name,'*.',exper.inputFileExt]));
if isempty(fileList)
  error('No files found!')
else
  fprintf('Processing %d file(s) found in %s\n',length(fileList),fullfile(dirs.dataroot,dirs.dataDir));
end

% find the blink file
exper.blinkFile = dir(fullfile(dirs.dataroot,dirs.dataDir,'*blink*.mat'));
if isempty(exper.blinkFile)
  error('No blink file found!')
else
  fprintf('Using blink file: %s\n',fullfile(dirs.dataroot,dirs.dataDir,exper.blinkFile.name));
end

% find where in the cell array to put the input file name
for i = 1:length(inArg)
  if ischar(inArg{i})
    if strcmp(inArg{i},'files')
      fileIndex = i + 1;
      break
    end
  end
end
% find where in the cell array to put the blink file name
for i = 1:length(inArg)
  if ischar(inArg{i})
    if strcmp(inArg{i},'blinkFile')
      blinkIndex = i + 1;
      break
    end
  end
end
inArg{blinkIndex} = fullfile(dirs.dataroot,dirs.dataDir,exper.blinkFile.name);

%% Process the data
if runLocally == 0
  %% Dream: create one task for each subject (i.e., submit one to each node)
  for i = 1:length(fileList)
    fprintf('Processing %s...\n',fileList(i).name);
    
    % set the input file name
    inArg{fileIndex} = {fullfile(dirs.dataroot,dirs.dataDir,fileList(i).name)};
    
    % Dream: create one task for each subject
    createTask(job,@ep_artifactCorrection,0,inArg);
  end
  
  runJob(job, timeOut, fullfile(dirs.dataroot,dirs.dataDir,[exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
else
  %% run ep_artifactCorrection locally, submitting all subjects at once
  
  % Store all the file names in 1 big cell array
  inputFiles = cell(1,length(fileList));
  for i = 1:length(fileList)
    inputFiles{i} = fullfile(dirs.dataroot,dirs.dataDir,fileList(i).name);
  end
  % set the input files
  inArg{1,fileIndex} = inputFiles;
  
  % create a log of EP's command window output
  thisRun = [exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS')];
  diary(fullfile(dirs.dataroot,dirs.dataDir,[thisRun,'.log']));
  
  % run ep_artifactCorrection
  ep_artifactCorrection(inArg);
  
  % turn off the diary
  diary off
  
%   %% run ep_artifactCorrection locally, submitting one subject at a time
%   for i = 1:length(fileList)
%     fprintf('Processing %s...\n',fileList(i).name);
%     
%     % set the input file name
%     inArg{fileIndex} = {fullfile(dataroot,dirs.dataDir,fileList(i).name)};
%     
%     % run ep_artifactCorrection
%     ep_artifactCorrection(inArg);
%   end
end

function job = newJob(dirs)
% newJob Creates a new PCT job and sets job's dependencies
%
%   exper -- data structure with fields for subject, session, conditions,
%              data locations, etc

% Set up scheduler, job
sched = findResource();
job = createJob(sched);
% define the directories to add to worker sessions' matlab path
myMatlabDir = fullfile(dirs.homeDir,'Documents','MATLAB');
p = path();
set(job, 'PathDependencies', {dirs.homeDir, myMatlabDir, pwd(), p, fullfile(dirs.dataroot,dirs.dataDir)});

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
  if ~errors
    destroy(job);
  else
    error([mfilename ':logOutput'],'%s had %d errors',job.Name, errors)
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

% function add_eeg_path
% 
%   %% add specific toolboxes to the path - you could instead put this in
%   %% ~/Documents/MATLAB/startup.m
%   
%   % set paths for EP Toolkit, FieldTrip, EEGLAB, etc.
% homeDir = getenv('HOME');
% 
% myMatlabDir = fullfile(homeDir,'Documents','MATLAB');
% 
% %% set up eeglab path
% eeglab_dir = dir(fullfile(myMatlabDir,'eeglab*'));
% if ~isempty(eeglab_dir)
%   eeglab_dir = fullfile(myMatlabDir,eeglab_dir.name);
%   % add top folder and all subfolders
%   addpath(genpath(eeglab_dir));
% 
%   % remove eeglab's external directory if it was added
%   eeglab_ext_dir = fullfile(eeglab_dir,'external');
%   if ~isempty(eeglab_ext_dir)
%     rmpath(genpath(eeglab_ext_dir));
%   end
%   
%   % % remove eeglab's fieldtrip directory if it was added
%   % eeglab_ft_dir = dir(fullfile(eeglab_dir,'external','fieldtrip*'));
%   % if ~isempty(eeglab_ft_dir)
%   %   eeglab_ft_dir = fullfile(myMatlabDir,eeglab_ft_dir.name);
%   %   rmpath(genpath(eeglab_ft_dir));
%   % end
% end
% 
% %% set up fieldtrip path
% ft_dir = dir(fullfile(myMatlabDir,'fieldtrip*'));
% if ~isempty(ft_dir)
%   ft_dir = fullfile(myMatlabDir,ft_dir.name);
%   % add only the top folder
%   addpath(ft_dir);
%   % run ft_defaults to add the subdirectories that FT needs
%   ft_defaults
%   
%   % % add the peer directory
%   % addpath(fullfile(ft_dir,'peer'));
%   
%   % % add the SPM directory
%   % addpath(fullfile(ft_dir,'external','spm8'));
%   
%   % % remove fieldtrip's external directory
%   % ft_ext_dir = fullfile(ft_dir,'external');
%   % if ~isempty(ft_ext_dir)
%   %   rmpath(genpath(ft_ext_dir));
%   % end
%   
%   % remove fieldtrip's eeglab directory if it was added
%   ft_eeglab_dir = dir(fullfile(ft_dir,'external','eeglab*'));
%   if ~isempty(ft_eeglab_dir)
%     ft_eeglab_dir = fullfile(myMatlabDir,ft_eeglab_dir.name);
%     rmpath(genpath(ft_eeglab_dir));
%   end
% end
% 
% %% set up EP_Toolkit path
% ep_dir = dir(fullfile(myMatlabDir,'EP_Toolkit*'));
% if ~isempty(ep_dir)
%   ep_dir = fullfile(myMatlabDir,ep_dir.name);
%   % add top folder and all subfolders
%   addpath(genpath(ep_dir));
% end
% 
% %% add my experiment, fieldtrip, and RM ANOVA scripts
% addpath(genpath(fullfile(myMatlabDir,'mat_mvm')));
% 
% %% add my other analysis scripts
% %addpath(genpath(fullfile(myMatlabDir,'recogmodel_mvm')));
% %addpath(fullfile(homeDir,'eeg'));
% 
% %% remove CVS and .svn directories from path
% entries = regexp(path, ['[^',pathsep,']*',pathsep], 'match');
% for i = 1:length(entries)
%   entry = char(entries{i});
%   if ~isempty(strfind(entry, '.svn'))
%     rmpath(entry);
%   end
%   if ~isempty(strfind(entry, 'CVS'))
%     rmpath(entry);
%   end
% end
