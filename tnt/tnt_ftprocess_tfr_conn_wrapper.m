function tnt_ftprocess_tfr_conn_wrapper(whichStages)
% tnt_ftprocess_tfr_pow_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub tnt_ftprocess_tfr_pow_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There is only one stage:
%  stage1 = call wrapper that calls create_ft_struct (which calls seg2ft,
%  which calls ft_freqanalysis) and saves one file per subject
%
% Input:
%  whichStages: the stage number(s) to run (default = 1)
%
% Output:
%  time-frequency data

% add eeg_directories to the path
%add_eeg_path

% check/handle arguments
error(nargchk(0,1,nargin))
STAGES = 1;
if nargin == 1
  STAGES = whichStages;
end

runLocally = 0;

% load analysisDetails instead of using specifications below
if runLocally == 0
  adFile = '/data/projects/curranlab/TNT_matt/eeg/-1000_1700/ft_data/B_NT_TH_eq1/conn_scd_mtmconvol_hanning_fourier_-500_980_3_9/analysisDetails.mat';
  [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
elseif runLocally == 1
  adFile = '/Volumes/curranlab/TNT/TNT_matt/eeg/-1000_1700/ft_data/B_NT_TH_eq1/conn_scd_mtmconvol_hanning_fourier_-500_980_3_9/analysisDetails.mat';
  [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
end

% % initialize the analysis structs
% exper = struct;
% files = struct;
% dirs = struct;
% ana = struct;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% MODIFY THIS STUFF
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Experiment-specific setup
% 
% exper.name = 'TNT';
% 
% exper.sampleRate = 250;
% 
% % pre- and post-stimulus times to read, in seconds (pre is negative)
% exper.prepost = [-1.0 1.7];
% 
% % equate the number of trials across event values?
% exper.equateTrials = 1;
% 
% % type of NS file for FieldTrip to read; raw or sbin must be put in
% % dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
% %exper.eegFileExt = 'raw';
% 
% % types of events to find in the NS file; these must be the same as the
% % events in the NS files
% exper.eventValues = sort({'B1of2','B2of2','NT1of2For','NT1of2Rec','NT2of2For','NT2of2Rec','T1of2For','T1of2Rec','T2of2For','T2of2Rec'});
% %exper.eventValues = sort({'NT1of2For','NT1of2Rec','NT2of2For','NT2of2Rec','T1of2For','T1of2Rec','T2of2For','T2of2Rec'});
% %exper.eventValues = sort({'NT1of2Rec','T1of2Rec'});
% 
% % combine some events into higher-level categories
% exper.eventValuesExtra.toCombine = {...
%   {'B1of2','B2of2'},...
%   {'NT1of2For','NT2of2For','NT1of2Rec','NT2of2Rec'},...
%   {'T1of2For','T2of2For','T1of2Rec','T2of2Rec'}...
% %   {'NT1of2For','NT1of2Rec'},{'NT2of2For','NT2of2Rec'}...
% %   {'T1of2For','T1of2Rec'},{'T2of2For','T2of2Rec'}...
% %   {'NT1of2For','NT2of2For'},{'NT1of2Rec','NT2of2Rec'}...
% %   {'T1of2For','T2of2For'},{'T1of2Rec','T2of2Rec'}...
%   };
% exper.eventValuesExtra.newValue = {...
%   {'B'},...
%   {'NT'},...
%   {'TH'}...
% %   {'NT1'},{'NT2'}...
% %   {'TH1'},{'TH2'}...
% %   {'NTF'},{'NTR'}...
% %   {'THF'},{'THR'}...
%   };
% 
% % keep only the combined (extra) events and throw out the original events?
% exper.eventValuesExtra.onlyKeepExtras = 1;
% 
% exper.subjects = {
%   'TNT 06';
%   'TNT 07';
%   'TNT 08';
%   'TNT 09';
%   'TNT 11';
%   'TNT 13';
%   'TNT 14';
%   'TNT 15';
%   'TNT 17';
%   'TNT 19';
%   'TNT 20';
%   'TNT 21';
%   'TNT 22';
%   'TNT 23';
%   'TNT 25';
%   'TNT 26';
%   'TNT 27';
%   'TNT 28';
%   'TNT 30';
%   'TNT 32';
%   'TNT 33';
%   'TNT 35';
%   'TNT 39';
%   'TNT 41';
%   'TNT 42';
%   'TNT 44';
%   'TNT 45';
%   'TNT 46';
%   'TNT 47';
%   'TNT 48';
%   'TNT 49';
%   'TNT 50';
%   'TNT 51';
%   'TNT 52';
%   'TNT 53';
%   'TNT 54';
%   };
% 
% % The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
% exper.sessions = {'session_0'};
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% POSSIBLY MODIFY THIS STUFF
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% set up file and directory handling parameters
% 
% % directory where the data to read is located
% dirs.dataDir = fullfile('TNT_matt','eeg',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));
% % Possible locations of the data files (dataroot)
% dirs.serverDir = fullfile('/Volumes','curranlab','Data');
% dirs.serverLocalDir = fullfile('/Volumes','RAID','curranlab','Data');
% dirs.localDir = fullfile(getenv('HOME'),'data');
% dirs.dreamDir = fullfile('/data/projects/curranlab',dirs.dataDir);
% 
% % pick the right dirs.dataroot
% if exist(dirs.serverDir,'dir')
%   dirs.dataroot = dirs.serverDir;
%   runLocally = 1;
% elseif exist(dirs.serverLocalDir,'dir')
%   dirs.dataroot = dirs.serverLocalDir;
%   runLocally = 1;
% elseif exist(dirs.localDir,'dir')
%   dirs.dataroot = dirs.localDir;
%   runLocally = 1;
% elseif exist(dirs.dreamDir,'dir')
%   dirs.dataroot = dirs.dreamDir;
%   runLocally = 0;
% else
%   error('Data directory not found.');
% end
% 
% % Use the FT chan locs file
% files.elecfile = 'GSN-HydroCel-129.sfp';
% files.locsFormat = 'besa_sfp';
% ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);
% 
% %% set up analysis parameters
% 
% ana.segFxn = 'seg2ft';
% ana.ftFxn = 'ft_freqanalysis';
% 
% % any preprocessing?
% cfg_pp = [];
% 
% cfg_proc = [];
% cfg_proc.output = 'pow';
% cfg_proc.pad = 'maxperlen';
% 
% %cfg_proc.output = 'powandcsd';
% % % channelcmb is set up as {'channel','cohref'} pairs
% %cfg_proc.channelcmb = {'E11','E62';'E11','E52';'E11','E92';'E11','E45';'E11','E108'}; % Fz, Pz; Fz, P3/P4; Fz, T7/T8
% % % cfg_proc.channelcmb = {'E3','E60';'E3','E62'};
% % % cfg_proc.channelcmb = {'E3','E60';'E3','E62';'E11','E60';'E11','E62'};
% %cfg_proc.channel = unique(cfg_proc.channelcmb);
% % % do not need to keep trials for powandcsd
% %cfg_proc.keeptrials = 'no';
% %cfg_proc.keeptapers = 'no';
% 
% 
% cfg_proc.output = 'fourier';
% cfg_proc.channelcmb = {'all','all'};
% % need to keep trials for fourier; not for powandcsd
% cfg_proc.keeptrials = 'yes';
% cfg_proc.keeptapers = 'yes';
% 
% 
% % % MTM FFT
% % cfg_proc.method = 'mtmfft';
% % cfg_proc.taper = 'dpss';
% % %cfg_proc.foilim = [3 50];
% % freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% % %cfg_proc.foi = 3:freqstep:50;
% % cfg_proc.foi = 3:freqstep:9;
% % cfg_proc.tapsmofrq = 5;
% % cfg_proc.toi = -0:0.04:1.0;
% 
% % multi-taper method
% cfg_proc.method = 'mtmconvol';
% cfg_proc.taper = 'hanning';
% %cfg_proc.taper = 'dpss';
% %cfg_proc.toi = -0.8:0.04:3.0;
% cfg_proc.toi = -0.5:0.04:1.0;
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% %cfg_proc.foi = 3:freqstep:50;
% cfg_proc.foi = 3:freqstep:9;
% %cfg_proc.foi = 3:1:9;
% %cfg_proc.foi = 2:2:30;
% cfg_proc.t_ftimwin = 4./cfg_proc.foi;
% % tapsmofrq is not used for hanning taper; it is used for dpss
% %cfg_proc.tapsmofrq = 0.4*cfg_proc.foi;
% 
% % % wavelet
% % cfg_proc.method = 'wavelet';
% % cfg_proc.width = 5;
% % %cfg_proc.toi = -0.8:0.04:3.0;
% % cfg_proc.toi = -0.3:0.04:1.0;
% % % evenly spaced frequencies, but not as many as foilim makes
% % freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% % %cfg_proc.foi = 3:freqstep:50;
% % cfg_proc.foi = 3:freqstep:9;
% % %cfg_proc.foilim = [3 9];
% 
% ana.ftype = cfg_proc.output;
% 
% % set the save directories
% [dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'conn');

%% set up for running stages and specifics for Dream

% name(s) of the functions for different stages of processing
stageFun = {@stage1};
timeOut  = {2}; % in HOURS

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
%thisRun = [exper.name,'_overview_',datestr(now,7) datestr(now,3) datestr(now,10)];
diary(fullfile(dirs.saveDirProc,[thisRun '.log']));
tStart = tic;
fprintf('START TIME: %s\n',datestr(now,13));
for i = STAGES
  tS = tic;
  fprintf('STAGE%d START TIME: %s\n',i, datestr(now,13));
  
  % execute the processing stage
  stageFun{i}(ana,exper,dirs,runLocally,timeOut{i});
  
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

function stage1(ana,exper,dirs,runLocally,timeOut)
% stage1: process the input files with FieldTrip based on the analysis
% parameters

%% Process the data
if runLocally == 0
  %% Dream: create one task for each subject (i.e., submit one to each node)
  
  % start a new job
  job = newJob(dirs);
  
  cfg_ft = [];
  cfg_ft.method = 'coh';
  %cfg_ft.method = 'plv';
  
  % do some pairwise combinations
  chans_to_pair = {'E1'; 'E9'; 'E15'; 'E22'; 'E32'; 'E11'; 'E33'; 'E24'; 'E124'; 'E122'; 'E120'; 'E116'; 'E111'; 'E6'; 'E29'; 'E34'; 'E43'; 'E108'; 'E103'; 'E87'; 'E55'; 'E37'; 'E41'; 'E45'; 'E95'; 'E97'; 'E84'; 'E86'; 'E62'; 'E72'; 'E53'; 'E66'; 'E51'; 'E64'; 'E89'; 'E75'; 'E69'; 'E81'; 'E52'; 'E92'};
  cfg_ft.channelcmb = nchoosek(chans_to_pair,2);
  
  % % do all pairwise combinations
  % chans_to_pair = cell(129,1);
  % for i = 1:length(chans_to_pair)
  %   chans_to_pair{i} = sprintf('E%d',i);
  % end
  % chans_to_pair{end} = 'Cz';
  % cfg_ft.channelcmb = nchoosek(chans_to_pair,2);
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      for evVal = 1:length(exper.eventValues)
        
        cfg_ft.inputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',ana.ftype,exper.eventValues{evVal}));
        cfg_ft.outputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',cfg_ft.method,exper.eventValues{evVal}));
        
        %if ~exist(cfg_ft.outputfile,'file')
        %fprintf('Processing %s, %s, %s...\n',exper.subjects{sub},exper.sessions{ses},exper.eventValues{evVal});
        %ft_connectivityanalysis(cfg_ft);
        %fprintf('Done.\n');
        
        inArg = {cfg_ft};
        createTask(job,@ft_connectivityanalysis,0,inArg);
        
        %else
        %  fprintf('ALREADY EXISTS: %s\n',cfg_ft.outputfile);
        %end
      end
    end
  end

  runJob(job, timeOut, fullfile(dirs.saveDirProc,[exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
  
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
  
  cfg_ft = [];
  cfg_ft.method = 'coh';
  %cfg_ft.method = 'plv';
  
  %cfg_ft.channelcmb = {'E11','E62';'E11','E52';'E11','E92';'E11','E45';'E11','E108'}; % Fz, Pz; Fz, P3/P4; Fz, T7/T8
  %cfg_ft.channelcmb = {'E3','E60';'E3','E62'};
  %cfg_ft.channelcmb = {'E3','E60';'E3','E62';'E11','E60';'E11','E62'};
  
  % do some pairwise combinations
  chans_to_pair = {'E1'; 'E6'; 'E9'; 'E11'; 'E15'; 'E22'; 'E24'; 'E29'; 'E32'; 'E33'; 'E34'; 'E37'; 'E41'; 'E43'; 'E45'; 'E51'; 'E52'; 'E53'; 'E55'; 'E62'; 'E64'; 'E66'; 'E69'; 'E72'; 'E75'; 'E81'; 'E84'; 'E86'; 'E87'; 'E89'; 'E92'; 'E95'; 'E97'; 'E103'; 'E108'; 'E111'; 'E116'; 'E120'; 'E122'; 'E124'};
  cfg_ft.channelcmb = nchoosek(chans_to_pair,2);
  
  % % do all pairwise combinations
  % chans_to_pair = cell(129,1);
  % for i = 1:length(chans_to_pair)
  %   chans_to_pair{i} = sprintf('E%d',i);
  % end
  % chans_to_pair{end} = 'Cz';
  % cfg_ft.channelcmb = nchoosek(chans_to_pair,2);
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      for evVal = 1:length(exper.eventValues)
        
        cfg_ft.inputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',ana.ftype,exper.eventValues{evVal}));
        cfg_ft.outputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',cfg_ft.method,exper.eventValues{evVal}));
        
        %if ~exist(cfg_ft.outputfile,'file')
        %  fprintf('Processing %s, %s, %s...\n',exper.subjects{sub},exper.sessions{ses},exper.eventValues{evVal});
        ft_connectivityanalysis(cfg_ft);
        fprintf('Done.\n');
        %else
        %  fprintf('ALREADY EXISTS: %s\n',cfg_ft.outputfile);
        %end
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

% function add_eeg_path
% 
% %% add specific toolboxes to the path - you could instead put this in ~/Documents/MATLAB/startup.m
% 
% %% initialize
% homeDir = getenv('HOME');
% 
% myMatlabDir = fullfile(homeDir,'Documents','MATLAB');
% 
% %% set up eegToolbox path
% eegToolboxDir = dir(fullfile(myMatlabDir,'eegToolbox*'));
% if ~isempty(eegToolboxDir)
%   eegToolboxDir = fullfile(myMatlabDir,eegToolboxDir.name);
%   % add top folder and all subfolders
%   addpath(genpath(eegToolboxDir));
% end
% 
% %% set up eeglab path
% eeglabDir = dir(fullfile(myMatlabDir,'eeglab*'));
% if ~isempty(eeglabDir)
%   eeglabDir = fullfile(myMatlabDir,eeglabDir.name);
%   % add top folder and all subfolders
%   addpath(genpath(eeglabDir));
% 
%   % remove eeglab's external directory if it was added
%   eeglabExtDir = fullfile(eeglabDir,'external');
%   if ~isempty(eeglabExtDir)
%     rmpath(genpath(eeglabExtDir));
%   end
%   
%   % % remove eeglab's fieldtrip directory if it was added
%   % eeglabFtDir = dir(fullfile(eeglabDir,'external','fieldtrip*'));
%   % if ~isempty(eeglabFtDir)
%   %   eeglabFtDir = fullfile(myMatlabDir,eeglabFtDir.name);
%   %   rmpath(genpath(eeglabFtDir));
%   % end
% end
% 
% %% set up fieldtrip path
% ftDir = dir(fullfile(myMatlabDir,'fieldtrip*'));
% if ~isempty(ftDir)
%   ftDir = fullfile(myMatlabDir,ftDir.name);
%   % add only the top folder
%   addpath(ftDir);
%   % run ft_defaults to add the subdirectories that FT needs
%   ft_defaults
%   
%   % add the peer directory
%   addpath(fullfile(ftDir,'peer'));
%   
%   % add the SPM directory
%   %addpath(fullfile(ftDir,'external','spm8'));
%   
%   % % remove fieldtrip's external directory
%   % ftExtDir = fullfile(ftDir,'external');
%   % if ~isempty(ftExtDir)
%   %   rmpath(genpath(ftExtDir));
%   % end
%   
%   % remove fieldtrip's eeglab directory if it was added
%   ftEeglabDir = dir(fullfile(ftDir,'external','eeglab*'));
%   if ~isempty(ftEeglabDir)
%     ftEeglabDir = fullfile(myMatlabDir,ftEeglabDir.name);
%     rmpath(genpath(ftEeglabDir));
%   end
% end
% 
% %% set up EP_Toolkit path
% epDir = dir(fullfile(myMatlabDir,'EP_Toolkit*'));
% if ~isempty(epDir)
%   epDir = fullfile(myMatlabDir,epDir.name);
%   % add top folder and all subfolders
%   addpath(genpath(epDir));
% end
% 
% %% add my other analysis scripts/files
% addpath(genpath(fullfile(myMatlabDir,'recogmodel_mvm')));
% 
% %% add my experiment, fieldtrip, and RM ANOVA scripts
% addpath(genpath(fullfile(myMatlabDir,'mat_mvm')));
% 
% %% put the ~/Documents/MATLAB folder at the top of the path
% addpath(myMatlabDir);
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
% 
