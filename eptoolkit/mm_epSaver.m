function mm_epSaver(experName,prepost,eegDir,inputFormat,outputFormat)
%MM_EPSAVER: Use EP's read/write to save EEG data
%
% mm_epSaver(experName,prepost,eegDir,inputFormat,outputFormat)
%
% This is a simple input/output wrapper to, e.g., save ep_mat format to an
% EGI-readable format (egi_egis, egi_sbin) using ep_writeData.m
%
% Input:
%
%   experName    = experiment name (required; e.g., 'COSI')
%
%   prepost      = pre- and post-stimulus timing, in seconds
%                  (required; e.g., [-1.0 2.0])
%                  NB: NEGATIVE AND POSITIVE VALUES IN SECONDS!
%
%   eegDir       = where the 'eppp' directory is, under the experName dir
%                  (default: fullfile('eeg','eppp'))
%
%   inputFormat  = format to read (default: 'ep_mat')
%
%   outputFormat = format to save (default: 'egi_egis')
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original files should be located in:
%    dataroot/[experName]/[eeg/eppp]/[pre_post_ms]/3_ep_[inputFormatExtension]
%    e.g.,
%    dataroot/COSI/eeg/eppp/-1000_2000/3_ep_mat
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output:
%
%   Files are saved to 
%     dataroot/[experName]/[eeg/eppp]/[pre_post_ms]/3_ep_[outputFormatExtension]_out
%     e.g.,
%     dataroot/COSI/eeg/eppp/-1000_2000/3_ep_egis_out
%
% See also: EP_WRITEDATA
%


%experName = 'COSI2';

% experName = 'SOCO';
% prepost = [-1 2];
% eegDir = fullfile('eeg','eppp');

% experName = 'FRCE';
% prepost = [-1.25 2.25];
% eegDir = fullfile('EEG','Sessions','cueing paradigm','relabeled','eppp');

if ~exist('experName','var') || isempty(experName)
  error('Must set experName string (e.g., ''EXPR'')');
end

if ~exist('prepost','var') || isempty(prepost)
  %prepost = [-1.0 2.0];
  error('Must set prepost (e.g., [-1.0 2.0])');
end

if ~exist('eegDir','var') || isempty(eegDir)
  eegDir = fullfile('eeg','eppp');
  fprintf('Setting eegDir to default: %s\n',eegDir);
end

if ~exist('inputFormat','var') || isempty(inputFormat)
  inputFormat = 'ep_mat';
  fprintf('Setting inputFormat to default: %s\n',inputFormat);
end

if ~exist('outputFormat','var') || isempty(outputFormat)
  outputFormat = 'egi_egis';
  fprintf('Setting outputFormat to default: %s\n',outputFormat);
end


% Input file extension
if strcmp(inputFormat,'egi_egis')
  inputFileExt = 'egis';
elseif strcmp(inputFormat,'egi_sbin')
  inputFileExt = 'raw';
elseif strcmp(inputFormat,'ep_mat')
  inputFileExt = 'ept';
else
  error('Input filetype %s is unknown.',inputFormat);
end

% Location of the data files (dataroot)s
dirs.baseDir = fullfile(eegDir,sprintf('%d_%d',prepost(1)*1000,prepost(2)*1000));
dirs.dataDir = fullfile(dirs.baseDir,sprintf('3_ep_%s',inputFileExt));

dirs.homeDir = getenv('HOME');

% 2 potential Curran server dataroots
dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data',experName);
dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data',experName);
dirs.localDir = fullfile(dirs.homeDir,'data',experName);
% Dream dataroot
dirs.dreamDir = fullfile(filesep,'data','projects','curranlab',experName);

% pick the right dirs.dataroot
if exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
elseif exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
elseif exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
elseif exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
else
  error('Data directory not found.');
end

% see which files are in the data directory
if ~exist(fullfile(dirs.dataroot,dirs.dataDir),'dir')
  error('%s does not exist',fullfile(dirs.dataroot,dirs.dataDir));
else
  fileList = dir(fullfile(dirs.dataroot,dirs.dataDir,[experName,'*.',inputFileExt]));
  if isempty(fileList)
    error('No files found in %s!',fullfile(dirs.dataroot,dirs.dataDir));
  else
    fprintf('Processing %d file(s) found in %s\n',length(fileList),fullfile(dirs.dataroot,dirs.dataDir));
  end
end

% Output file extension
if strcmp(outputFormat,'egi_egis')
  outputFileExt = 'egis';
elseif strcmp(outputFormat,'egi_sbin')
  outputFileExt = 'raw';
elseif strcmp(outputFormat,'ep_mat')
  outputFileExt = 'ept';
else
  error('Output filetype %s is unknown.',outputFormat);
end

% set the output directory
dirs.outputDir = fullfile(dirs.dataroot,dirs.baseDir,sprintf('3_ep_%s_out',outputFileExt));
if ~exist(dirs.outputDir,'dir')
  mkdir(dirs.outputDir);
end

% note the originating directory so we can go back to it
dirs.currDir = pwd;

% go through the file list, read it in and save it out
for i = 1:length(fileList)
  % get the file name
  [filePath,fileName,fileExt] = fileparts(fileList(i).name);
  
  inArgs = {
    'file',[fileName,fileExt],...
    'format',inputFormat,...
    };
  
  % cd to the input directory
  cd(fullfile(dirs.dataroot,dirs.dataDir));
  
  % read in the file
  fprintf('Reading %s...\n',fullfile(dirs.dataroot,dirs.dataDir,[fileName,'.',inputFileExt]));
  [EPdata] = ep_readData(inArgs);
  fprintf('Done.\n');
  
  fprintf('You should check on the Artifact_Correction_Log for this subject regarding globally bad channels, shorted channels, and bad subject information.\n');
  
  % In order to store data about bad channels, we might want to also have
  % the 'Artifact_Correction_Log fileName date.txt' files in this
  % directory.
  
  % Things we need to be aware of:
  % 1. Globally bad channels that were not interpolated
  % 2. Shorted channels
  % 3. Bad subjects
  
  % If a globally bad channel was successfully interpolated, is it still
  % listed under "Global bad Channels"? This would mean that any channel
  % listed under the "Global bad Channels" either was not or could not be
  % interpolated.
  
  % What to do about shorted channels? (Channels that are perfectly or
  % highly correlated)
  
  % What to do about bad subjects?
  
  
  % cd to the output directory
  cd(dirs.outputDir);
  
  if ~exist(fullfile(dirs.outputDir,[fileName,outputFileExt]),'file')
    % write out the file
    fprintf('Saving %s...\n',fullfile(dirs.dataroot,dirs.outputDir,[fileName,'.',outputFileExt]));
    ep_writeData(EPdata,fileName,outputFormat);
    fprintf('Done.\n');
  else
    error('%s already exists.',fullfile(dirs.outputDir,[fileName,'.',outputFileExt]));
    %fprintf('%s already exists. Moving on.\n',fullfile(dirs.outputDir,[fileName,outputFileExt]));
    %continue
  end
  
end

% go back to the original directory
cd(dirs.currDir);

%end
