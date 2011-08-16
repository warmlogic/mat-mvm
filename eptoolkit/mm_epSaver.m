%function mm_epSaver(experName,prepost,inputFormat,outputFormat)
%MM_EPSAVER: Use EP's read/write to save EEG data
%
% This is a simple input/output wrapper to, e.g., save ep_mat format to an
% EGI-readable format (egi_egis, egi_sbin)
%
% Input:
%
%   experName    = experiment name ('COSI')
%   prepost      = pre- and post-stimulus timing (in seconds: [-1 2])
%   inputFormat  = format to read ('ep_mat')
%   outputFormat = format to save ('egi_egis')
%
%   Original files should be located in:
%     dataroot/experName/eeg/eppp/pre_post_ms/3_ep_[inputFormatExtension]
%     e.g.,
%     dataroot/COSI/eeg/eppp/-1000_2000/3_ep_mat
%
% Output:
%
%   Files are saved to 
%     dataroot/experName/eeg/eppp/pre_post_ms/3_ep_[outputFormatExtension]_out
%     e.g.,
%     dataroot/COSI/eeg/eppp/-1000_2000/3_ep_egis_out
%

experName = 'COSI';
prepost = [-1 2];
inputFormat = 'ep_mat';
outputFormat = 'egi_egis';

% Input file extension
if strcmp(inputFormat,'egi_egis')
  inputFileExt = 'egis';
elseif strcmp(inputFormat,'egi_sbin')
  inputFileExt = 'raw';
elseif strcmp(inputFormat,'ep_mat')
  inputFileExt = 'mat';
else
  error('Input filetype %s is unknown.',inputFormat);
end

% Location of the data files (dataroot)s
dirs.baseDir = fullfile('eeg','eppp',sprintf('%d_%d',prepost(1)*1000,prepost(2)*1000));
dirs.dataDir = fullfile(dirs.baseDir,sprintf('3_ep_%s',inputFileExt));

dirs.homeDir = getenv('HOME');

% 2 potential Curran server dataroots
dirs.serverDir = fullfile('/Volumes/curranlab/Data',experName);
dirs.serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data',experName);
dirs.localDir = fullfile(dirs.homeDir,'data',experName);
% Dream dataroot
dirs.dreamDir = fullfile('/data/projects/curranlab',experName);

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
  outputFileExt = 'mat';
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
  [~,fileName,fileExt] = fileparts(fileList(i).name);
  
  inArgs = {
    'file',[fileName,'.',fileExt],...
    'format',inputFormat,...
    };
  
  % cd to the input directory
  cd(fullfile(dirs.dataroot,dirs.dataDir));
  
  % read in the file
  [EPdata] = ep_readData(inArgs);
  
  % cd to the output directory
  cd(dirs.outputDir);
  
  % write out the file
  ep_writeData(EPdata,fileName,outputFormat);
  
end

% go back to the original directory
cd(dirs.currDir);

%end
