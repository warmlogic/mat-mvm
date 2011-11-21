function ns_average_rep(dataroot,inputTimesMS,inputFormat,timesToAvgMS,outputFormat,numberOfReps)
%ns_average_rep: Replicates average of specified samples in an EEG file 
%
% Usage:
%   ns_average_rep(dataroot,inputTimesMS,inputFormat,timesToAvgMS,outputFormat,numberOfReps);
%
% Requires ERP PCA Toolkit: http://sourceforge.net/projects/erppcatoolkit/
% (used for reading and writing Net Station files)
%
% Input files must be exported from Net Station as either EGIS (single
% subject single trials or average; extension = .egis) or SBIN/RAW
% (extension = .raw)
%
% Input:
%
%  Required:
%   dataroot     = full path of directory containing the data files
%                  (e.g., '/Volumes/curranlab/Data/PCOR/PCORf/NetStation Data/Short Epochs/ICA_Blink/avg/fMRIonly_ex06&29')
%   inputTimesMS = the segment length of the files. Used to index the
%                  selection of times to average.
%                  (e.g., [-200 1000])
%   inputFormat  = the format of the input files ('egi_egis' for single
%                  subject single trials, 'egi_egia' for average, or
%                  'egi_sbin' raw files)
%   timesToAvgMS = the times to average across (and repeat), in ms
%                  (e.g., [500 800])
%  Optional:
%   outputFormat = the format to write the files to (default = inputFormat)
%   numberOfReps = number of repetitions of the averaged data
%                  (default = the number of samples in timesToAvgMS)
%                  (minimum = 5 because it seems like NS's only requirement
%                  for reading the output EGIS file is that it is larger
%                  than 12KB)
%
% Output:
%
%   Saves files (format = outputFormat) in dataroot, with timesToAvgMS
%   added to the name. Any periods in the input file name will be converted
%   to underscores becasue of how EP Toolkit writes out files.
%   (e.g., if timesToAvgMS = [500 800], new file is origname_500_800.ext)
%
% NB: Before processing the output file further, you will need to either
% manually open the output file in NS to apply the correct electrode
% montage or use the Export to NS Format tool with "Missing SLF" defined
% correctly.
%

% make sure there are enough repetitions
if exist('numberOfReps','var')
  if numberOfReps < 5
    error('numberOfReps must be greater than 5 (this is an approximate guess; NS seems to need the output file(s) to be larger than 12KB).');
  end
end

% set the dataroot
dirs.dataroot = dataroot;

% Input file extension
if strcmp(inputFormat,'egi_egis')
  inputFileExt = 'egis';
elseif strcmp(inputFormat,'egi_egia')
  inputFileExt = 'egis';
elseif strcmp(inputFormat,'egi_sbin')
  inputFileExt = 'raw';
else
  error('Input filetype %s is unknown.',inputFormat);
end

% see which files are in the data directory
if ~exist(fullfile(dirs.dataroot),'dir')
  error('%s does not exist',fullfile(dirs.dataroot));
else
  fileList = dir(fullfile(dirs.dataroot,['*.',inputFileExt]));
  if isempty(fileList)
    error('No files found in %s!',fullfile(dirs.dataroot));
  else
    fprintf('Processing %d file(s) found in %s\n',length(fileList),fullfile(dirs.dataroot));
  end
end

% Output file extension
if ~exist('outputFormat','var')
  outputFormat = inputFormat;
end
if strcmp(outputFormat,'egi_egis')
  outputFileExt = 'egis';
elseif strcmp(outputFormat,'egi_egia')
  outputFileExt = 'egis';
elseif strcmp(outputFormat,'egi_sbin')
  outputFileExt = 'raw';
else
  error('Output filetype %s is unknown.',outputFormat);
end

% set the output directory
dirs.outputDir = dirs.dataroot;
% dirs.outputDir = fullfile(dirs.dataroot,sprintf('ave_%s_out',outputFileExt));
% if ~exist(dirs.outputDir,'dir')
%   mkdir(dirs.outputDir);
% end

% note the originating directory so we can go back to it
dirs.currDir = pwd;

% go through the file list, read it in, average, and save it out
for i = 1:length(fileList)
  % get the file name
  [filePath,fileName,fileExt] = fileparts(fileList(i).name);
  
  inArgs = {
    'file',[fileName,fileExt],...
    'format',inputFormat,...
    };
  
  % cd to the input directory
  cd(dirs.dataroot);
  
  % read in the file
  fprintf('Reading %s...\n',fullfile(dirs.dataroot,[fileName,'.',inputFileExt]));
  [EPdata] = ep_readData(inArgs);
  fprintf('Done.\n');
  
  % get the recording sample rate
  sampleRate = EPdata.Fs;
  
  % figure out which times we want
  inputTimeNamesMS = (inputTimesMS(1):(1000/sampleRate):(inputTimesMS(2)-1));
  timesToAvgIndex = logical(inputTimeNamesMS >= timesToAvgMS(1) & inputTimeNamesMS < timesToAvgMS(2));
  
  % specify the number of repetitions if it's not defined; also note the
  % time names to save out, linearly spaced between the first and last
  % times if using a custom number of repetitions
  if ~exist('numberOfReps','var')
    numberOfReps = sum(timesToAvgIndex);
    outputTimeNamesMS = EPdata.timeNames(timesToAvgIndex);
  else
    firstTime = find(inputTimeNamesMS >= timesToAvgMS(1));
    lastTime = find(inputTimeNamesMS < timesToAvgMS(2));
    outputTimeNamesMS = EPdata.timeNames(ceil(linspace(firstTime(1),lastTime(end),numberOfReps)));
  end
  
  % average the EEG data across specified timepoints and repeat it 
  fprintf('Averaging over %d--%d ms (%d samples)...',timesToAvgMS(1),timesToAvgMS(2),sum(timesToAvgIndex));
  rep_data = repmat(mean(EPdata.data(:,timesToAvgIndex,:),2),[1,numberOfReps,1]);
  fprintf('Done.\n');
  
  % put the average data back in the EPdata struct
  EPdata.data = rep_data;
  % set the timeNames to the new averaged subset (or something close)
  EPdata.timeNames = outputTimeNamesMS;
  
  % cd to the output directory
  cd(dirs.outputDir);
  
  % name the output file
  outputFileName = sprintf('%s_%d_%d',fileName,timesToAvgMS(1),timesToAvgMS(2));
  % remove any periods from the file name because EP doesn't like that
  outputFileName = strrep(outputFileName,'.','_');
  
  if ~exist(fullfile(dirs.outputDir,[outputFileName,outputFileExt]),'file')
    % write out the file
    fprintf('Saving %s...\n',fullfile(dirs.dataroot,[outputFileName,'.',outputFileExt]));
    ep_writeData(EPdata,outputFileName,outputFormat);
    fprintf('Done.\n');
  else
    error('%s already exists.',fullfile(dirs.outputDir,[outputFileName,'.',outputFileExt]));
    %fprintf('%s already exists. Moving on.\n',fullfile(dirs.outputDir,[fileName,outputFileExt]));
    %continue
  end
  
end

% go back to the original directory
cd(dirs.currDir);

end
