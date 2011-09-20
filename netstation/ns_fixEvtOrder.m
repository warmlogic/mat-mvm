function ns_fixEvtOrder(newExt,dataroot,prepost)
%NS_FIXEVTORDER
%
% ns_fixEvtOrder(newExt,dataroot,prepost)
%
% Inputs:
%
%  newExt:      the extension of the new evt file (default: '_e.evt', to
%               match the output file extension from EP Toolkit)
%
%  dataroot:    the location of the original .evt files; the new .evt files
%               are saved here.
%               (default: pwd), meaning that by default, all .evt files in
%               the current working dir are processed
%
%  prepost:     the segmentation time (in ms) before and after each event.
%               e.g., [1000 2000] is 1000ms before and 2000ms after event.
%               This only applies when events were exported using relative
%               time mode, for which the actual event onset time is
%               recorded. Onset time is not needed for segmented files, so
%               if you just use Epoch Time in the event export tool you can
%               ignore this setting. If using relative time, all events
%               must have the same segmentation time.
%               (default: []; will error if using a relative time file)
%
% IMPORTANT: Make sure you don't have the original evt file open (e.g., in
%            a spreadsheet app) or Matlab won't be able to read it (the
%            error will say: 'Invalid file identifier.')
%
% This script sorts Net Station evt files (exported from any segmented file
% using the Event Export tool) by category and then orders the events by
% occurrence time.  This is useful for when NS incorrectly labels the time
% in evt file created from the seg file such that the onset time refers to
% the event's occurence in the continuous session file.  Instead, the onset
% time should refer to the event's location in the segmented file, where
% the events are separate in time but are contiguous in the file.  This
% script corrects the onset time to the latter situation sorted first by
% category (alphabetically) and then by onset time.
%
% To create the evt file(s) that this script reads:
% Make an Event Export tool in NS; check the "Export Category Name" box and
% set Time Mode to Epoch Time (Relative Time also works, but be sure to set
% the prepost variable correctly).  It does not matter whether you sort by
% category or by onset.  Of course, you also need to choose the event types
% to export. Drag the segmented NS file created from the original
% continuous recording file onto the tool to find the event types that you
% want.
%
% The new evt file can be used to markup a NS file created from the
% SBIN/RAW or EGIS file created by EP Toolkit (this is my only reason for
% creating this script). Use the Markup From File tool in NS; note that the
% NS and evt files must have the filenames.
%

% By Matt Mollison (2010-2011), matt.mollison@gmail.com

%% Can be set by the user

if nargin < 3
  % where the evt files are located
  dataroot = pwd;
  if nargin < 2
    % extension for the new .evt file
    newExt = '_e.evt';
    if nargin < 1
      % pre- and post-event segmentation times
      prepost = [];
    end
  end
end

%% Read in the old evt and write out the new evt

% Find the evt files
fprintf('Searching for .evt files in %s...',dataroot);
evt = dir(fullfile(dataroot,'*.evt'));
if ~isempty(evt)
  fprintf('found %d.\n',length(evt));
else
  error('No .evt files found!');
end

for i = 1:length(evt)
  fprintf('Processing %s...',evt(i).name);
  
  % initialize the timeMode
  timeModeStr = [];
  
  %% Read in the evt file
  [path,name,ext] = fileparts(evt(i).name);
  
  % figure out how many columns there are
  fid = fopen(fullfile(dataroot,evt(i).name),'r');
  % get the file name line
  fgetl(fid);
  % get the time mode line
  fgetl(fid);
  % get the header line
  tline = fgetl(fid);
  % headers
  if strcmp(tline(end),sprintf('\t'))
    tline = tline(1:end-1);
  end
  hdr = regexp(tline,'\t','split');
  cols.onset = find(strcmp(hdr,'Onset'));
  cols.categ = find(strcmp(hdr,'Category'));
  
  % get the first data line
  tline = fgetl(fid);
  % close the file
  fclose(fid);
  % if the last character is a tab, remove it
  if strcmp(tline(end),sprintf('\t'))
    tline = tline(1:end-1);
  end
  % since it's tab delimited, split it and count the length
  numCols = length(regexp(tline,'\t','split'));
  
  fid = fopen(fullfile(dataroot,evt(i).name),'r');
  data = textscan(fid,repmat('%s',1,numCols),'Headerlines',3,'Delimiter','\t');
  fclose(fid);
  
  %% 1. sort by category
  if ~issorted(data{cols.categ})
    [categ,categInd] = sort(data{cols.categ});
    for j = 1:numCols
      data{j} = data{j}(categInd);
    end
  end
  
  %% 2. sort by onset time or epoch number within each category
  categs = unique(data{cols.categ});
  for j = 1:length(categs)
    thisCatInd = [];
    for k = 1:length(data{cols.categ})
      if strcmp(data{cols.categ}(k),categs{j})
        thisCatInd = [thisCatInd; k];
      end
    end
    
    thisCatOnsets = data{cols.onset}(thisCatInd);
    
    % if the onset times have epoch numbers inside square brackers, then
    % we're in epoch time mode
    if ~isempty(strfind(thisCatOnsets{1},'['))
      if isempty(timeModeStr)
        timeModeStr = 'Time Mode: Epoch Time';
      end
      
      % order by epoch numbers
      epochNumbers = [];
      for k = 1:length(thisCatOnsets)
        thisEpochNum = str2double(thisCatOnsets{k}(strfind(thisCatOnsets{k},'[')+1:strfind(thisCatOnsets{k},']')-1));
        epochNumbers = [epochNumbers; thisEpochNum];
      end
      [epoch,epochInd] = sort(epochNumbers);
      
      % sort the events within this category
      for k = 1:numCols
        data{k}(thisCatInd) = data{k}(epochInd + min(thisCatInd) - 1);
      end
      
    else
      % otherwise the onset times are _hh:mm:ss.sss; order by these times
      if isempty(timeModeStr)
        timeModeStr = 'Time Mode: Relative Time';
      end
      
      onsetTimesMS = [];
      for k = 1:length(thisCatOnsets)
        hrs = str2double(thisCatOnsets{k}(2:3));
        mins = str2double(thisCatOnsets{k}(5:6));
        secs = str2num(thisCatOnsets{k}(8:end));
        
        % convert to ms so we can sort
        h_ms = hrs * (1000 * 60 * 60);
        m_ms = mins * (1000 * 60);
        s_ms = secs * 1000;
        ms = h_ms + m_ms + s_ms;
        
        onsetTimesMS = [onsetTimesMS; ms];
      end
      [onsetMS,onsetMSInd] = sort(onsetTimesMS);
      
      % sort the events within this category
      for k = 1:numCols
        data{k}(thisCatInd) = data{k}(onsetMSInd + min(thisCatInd) - 1);
      end
      
    end % epoch/onset
  end % categs
  
  %% renumber epochs or onset times
  if ~isempty(strfind(data{cols.onset}{1},'['))
    % Epoch Time mode
    for k = 1:length(data{cols.onset})
      left = data{cols.onset}{k}(1:strfind(data{cols.onset}{k},'['));
      right = data{cols.onset}{k}(strfind(data{cols.onset}{k},']'):end);
      data{cols.onset}{k} = sprintf('%s%s%s',left,num2str(k),right);
    end
  else
    % Relative Time mode
    if isempty(prepost)
      error('Need to set the prepost variable. See HELP NS_FIXEVTORDER for more information.');
    end
    
    % onset times are in NS's _hh:mm:ss.sss format
    mstime = prepost(1);
    for k = 1:length(data{cols.onset})
      hrs = sprintf('%02.0f',floor(mstime / (1000 * 60 * 60)));
      mins = sprintf('%02.0f',floor(mod(mstime,(1000 * 60 * 60)) / (1000 * 60)));
      secs = sprintf('%06.3f',mod(mod(mstime,(1000 * 60 * 60)),(1000 * 60)) / 1000);
      timeStr = sprintf('_%s:%s:%s',hrs,mins,secs);
      
      data{cols.onset}{k} = sprintf('%s',timeStr);
      
      mstime = mstime + prepost(1) + prepost(2);
    end
  end
  
  %% save out the new evt file
  newEvtFile = [name,newExt];
  fprintf('Saving %s...',newEvtFile);
  
  outfile = fopen(fullfile(dataroot,newEvtFile),'wt');
  
  % put in the basic info at the top
  fprintf(outfile,'%s\r',name);
  fprintf(outfile,'%s\r',timeModeStr);
  fprintf(outfile,'Code\tLabel\tType\tTrack\tOnset\tDuration\tCategory\r');
  % and add the rest of the data
  for j = 1:length(data{cols.categ})
    thisline = sprintf('data{1}{%d}',j);
    for k = 2:numCols
      thisline = [thisline sprintf(',data{%d}{%d}',k,j)];
    end
    formatStr = repmat('%s\t',1,numCols);
    eval(sprintf('fprintf(outfile,''%s%s'',%s);',formatStr(1:end-2),'\r',thisline));
  end
  fclose(outfile);
  fprintf('Done.\n');
  
end
