function ns_fixEvtOrder(numCols,dataroot)
%NS_FIXEVTORDER
%
% ns_fixEvtOrder(numCols,dataroot)
%
% no inputs necessary. processes all .evt files in the current working dir
% by default. otherwise, specify where the .evt files are located.
%
% numCols is the number of columns in each original .evt file
%
% This script sorts a NetStation evt file (exported from a segmented file)
% by category and then orders the events by occurrence time. This is useful
% for when NS incorrectly labels the time in evt file created from the seg
% file such that the onset time refers to the event's occurence in the
% continuous session file. Instead, the onset time should refer to the
% event's location in the segmented file, where the events are separate in
% time but are contiguous in the file. This script corrects the onset time
% to the latter situation sorted first by category (alphabetically) and
% then by onset time.
%
% To create the evt file that this script reads in, make an Event Export
% tool in NS; check the "Export Category Name" box and set Time Mode to
% either Relative or Epoch Time. It does not matter whether you sort by
% category or by onset. Of course, you also need to choose the event types
% to export. Drag the segmented NS file created from the original
% continuous recording file onto the tool to find the event types that you
% want.
%
% The new evt file can be used to markup a nsf file created from the sbin
% or egis file created by EP_Toolkit (this is my only reason for creating
% this script).
%
% In this script the user must set:
%  1) the number of columns in the original evt file
%  2) where the original evt files (created from the seg files) are located
%     (default is the current directory)
%  3) the column numbers in the evt file for Category and Onset
%  4) the segmentation times before and after an event (all events must
%     have the same length)
%  5) the extension for the new evt file
%
%  NB:
%  1 and 2 can be set in the function arguments
%  3, 4, and 5 are hardcoded in the script
%

%% To be set by the user

% 1) where the evt files are located
if nargin < 2
  dataroot = pwd;
  if nargin < 1
    %numCols = 55;
    numCols = 43;
  end
end

% 2) the number of columns in the original evt file

% 3) the columns for Onset and Category
cols.onset = 5;
cols.categ = 7;

% 4) the amount of ms segmented before and after the event, respectively
prepost = [1000 2000];

% 5) the extension to be used on the new evt file
evtExt = '_e.evt';

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
  
  fid = fopen(fullfile(dataroot,evt(i).name),'r');
  data = textscan(fid,repmat('%s',1,numCols),'Headerlines',3,'Delimiter','\t');
  fclose(fid);
  
  %% 1. sort by category
  if ~issorted(data{cols.categ})
    [categ categInd] = sort(data{cols.categ});
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
      [epoch epochInd] = sort(epochNumbers);
      
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
      [onsetMS onsetMSInd] = sort(onsetTimesMS);
      
      % sort the events within this category
      for k = 1:numCols
        data{k}(thisCatInd) = data{k}(onsetMSInd + min(thisCatInd) - 1);
      end
      
    end % epoch/onset
  end % categs
  
  %% renumber epochs or onset times
  if ~isempty(strfind(data{cols.onset}{1},'['))
    % epochs
    for k = 1:length(data{cols.onset})
      left = data{cols.onset}{k}(1:strfind(data{cols.onset}{k},'['));
      right = data{cols.onset}{k}(strfind(data{cols.onset}{k},']'):end);
      data{cols.onset}{k} = sprintf('%s%s%s',left,num2str(k),right);
    end
  else
    % onset times in NS's _hh:mm:ss.sss format
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
  newEvtFile = [name,evtExt];
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
