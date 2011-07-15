function ns_rejectEventsBCI(dataroot,subject,session,badFilters,nChan,sortField)
%NS_REJECTEVENTSBCI: Reject events by writing out a new bci file
%
% ns_rejectEventsBCI(dataroot,subject,session,badFilters,nChan,sortField)
%
% badFilters.expr     = expression for the events to reject
% badFilters.varargin = optional args to be used in the expr
% nChan               = the number of channels used in this recording
%                       (default: 129)
% sortField           = the field name by which to sort the events; this
%                       needs to be the same as the NS categories in the
%                       bci file; ns_addArtifactInfo will add this field,
%                       though then you have to run it again after this fxn
%
% EXAMPLE:
%   badFilters.expr     = 'rt > 1000 & ismember(type,varargin{1})'; % first
%                          argument of ismember is event struct field name
%   badFilters.varargin = {'TEST_TARGET'}; % cell array of event types
%
% See the eeg_toolbox's filterStruct.m function for more information
%
% NB: This does not modify the events.mat file! You need to re-reun
% NS_ADDARTIFACTINFO to do that.
%
% See also: NS_ADDARTIFACTINFO, FILTERSTRUCT

if nargin < 6
  sortField = 'nsCategory';
  if nargin < 5
    nChan = 129;
  end
end

if ~isfield(badFilters,'varargin')
  badFilters.varargin = {};
end

% read everything in as strings so it's easier to write it back out
format_str = ['%s%s%s%s',repmat('%s',[1,nChan*2]),'%s'];
%format_str = ['%s%d8%d8%s',repmat('%d8',[1,nChan*2]),'%s'];

nCols = 5 + (nChan*2);

% % debug
% dataroot = '/Volumes/curranlab/Data/COSI/eeg/eppp/-1000_2000';
% subject = 'COSI001';
% session = 'session_0';

% assuming behavioral events are stored in:
% experiment/eeg/behavioral/subject/session/events/events.mat
behroot = fullfile(dataroot(1:strfind(dataroot,'eeg')+2),'behavioral');

% assumes bci files are stored in dataroot/session/ns_bci
bcipath = fullfile(dataroot,session,'ns_bci');

backupDir = 'orig';
if ~exist(fullfile(bcipath,backupDir),'dir')
  mkdir(fullfile(bcipath,backupDir));
end

% load in the newest PyEPL events
eventsDir = fullfile(behroot,subject,session,'events');
fprintf('Loading events for %s, %s...\n',subject,session);
if exist(fullfile(eventsDir,'events.mat'),'file')
  events = loadEvents(fullfile(eventsDir,'events.mat'));
  fprintf('Done.\n');
else
  fprintf('events.mat for %s %s does not exist. Moving on.\n',subject,session);
  return
end

fprintf('Loading bci file for %s, %s...\n',subject,session);

% define the metadata NS export file with the session summary
summaryFile = dir(fullfile(bcipath,[subject,'*.bci']));

% make sure we found only one bci file (e.g., multiple sessions)
if length(summaryFile) > 1
  if isfield(events,'nsFile')
    summaryFile = dir(fullfile(bcipath,[events(1).nsFile,'.bci']));
  else
    error('More than one bci file, and no nsFile field to denote which to choose: %s',fullfile(bcipath,[subject,'*.bci']));
  end
end

if ~isempty(summaryFile)
  % read in session summary file
  fid = fopen(fullfile(bcipath,summaryFile.name),'r');
  sesSummary = textscan(fid,format_str,'Headerlines',1,'delimiter','\t');
  fclose(fid);
  
  % copy the bci file to the old file
  unix(sprintf('cp ''%s'' ''%s''',fullfile(bcipath,summaryFile.name),fullfile(bcipath,backupDir,sprintf('%s.old',summaryFile.name))));
else
  warning([mfilename,':no_bci_file'],'MISSING FILE: %s',fullfile(bcipath,[subject,'*.bci']));
  return
end

% find only the categories that we have
[events] = filterStruct(events,'ismember(nsCategory,varargin{1})',badFilters.eventValues);

% sort the events
if isfield(events,sortField)
  [~,evIndex] = sort({events.(sortField)});
  events = events(evIndex);
else
  error('Field %s does not exist.',sortField);
end

% find the bad events
[~,ind] = filterStruct(events,badFilters.expr,badFilters.varargin);

% mark them as bad
alreadyBadCounter = 0;
sesSummary{4}(ind) = {'bad'};
for i = 1:length(ind)
  % also check that the events and the bci entries are in the right order
  if ~strcmp(sesSummary{1}(i),events(i).nsCategory)
    error('The events are not in the same order as the bci file. Something is really wrong.');
  end
  
  if ind(i) == 1
    if strcmp(sesSummary{nCols}(i),'No Faults')
      sesSummary{nCols}(i) = {'manual'};
    else
      sesSummary{nCols}(i) = {sprintf('%s,%s',sesSummary{nCols}{i},'manual')};
      alreadyBadCounter = alreadyBadCounter + 1;
    end
  end
end
fprintf('%d bad events marked based on these criteria: %s\n',sum(ind),badFilters.expr);
fprintf('%d would have already been bad due to artifacts.\n',alreadyBadCounter);
fprintf('kept %d out of %d events.\n',sum(~ind),length(ind));

% start the new bci file
bcifilenew = fullfile(bcipath,sprintf('%s',summaryFile.name));
outfile = fopen(bcifilenew,'wt');

% print the header
fprintf(outfile,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Category','StartTime','Seg#','Status','ChannelInfo(1 to nChannels)','. . . ','Faults');

% write out each line of the bci file
for i = 1:length(sesSummary{1})
  lineStr = [];
  for j = 1:nCols
    if ~isempty(cell2mat(sesSummary{j}(i)))
      lineStr = sprintf('%s%s',lineStr,cell2mat(sesSummary{j}(i)));
      if j < nCols
        if ~isempty(cell2mat(sesSummary{j+1}(i)))
          lineStr = sprintf('%s\t',lineStr);
        end
      end
    end
  end
  
  fprintf(outfile,'%s\n',lineStr);
end % sesSummary

fclose(outfile);
fprintf('Saved new bci file, rejecting events based on these criteria: %s\n',badFilters.expr);

end
