function [events,goodEv] = ns_addArtifactInfo(dataroot,subject,session,nsEvFilters,overwrite)
%NS_ADDARTIFACTINFO Add NS artifact information to a PyEPL event structure
%for doing behavioral analyses on the artifact-free subset of events
%
% [events,goodEv] = ns_addArtifactInfo(dataroot,subject,session,nsEvFilters,overwrite)
%
% Expects rereferenced data with 129 channels
%
% only processes these events: RCR, RHSC, RHSI
%
% create the bci file from the pare average rereferenced files (step 7)
% using the export metadata NS function. 'segment information' exports the
% bci file.
%
% Can I also bring the more-detailed information into FT?  Maybe create an
% event struct with events sorted alphabetically by category (and within
% that sort by time) of only the good events.

if nargin < 5
  overwrite = 0;
end

nChan = 129;
format_str = ['%s%d8%d8%s',repmat('%d8',[1,nChan*2]),'%s'];
%format_str = ['%s%s%s%s',repmat('%s',[1,nChan*2]),'%s'];

% % debug
% dataroot = '/Volumes/curranlab/Data/RSRC/eeg/eppp_nobackup/-1000_2000';
% % debug
% subject = 'RSRC021';

% % debug
% dataroot = '/Volumes/RAID/curranlab/Data/SLORK/eeg/eppp_nobackup/-1000_2000';
% % debug
% subject = 'SLORK002';
% moreFilter = {'ismember(testType,varargin{1}),{''side''}','ismember(testType,varargin{1}),{''task''}'};

% assuming events are stored in
% eperiment/eeg/behavioral/subject/session/events/events.mat
behroot = fullfile(dataroot(1:strfind(dataroot,'eeg')+2),'behavioral');

fprintf('Getting NS artifact info for %s, %s...\n',subject,session);

% define the metadata NS export file with the session summary
summaryFile = dir(fullfile(dataroot,'ns_bci',[subject,'*.bci']));

if ~isempty(summaryFile)
  summaryFile = fullfile(dataroot,'ns_bci',summaryFile.name);
  % read in session summary file
  fid = fopen(summaryFile,'r');
  sesSummary = textscan(fid,format_str,'Headerlines',1,'delimiter','\t');
  fclose(fid);
else
  error('MISSING FILE: %s',fullfile(dataroot,'ns_bci',[subject,'*.bci']));
end

% create event struct for session summary
artEv_ns = struct('category',sesSummary{1},'status',sesSummary{4},'badChan',[],'reason',sesSummary{length(sesSummary)});
for ev = 1:length(artEv_ns)
  % if we find a bad event
  if strcmp(artEv_ns(ev).status,'bad')
    % find the channels that were bad
    for c = 1:nChan
      if sesSummary{5+(c*2-1)}(ev) == 0
        artEv_ns(ev).badChan = [artEv_ns(ev).badChan c];
      end
    end
  end
end

% make it so the event value cell array is easier to access
eventValues = nsEvFilters.eventValues;

% figure if there are extra fields are so we can include them in filtering;
% this relies on the assumption that there are always fields for 'type' and
% 'filters'; currently any extra fields must be strings
fn = fieldnames(nsEvFilters.(eventValues{1}));
if sum(~ismember(fieldnames(nsEvFilters.SCR),'filters') & ~ismember(fieldnames(nsEvFilters.SCR),'type')) > 0
  extraFields = fn(~ismember(fieldnames(nsEvFilters.SCR),'filters') & ~ismember(fieldnames(nsEvFilters.SCR),'type'));
else
  extraFields = [];
end

% separate the NS event categories
for evVal = 1:length(eventValues)
  nsEvents.(eventValues{evVal}) = filterStruct(artEv_ns,'ismember(category,varargin{1})',eventValues{evVal});
end

% load in the newest PyEPL events
eventsDir = fullfile(behroot,subject,session,'events');
fprintf('Loading events...');
events = loadEvents(fullfile(eventsDir,'events.mat'));
fprintf('Done.\n');

% exit out if we don't want to overwrite
if overwrite == 0
  if isfield(events,'nsArt')
    fprintf('NS artifact information has already been added to this struct. Moving on.\n');
    return
  end
elseif overwrite == 1
  if isfield(events,'nsArt')
    fprintf('Removing NS artifact information from this struct.\n');
    events = rmfield(events,{'nsArt','nsBadChan','nsBadReason'});
  end
end

% initialize the count of the number of NS events we've gone through
for evVal = 1:length(eventValues)
  countAll.(eventValues{evVal}) = 0;
end

for i = 1:length(events)
  % initialize to see if we match the filter in one of the event types
  matchedFilter = 0;
  
  for evVal = 1:length(eventValues)
    % start the string to evaluate using the type field, because that will
    % always be there
    evValStr = sprintf('strcmp(''%s'',''%s'')',events(i).type,nsEvFilters.(eventValues{evVal}).type);
    % add any extra fields if necessary
    if ~isempty(extraFields)
      for ef = 1:size(extraFields,1)
        evValStr = sprintf('%s && strcmp(''%s'',''%s'')',evValStr,events(i).(extraFields{ef}),nsEvFilters.(eventValues{evVal}).(extraFields{ef}));
      end
    end
    
    % see if it matches this eventValue (and any extra fields)
    if eval(evValStr)
      % construct the filter requirements statement
      filtStr = sprintf('events(i).%s',nsEvFilters.(eventValues{evVal}).filters{1});
      % continue constructing
      if length(nsEvFilters.(eventValues{evVal}).filters) > 1
        for filt = 2:length(nsEvFilters.(eventValues{evVal}).filters)
          filtStr = cat(2,filtStr,sprintf(' && events(i).%s',nsEvFilters.(eventValues{evVal}).filters{filt}));
        end
      end
      % then see if it matches the filter requirements
      if eval(filtStr)
        % mark that we matched the filter
        matchedFilter = 1;
        % increase the count for this eventValue
        countAll.(eventValues{evVal}) = countAll.(eventValues{evVal}) + 1;
        if strcmp(nsEvents.(eventValues{evVal})(countAll.(eventValues{evVal})).status,'bad')
          % mark it as an artifact
          events(i).nsArt = 1;
          events(i).nsBadChan = nsEvents.(eventValues{evVal})(countAll.(eventValues{evVal})).badChan;
          events(i).nsBadReason = nsEvents.(eventValues{evVal})(countAll.(eventValues{evVal})).reason;
        else
          % mark it as good
          events(i).nsArt = 0;
          events(i).nsBadChan = [];
          events(i).nsBadReason = '';
        end
      end % if matches filter requirements
    end % if matches event type
  end % for each event type
  
  % if this doesn't match any of the filter requirements for any of the
  % events, then we don't want to see this event
  if matchedFilter == 0
    events(i).nsArt = -1;
    events(i).nsBadChan = [];
    events(i).nsBadReason = '';
  end
end

%badcEv = filterStruct(events,'ismember(nsBadReason,varargin{1})','badc');
%badc = unique(getStructField(badcEv,'nsBadChan'));

%eyemeyebEv = filterStruct(events,'ismember(nsBadReason,varargin{1})','eyem,eyeb');
%eyemEv = filterStruct(events,'ismember(nsBadReason,varargin{1})','eyem');
%eyebEv = filterStruct(events,'ismember(nsBadReason,varargin{1})','eyeb');

% grab only the good events
goodEv = {};
for evVal = 1:length(eventValues)
  % construct the filter
  filtStr = [sprintf('nsArt == 0 & ismember(type,varargin{1})'),sprintf(repmat(' & %s',1,length(nsEvFilters.(eventValues{evVal}).filters)),nsEvFilters.(eventValues{evVal}).filters{:})];
  varargStr = sprintf('{''%s''}',nsEvFilters.(eventValues{evVal}).type);
  if ~isempty(extraFields)
    for ef = 1:size(extraFields,1)
      filtStr = sprintf('%s & ismember(%s,varargin{%d})',filtStr,extraFields{ef},ef+1);
      varargStr = sprintf('%s,{''%s''}',varargStr,nsEvFilters.(eventValues{evVal}).(extraFields{ef}));
    end
  end
  % filter the struct so we only get good events
  goodEvAll.(eventValues{evVal}) = eval(sprintf('filterStruct(events,''%s'',%s)',filtStr,varargStr));
  % concatenate good event structs together
  goodEv = cat(2,goodEv,struct2cell(goodEvAll.(eventValues{evVal})));
end
goodEv = cell2struct(goodEv,fieldnames(events),1);

% always save backup events
oldEv = flipud(dir(fullfile(eventsDir,'*events*mat*')));
if ~isempty(oldEv)
  for oe = 1:length(oldEv)
    backupfile = [fullfile(eventsDir,oldEv(oe).name),'.old'];
    fprintf('Making backup: copying %s to %s\n',fullfile(eventsDir,oldEv(oe).name),backupfile);
    
    unix(sprintf('cp %s %s',fullfile(eventsDir,oldEv(oe).name),backupfile));
    
    %saveEvents(events,backupfile);
  end
  fprintf('Done.\n');
end

% save the events
fprintf('Saving events with NS artifact information...\n');
saveEvents(events,fullfile(eventsDir,'events.mat'));
fprintf('Done.\n');
% save the good events
fprintf('Saving only good events...\n');
saveEvents(goodEv,fullfile(eventsDir,'events_good.mat'));
fprintf('Done.\n');
