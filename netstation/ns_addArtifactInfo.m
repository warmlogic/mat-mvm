function [events,goodEv] = ns_addArtifactInfo(dataroot,subject,session,evFilters,nChan,overwriteArtFields)
%NS_ADDARTIFACTINFO: Add NS artifact information to a PyEPL event structure
%for doing behavioral analyses on the artifact-free subset of events
%
% [events,goodEv] = ns_addArtifactInfo(dataroot,subject,session,evFilters,nChan,overwriteArtFields)
%
% Input:
% nChan               = the number of channels used in this recording
%                       (default: 129)
%
% TODO: Documentation incomplete
%
%
% Expects rereferenced data with nChan channels (in the bci file)
%
% You must first create the bci file from the pare average rereferenced
% files using the export metadata NS function. 'segment information'
% exports the bci file.
%
% dataroot is: /path/to/experiment/eeg/processingType/time/ (I use
% something like that; anything else can be used, it just needs to contain
% the session folders)
%
% This function assumes bci files are stored in: dataroot/session/ns_bci/
%
% This function assumes that subject events are stored in:
% /path/to/experiment/eeg/behavioral/subject/session/events/events.mat
% (gets automatically created from dataroot)
%
% TODO: Can I also bring the more-detailed information into FT?  Maybe
% create an event struct with events sorted alphabetically by category (and
% within that sort by time) of only the good events.
%
% See inline code for some better explanations of what this function does

if nargin < 6
  overwriteArtFields = 0;
  if nargin < 5
    nChan = 129;
  end
end

format_str = ['%s%d8%d8%s',repmat('%d8',[1,nChan*2]),'%s'];

% % debug
% dataroot = '/Volumes/curranlab/Data/COSI/eeg/eppp/-1000_2000';
% subject = 'COSI001';
% session = 'session_0';

% assuming behavioral events are stored in:
% experiment/eeg/behavioral/subject/session/events/events.mat
behroot = fullfile(dataroot(1:strfind(dataroot,'eeg')+2),'behavioral');

% assumes bci files are stored in dataroot/session/ns_bci
bcipath = fullfile(dataroot,session,'ns_bci');

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

fprintf('Getting NS artifact info for %s, %s...\n',subject,session);

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
else
  warning([mfilename,':no_bci_file'],'MISSING FILE: %s',fullfile(bcipath,[subject,'*.bci']));
  return
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
eventValues = evFilters.eventValues;

% determine whether there are extra fields are so we can include them in
% filtering; this relies on the assumption that there are always fields for
% 'type' and 'filters'; currently any extra fields must be strings
%
% The purpose of extraFilters is to provide additional filters. The name of
% the extraFilter fild and its value (a string) correspond to this info
% in a struct. For example, if you wanted to separate two conditions, color
% and side, and the field in events.mat is 'cond', the extraFilter field
% name is cond and its value is 'color' for one condition, and cond='side'
% for the other.
fn = fieldnames(evFilters.(eventValues{1}));
if sum(~ismember(fieldnames(evFilters.(eventValues{1})),'filters') & ~ismember(fieldnames(evFilters.(eventValues{1})),'type')) > 0
  extraFilters = fn(~ismember(fieldnames(evFilters.(eventValues{1})),'filters') & ~ismember(fieldnames(evFilters.(eventValues{1})),'type'));
else
  extraFilters = [];
end

% separate the NS event categories
for evVal = 1:length(eventValues)
  nsEvents.(eventValues{evVal}) = filterStruct(artEv_ns,'ismember(category,varargin{1})',eventValues{evVal});
end

% exit out if we don't want to overwriteArtFields
if isfield(events,'nsArt')
  if overwriteArtFields == 0
    fprintf('NS artifact information has already been added to this struct. Moving on.\n');
    return
  elseif overwriteArtFields == 1
    fprintf('Removing NS artifact information from this struct; will overwriteArtFields them with current information.\n');
    events = rmfield(events,{'nsArt','nsBadChan','nsBadReason'});
    if isfield(events,'nsCategory')
      events = rmfield(events,'nsCategory');
    end
  end
else
  fprintf('No NS artifact information exists. Adding it.\n');
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
    evValStr = sprintf('strcmp(''%s'',''%s'')',events(i).type,evFilters.(eventValues{evVal}).type);
    % add any extra fields if necessary
    if ~isempty(extraFilters)
      for ef = 1:size(extraFilters,1)
        evValStr = sprintf('%s && strcmp(''%s'',''%s'')',evValStr,events(i).(extraFilters{ef}),evFilters.(eventValues{evVal}).(extraFilters{ef}));
      end
    end
    
    % see if it matches this eventValue (and any extra fields)
    if eval(evValStr)
      % construct the filter requirements statement
      filtStr = sprintf('events(i).%s',evFilters.(eventValues{evVal}).filters{1});
      % continue constructing
      if length(evFilters.(eventValues{evVal}).filters) > 1
        for filt = 2:length(evFilters.(eventValues{evVal}).filters)
          filtStr = cat(2,filtStr,sprintf(' && events(i).%s',evFilters.(eventValues{evVal}).filters{filt}));
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
        
        % add the NS Category
        events(i).nsCategory = eventValues{evVal};
        
      end % if matches filter requirements
    end % if matches event type
  end % for each event type
  
  % if this doesn't match any of the filter requirements for any of the
  % events, then we don't want to see this event
  if matchedFilter == 0
    events(i).nsArt = -1;
    events(i).nsBadChan = [];
    events(i).nsBadReason = '';
    events(i).nsCategory = '';
  end
end

% badcEv = filterStruct(events,'nsArt == 1 & ismember(nsBadReason,varargin{1})','badc');
% badc = unique(getStructField(badcEv,'nsBadChan'));
% 
% eyemeyebEv = filterStruct(events,'nsArt == 1 & ismember(nsBadReason,varargin{1})','eyem,eyeb');
% eyemEv = filterStruct(events,'nsArt == 1 & ismember(nsBadReason,varargin{1})','eyem');
% eyebEv = filterStruct(events,'nsArt == 1 & ismember(nsBadReason,varargin{1})','eyeb');

% grab only the good events
goodEv = {};
for evVal = 1:length(eventValues)
  % construct the filter
  filtStr = [sprintf('nsArt == 0 & ismember(type,varargin{1})'),sprintf(repmat(' & %s',1,length(evFilters.(eventValues{evVal}).filters)),evFilters.(eventValues{evVal}).filters{:})];
  varargStr = sprintf('{''%s''}',evFilters.(eventValues{evVal}).type);
  if ~isempty(extraFilters)
    for ef = 1:size(extraFilters,1)
      filtStr = sprintf('%s & ismember(%s,varargin{%d})',filtStr,extraFilters{ef},ef+1);
      if ischar(evFilters.(eventValues{evVal}).(extraFilters{ef}))
        varargStr = sprintf('%s,{''%s''}',varargStr,evFilters.(eventValues{evVal}).(extraFilters{ef}));
      elseif isnumeric(evFilters.(eventValues{evVal}).(extraFilters{ef}))
        error('The extraFilter %s (set to %s) is numeric. You must add it to the regular filter list.',extraFilters{ef},num2str(evFilters.(eventValues{evVal}).(extraFilters{ef})));
      else
        error('Do not know what to do with input for this extraFilter.');
      end
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
    
    unix(sprintf('cp ''%s'' ''%s''',fullfile(eventsDir,oldEv(oe).name),backupfile));
    
    %saveEvents(events,backupfile);
  end
  fprintf('Done.\n');
end

% save the events
fprintf('Saving events with NS artifact information...\n');
saveEvents(events,fullfile(eventsDir,'events.mat'));
fprintf('Done.\n');

% % save the good events
% fprintf('Saving only good events...\n');
% saveEvents(goodEv,fullfile(eventsDir,'events_good.mat'));

fprintf('Artifact addition complete.\n');
