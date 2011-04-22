function rsrc3_prepData_events(subjects,session,prep_eeg)
% rsrc3_prepData_events(subjects,session,prep_eeg)
%
% Purpose
%   Create behavioral events; if prep_eeg: run prep_egi_data, export
%   NetStation events
%
% Inputs
%   subjects: a cell of subject number strings
%
%   session: 1 or 2; corresponds to pairs of {'session_0','session_1'} and
%   {'session_2','session_3'}
%
%   prep_eeg: vector of two numbers, 1 or 0; whether or not to prepare the
%   EEG data for the session pairs (default = [0 1])
%
% Outputs
%
%   dataroot is either:
%     /Volumes/curranlab/Data (first choice if server is mounted)
%     or
%     /Users/matt/data (fallback local storage)
%
%   Events (struct and NetStation) will be saved in the dataroot; either:
%     dataroot/RSRC3/eeg/behavioral/subject/session/events
%     or
%     dataroot/RSRC3/eeg/behavioral/subject/session/events
%
% Assumptions
%
%   The behavioral data is located in:
%     dataroot/RSRC3/eeg/behavioral/subject/session
%
%   The NetStation raw file is located in:
%     dataroot/RSRC3/eeg/behavioral/subject/session/eeg
%

%rsrc3_prepData_events({'RSRC3001','RSRC3002','RSRC3003','RSRC3004','RSRC3005','RSRC3006','RSRC3007','RSRC3008','RSRC3009','RSRC3010','RSRC3011','RSRC3012','RSRC3013'},1,[0 1]);
%rsrc3_prepData_summary(1,0);

expName = 'RSRC3';

serverDir = fullfile('/Volumes/curranlab/Data',expName,'eeg/behavioral');
serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data',expName,'eeg/behavioral');
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  uroot = getenv('HOME');
  dataroot = fullfile(uroot,'data',expName,'eeg','behavioral');
end

if exist('prep_eeg','var')
  if length(prep_eeg) == 1
    error('prep_eeg must specify both sessions')
  end
end

if nargin < 3
  prep_eeg = [0 1];
end

if session == 1
  sessions = {'session_0','session_1'};
elseif session == 2
  sessions = {'session_2','session_3'};
end

event_names = {'study_events.mat','test_events.mat'};

% dr = dir(dataroot);
% s = filterStruct(dr,{'strncmp(name,''RSRC'',4)'});
% subjects = getStructField(s,'name');

% if prep_eeg == 1
%   % get bad channel info
%   [subBC,sesBC,badChan] = textread(fullfile(dataroot,'rsrc_badChan.txt'),'%s%s%s','delimiter','\t');
%   for i = 1:length(badChan)
%     badChan{i} = str2num(badChan{i});
%   end
%   infoStruct = struct('subject',subBC,'session',sesBC,'badChan',badChan);
% end

for sub = 1:length(subjects)
  fprintf('Getting data for %s...',subjects{sub});
  for ses = 1:length(sessions)
    fprintf('%s...\n',sessions{ses});
    
    % if prep_eeg == 1
    %   % find the bad channels for this subject and session
    %   sesStruct = filterStruct(infoStruct,'ismember(subject,varargin{1}) & ismember(session,varargin{2})',subjects{sub},sessions{ses});
    %   subSesBadChan = sesStruct.badChan;
    %   if isempty(sesStruct)
    %     error('no subject listing found for this session');
    %   end
    % end
    subSesBadChan = [];
    
    % set the subject events directory
    eventsOutdir_sub = fullfile(dataroot,subjects{sub},sessions{ses},'events');
    if ~exist(eventsOutdir_sub,'dir')
      mkdir(eventsOutdir_sub);
    end
    
    if ~exist(fullfile(dataroot,subjects{sub},sessions{1},'events'),'dir')
      mkdir(fullfile(dataroot,subjects{sub},sessions{1},'events'));
    end
    if ~exist(fullfile(dataroot,subjects{sub},sessions{2},'events'),'dir')
      mkdir(fullfile(dataroot,subjects{sub},sessions{2},'events'));
    end
    
    % save the study and test events separately; only pass the test events
    % to the EEG processing
    
    % set the subject events file
    eventsOutfile_sub = fullfile(eventsOutdir_sub,event_names{ses});
    if exist(eventsOutfile_sub,'file')
      %fprintf('%s already exists! Skipping this subject!\n',eventsOutfile_sub);
      %continue
      
      fprintf('%s already exists! Moving to next step (either prep_eeg or next subject)...\n',eventsOutfile_sub);
      %events = loadEvents(eventsOutfile_sub);
    else
      %if ~lockFile(eventsOutfile_sub)
      fprintf('Creating events...\n');
      % create the events
      [study_events,test_events] = rsrc3_createEvents(dataroot,subjects{sub},session);
      
      fprintf('Saving study and test events for %s...\n',fullfile(dataroot,subjects{sub}));
      
      saveEvents(study_events,fullfile(dataroot,subjects{sub},sessions{1},'events',event_names{1}));
      saveEvents(test_events,fullfile(dataroot,subjects{sub},sessions{2},'events',event_names{2}));
      
      % release the lockFile
      %releaseFile(eventsOutfile_sub);
    end
    
    if prep_eeg(ses) == 1
      fprintf('\nPrepping EEG data...\n');
      % get this subject's session dir
      sesDir = fullfile(dataroot,subjects{sub},sessions{ses});
      % % this is where the RAW file is currently stored
      % oldEEGDir = fullfile(sesDir,'eeg');
      % % this is where we want to move the EEG data to
      % newEEGRoot = fullfile(eegroot,subjects{sub},sessions{ses});
      % newEEGDir = fullfile(newEEGRoot,'eeg');
      % if exist(newEEGDir,'dir')
      %   sprintf('%s already exists, skipping!\n',newEEGDir);
      %   continue
      % else
      %   mkdir(newEEGRoot);
      % end
            
      subEegDir = fullfile(sesDir,'eeg','eeg.noreref');
      pfile = dir(fullfile(subEegDir,[subjects{sub},'*params.txt']));
      
      if ~exist(fullfile(subEegDir,pfile.name),'file')
        curDir = pwd;
        % cd to the session directory since prep_egi_data needs to be there
        cd(sesDir);
        % align the behavioral and EEG data
        prep_egi_data_CU(subjects{sub},sesDir,{eventsOutfile_sub},subSesBadChan,'mstime','HCGSN');
        % go back to the previous working directory
        cd(curDir);
      end
      
      % export the events for netstation; saves to the session's events dir
      rsrc3_events2ns(dataroot,subjects{sub},sessions{ses});
      
      % % change the location of EEG files
      % if ~exist(newEEGDir,'dir')
      %   % do the actual move
      %   unix(sprintf('mv %s %s',oldEEGDir,newEEGRoot));
      %   % resave the events
      %   events = loadEvents(fullfile(sesDir,'events/events.mat'),{oldEEGDir,newEEGDir});
      %   saveEvents(events,fullfile(sesDir,'events/events.mat'));
      % else
      %   error('Error: %s already exists, not moving any directories!\n',newEEGDir);
      % end

    end % prep_eeg
  end % ses
  fprintf('Done.\n');
end
