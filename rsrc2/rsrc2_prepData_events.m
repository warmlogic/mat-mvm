function rsrc2_prepData_events(subjects,prep_eeg)
% rsrc_prepData_events(subjects,prep_eeg)
%
% Purpose
%   Create behavioral events; if prep_eeg: run prep_egi_data, export
%   NetStation events
%
% Inputs
%   subjects: a cell of subject numbers
%
%   prep_eeg: whether or not to prepare the EEG data (should this
%   be a matrix that corresponds to subjects?)
%
% Outputs
%   Events (struct and NetStation) will be saved in:
%     /Users/matt/data/RSRC/behavioral/subject/session/events
%
% Assumptions
%   Each subject only ran one session (session_0)
%
%   The behavioral data is located in:
%     /Users/matt/data/RSRC/behavioral/subject/session
%
%   The NetStation raw file is located in:
%     /Users/matt/data/RSRC/behavioral/subject/session/eeg
%
%   The EEG data will be moved to:
%     /Users/matt/data/RSRC/eeg/subject/session/eeg/
%     (The events.eegfile field will be modified accordingly.)
%

expName = 'RSRC2';

serverDir = fullfile('/Volumes/curranlab/Data',expName,'eeg/behavioral');
serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data',expName,'eeg/behavioral');
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  uroot = getenv('HOME');
  dataroot = fullfile(uroot,'data',expName,'eeg');
end

if nargin < 2
  prep_eeg = 1;
end

sessions = {'session_0'};

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
    
    % set the subject events file
    eventsOutfile_sub = fullfile(eventsOutdir_sub,'events.mat');
    if exist(eventsOutfile_sub,'file')
      %fprintf('%s already exists! Skipping this subject!\n',eventsOutfile_sub);
      %continue
      fprintf('%s already exists! Loading subject events...\n',eventsOutfile_sub);
      events = loadEvents(eventsOutfile_sub);
    else
      %if ~lockFile(eventsOutfile_sub)
      fprintf('Creating events...\n');
      % create the events
      events = rsrc2_createEvents(dataroot,subjects{sub},sessions{ses});
      fprintf('Saving %s...\n',eventsOutfile_sub);
      % save each subject's events
      saveEvents(events,eventsOutfile_sub);
      
      % release the lockFile
      %releaseFile(eventsOutfile_sub);
    end
    
    if prep_eeg == 1
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
        prep_egi_data_CU(subjects{sub},sesDir,{fullfile(sesDir,'events/events.mat')},subSesBadChan,'mstime','HCGSN');
        % go back to the previous working directory
        cd(curDir);
      end
      
      % export the events for netstation; saves to the session's events dir
      rsrc2_events2ns(dataroot,subjects{sub},sessions{ses});
      
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
