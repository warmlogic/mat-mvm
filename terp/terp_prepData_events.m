function terp_prepData_events(subjects,prep_eeg)
% terp_prepData_events(subjects,prep_eeg)
%
% Purpose
%   Create behavioral events; if prep_eeg == 1: export Net Station events
%
% Inputs
%   subjects: a cell of subject numbers
%
%   prep_eeg: boolean, whether or not to prepare the EEG data
%
% Outputs
%   Events (struct and NetStation) will be saved in:
%     ~/data/TERP/beh/subject/events
%
% Assumptions
%   Each subject ran in both sessions (session_0 and session_1)
%
%   The behavioral data is located in:
%     ~/data/TERP/beh/subject/session
%
% NB: EEG processing and NS event creation are not done yet
%

expName = 'TERP';

serverDir = fullfile(filesep,'Volumes','curranlab','Data',expName,'beh');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data',expName,'beh');
localDir = fullfile(getenv('HOME'),'data',expName,'beh');
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
elseif exist(localDir,'dir')
  dataroot = localDir;
else
  error('No data directory found.');
end
saveDir = dataroot;

if nargin == 0
  subjects = {
    'TERP001';
    'TERP002';
    };
  
  prep_eeg = 0;
end

sessions = {'session_0','session_1'};

%matlabpool open

for sub = 1:length(subjects)
  fprintf('Getting data for %s (both sessions)...\n',subjects{sub});
  
  if prep_eeg == 1
    % find the bad channels for this subject and session
    %       sesStruct = filterStruct(infoStruct,'ismember(subject,varargin{1}) & ismember(session,varargin{2})',subjects{sub},sessions{ses});
    %       subSesBadChan = sesStruct.badChan;
    %       if isempty(sesStruct)
    %         error('no subject listing found for this session');
    %       end
    subSesBadChan = [];
  end
  
  % set the subject events directory
  eventsOutdir_sub = fullfile(saveDir,subjects{sub},'events');
  if ~exist(eventsOutdir_sub,'dir')
    mkdir(eventsOutdir_sub);
  end
  
  % set the subject events file
  eventsOutfile0_sub = fullfile(eventsOutdir_sub,'events_ses0.mat');
  eventsOutfile1_sub = fullfile(eventsOutdir_sub,'events_ses1.mat');
  if exist(eventsOutfile0_sub,'file') && exist(eventsOutfile1_sub,'file')
    %fprintf('%s already exists! Skipping this subject!\n',eventsOutfile_sub);
    %continue
    fprintf('%s and %s already exist! Next subject...\n',eventsOutfile0_sub,eventsOutfile1_sub);
    continue
  else
    %if ~lockFile(eventsOutfile_sub)
    fprintf('Creating events...\n');
    % create the events
    [events_ses0, events_ses1] = terp_createEvents(dataroot,subjects{sub},sessions);
    fprintf('Saving %s...\n',eventsOutfile0_sub);
    % save each subject's events
    saveEvents(events_ses0,eventsOutfile0_sub);
    fprintf('Saving %s...\n',eventsOutfile1_sub);
    saveEvents(events_ses1,eventsOutfile1_sub);
    
    % release the lockFile
    %releaseFile(eventsOutfile_sub);
  end
  
  %% prep the EEG data
  if prep_eeg == 1
    for ses = 1:length(sessions)
      sesNum = str2double(strrep(sessions{ses},'session_',''));
      fprintf('Prepping EEG data for %s...\n',sessions{ses});
      % get this subject's session dir
      subDir = fullfile(dataroot,subjects{sub});
      sesDir = fullfile(subDir,sessions{ses});
      
      subEegDir = fullfile(sesDir,'eeg','eeg.noreref');
      pfile = dir(fullfile(subEegDir,[subjects{sub},'*params.txt']));
      
      if ~exist(fullfile(subEegDir,pfile.name),'file')
        curDir = pwd;
        % cd to the session directory since prep_egi_data needs to be there
        cd(sesDir);
        % align the behavioral and EEG data
        %
        % TODO: make sure this works
        prep_egi_data_CU(subjects{sub},sesDir,{fullfile(subDir,'events',sprintf('events_ses%d.mat',sesNum))},subSesBadChan,'mstime','HCGSN');
        % go back to the previous working directory
        cd(curDir);
      end
      
      % export the events for netstation; saves to the session's events dir
      %
      % TODO: create/fix terp_events2ns so it doesn't rely on events in
      % session directory
      fprintf('terp_events2ns IS NOT DONE YET!!!\n');
      terp_events2ns(dataroot,subjects{sub},sessions{ses});
      
    end % ses
    
  end % prep_eeg
  fprintf('Done.\n');
end % sub

%matlabpool close
