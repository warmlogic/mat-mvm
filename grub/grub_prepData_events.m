function grub_prepData_events(subjects,prep_eeg)
% grub_prepData_events(subjects,prep_eeg)
%
% Purpose
%   Create behavioral events; if prep_eeg == 1: run prep_egi_data, export
%   NetStation events
%
% Inputs
%   subjects: a cell of subject numbers
%
%   prep_eeg: boolean, whether or not to prepare the EEG data
%
% Outputs
%   Events (struct and NetStation) will be saved in:
%     /Users/username/data/GRUB/subject/session/events
%
% Assumptions
%   Each subject only ran one session (session_0)
%
%   The behavioral data is located in:
%     /Users/username/data/GRUB/subject/session
%
%   The NetStation raw file is located in:
%     /Users/username/data/GRUB/subject/session/eeg
%
%   The EEG data will be moved to:
%     /Users/username/data/GRUB/eeg/subject/session/eeg/
%     (The events.eegfile field will be modified accordingly.)
%

serverDir = '/Volumes/curranlab/Data/GRUB/eeg/behavioral';
serverLocalDir = '/Volumes/RAID/curranlab/Data/GRUB/eeg/behavioral';
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  %uname = getenv('USER');
  uroot = getenv('HOME');
  dataroot = fullfile(uroot,'data/GRUB/eeg');
end
saveDir = dataroot;

if nargin == 0
  subjects = {...
      'GRUB002'...
      'GRUB003'...
      'GRUB004'...
      'GRUB005'...
      'GRUB006'...
      'GRUB007'...
      'GRUB008'...
      'GRUB009'...
      'GRUB010'...
      'GRUB011'...
      'GRUB012'...
      'GRUB014'...
      'GRUB015'...
      'GRUB016'...
      'GRUB017'...
      'GRUB018'...
      'GRUB019'...
      'GRUB020'...
      'GRUB021'...
      'GRUB022'...
      'GRUB023'...
      'GRUB024'...
      'GRUB025'...
      'GRUB026'};
  %    'GRUB001'... % exp was not set up correctly
  %    'GRUB013'... % only half the session was recorded
  
  prep_eeg = 0;
end

sessions = {'session_0'};

if prep_eeg == 1
  eegroot = fullfile(dataroot,'eeg');
  
  % get bad channel info
  fid = fopen(fullfile(dataroot,'grub_badChan.txt'));
  [subBC,sesBC,badChan] = textscan(fid,'%s%s%s','delimiter','\t');
  fclose(fid);
  for i = 1:length(badChan)
    badChan{i} = str2double(badChan{i});
  end
  infoStruct = struct('subject',subBC,'session',sesBC,'badChan',badChan);
end

for sub = 1:length(subjects)
  fprintf('Getting data for %s...',subjects{sub});
  for ses = 1:length(sessions)
    fprintf('%s...\n',sessions{ses});
    
    if prep_eeg == 1
      % find the bad channels for this subject and session
      sesStruct = filterStruct(infoStruct,'ismember(subject,varargin{1}) & ismember(session,varargin{2})',subjects{sub},sessions{ses});
      subSesBadChan = sesStruct.badChan;
      if isempty(sesStruct)
        error('no subject listing found for this session');
      end
    end
    
    % set the subject events directory
    eventsOutdir_sub = fullfile(saveDir,subjects{sub},sessions{ses},'events');
    if ~exist(eventsOutdir_sub,'dir')
      mkdir(eventsOutdir_sub);
    end
    
    % set the subject events file
    eventsOutfile_sub = fullfile(eventsOutdir_sub,'events.mat');
%     if ~lockFile(eventsOutfile_sub)
%       fprintf('%s already exists! Skipping this subject!\n',eventsOutfile_sub);
%       continue
%     end
    
    % create the events
    events = grub_createEvents(dataroot,subjects{sub},sessions{ses});
    
    fprintf('Saving %s...',eventsOutfile_sub);
    % save each subject's events
    saveEvents(events,eventsOutfile_sub);
    % release the lockFile
    %releaseFile(eventsOutfile_sub);
    
    rec_targEv = filterStruct(events,'rec_isTarg == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
    rec_lureEv = filterStruct(events,'rec_isTarg == 0 & ismember(type,varargin{1})',{'LURE_PRES'});
    fprintf('\n\nRecognition accuracy info:\n');
    rec_correct_targ = getStructField(rec_targEv,'rec_correct');
    rec_correct_lure = getStructField(rec_lureEv,'rec_correct');
    fprintf('Hits: %.4f\n',sum(rec_correct_targ == 1) / length(rec_targEv));
    fprintf('Misses: %.4f\n',sum(rec_correct_targ == 0) / length(rec_targEv));
    fprintf('Correct rejections: %.4f\n',sum(rec_correct_lure == 1) / length(rec_lureEv));
    fprintf('False alarms: %.4f\n',sum(rec_correct_lure == 0) / length(rec_lureEv));
    
    src_targEv = filterStruct(rec_targEv,'src_isTarg == 1 & rec_correct == 1');
    src_lureEv = filterStruct(rec_targEv,'src_isTarg == 0 & rec_correct == 1');
    fprintf('\nSource accuracy info:\n');
    src_correct_targ = getStructField(src_targEv,'src_correct');
    src_correct_lure = getStructField(src_lureEv,'src_correct');
    fprintf('Hits: %.4f\n',sum(src_correct_targ == 1) / length(src_targEv));
    fprintf('Misses: %.4f\n',sum(src_correct_targ == 0) / length(src_targEv));
    fprintf('Correct rejections: %.4f\n',sum(src_correct_lure == 1) / length(src_lureEv));
    fprintf('False alarms: %.4f\n',sum(src_correct_lure == 0) / length(src_lureEv));
    
    if prep_eeg == 1
      fprintf('\nPrepping EEG data...\n');
      % get this subject's session dir
      sesDir = fullfile(dataroot,subjects{sub},sessions{ses});
      % this is where the RAW file is currently stored
      oldEEGDir = fullfile(sesDir,'eeg');
      % this is where we want to move the EEG data to
      newEEGRoot = fullfile(eegroot,subjects{sub},sessions{ses});
      newEEGDir = fullfile(newEEGRoot,'eeg');
      if exist(newEEGDir,'dir')
        sprintf('%s already exists, skipping!\n',newEEGDir);
        continue
      else
        mkdir(newEEGRoot);
      end
            
      curDir = pwd;
      % cd to the session directory since prep_egi_data needs to be there
      cd(sesDir);
      % align the behavioral and EEG data
      prep_egi_data(subjects{sub},sesDir,{fullfile(sesDir,'events/events.mat')},subSesBadChan,'mstime','HCGSN');
      % go back to the previous working directory
      cd(curDir);
      
      % export the events for netstation; saves to the session's events dir
      grub_events2ns(dataroot,subjects{sub},sessions{ses});
      
      % change the location of EEG files
      if ~exist(newEEGDir,'dir')
        % do the actual move
        unix(sprintf('mv %s %s',oldEEGDir,newEEGRoot));
        % resave the events
        events = loadEvents(fullfile(sesDir,'events/events.mat'),{oldEEGDir,newEEGDir});
        saveEvents(events,fullfile(sesDir,'events/events.mat'));
      else
        error('Error: %s already exists, not moving any directories!\n',newEEGDir);
      end

    end
  end
  fprintf('Done.\n\n');
end
