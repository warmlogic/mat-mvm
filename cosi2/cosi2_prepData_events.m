function cosi2_prepData_events(subjects,prep_eeg)
% cosi2_prepData_events(subjects,prep_eeg)
%
% Purpose
%   Create behavioral events; if prep_eeg==true: export Net Station events
%
% Inputs
%   subjects: a cell of subject numbers
%
%   prep_eeg: boolean, whether or not to prepare the EEG data
%             (default: true)
%
% Outputs
%   Events (struct and NetStation) will be saved in:
%     /Users/username/data/COSI2/subject/session/events
%
% Assumptions
%   Each subject ran two sessions (session_0, session_1)
%
%   The behavioral data is located in:
%     /Volumes/curranlab/Data/COSI2/eeg/behavioral/subject/session
%

% NB:
% COSI2 subjects 8-12 have a few screwed up things in session.log.old, but
% they have been fixed in session.log, which is what gets read. 009 is
% the most messed up. The keys that were recorded (and the responses) for Y
% and B source responses were reversed. So subject 009 thought X was for B
% and / was for Y (because that's what they saw on the screen), but it
% recorded it as X=Y and /=B. This should be fixed to score the subject
% correctly. Subjects 8, 10, 11, and 12 had their color response for false
% alarms messed up (color response to a new item). I don't think this
% affects anything behaviorally, because we only look at correctly
% identified new items (correct rejections), but anyway: 8, 9, 10, and 12
% had all incorrect false alarm responses recorded as "Yellow", and 11 had
% them all recorded as "Blue" (the recorded keys were correct, so this
% should be easy to fix if we need to). 9 would have been the same as
% yellow but the colors got mixed up.
%
% Again, these have been fixed on the server using cosi2_fixLogs_8to12.m.

expName = 'COSI2';

serverDir = fullfile(filesep,'Volumes','curranlab','Data',expName,'eeg','behavioral');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data',expName,'eeg','behavioral');
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  dataroot = fullfile(getenv('HOME'),'data',expName,'eeg','behavioral');
end
%saveDir = dataroot;

if ~exist('subjects','var') || isempty(subjects)
  subjects = {
    % %     'COSI2001'; % short study durations began with COSI2001
    % %     % (500ms preview, 1000ms, 625+-125ms ISI; too fast)
    % % COSI2001 did not finish
    %     'COSI2002';
    %     'COSI2003';
    %     'COSI2004';
    %     'COSI2005';
    %     'COSI2006';
    %     'COSI2007';
    'COSI2008'; % longer study durations began with COSI2008
    %(500ms preview, 2000ms, 1125+-125ms ISI)
    'COSI2009';
    'COSI2010';
    'COSI2011'; % COSI2011 did not do session_1 (didn't like EEG)
    'COSI2012';
    'COSI2013'; % redo option and source-on-one-side began with COSI2013
    'COSI2014'; % COSI2014 did not do session_1 (didn't perform well)
    'COSI2015';
    'COSI2016';
    'COSI2017';
    'COSI2018';
    'COSI2019';
    'COSI2020';
    'COSI2021';
    'COSI2022';
    'COSI2023';
    'COSI2024';
    'COSI2025';
    'COSI2026';
    'COSI2027';
    'COSI2028';
    'COSI2029';
    'COSI2030';
    'COSI2031'; % Thought EEG reference electrode messed up. No session_1.
    'COSI2032';
    'COSI2033';
    'COSI2034';
    'COSI2035';
    'COSI2036';
    'COSI2037';
    'COSI2038'; % COSI2038: potentially bad session_1 (bathroom, sick)
    'COSI2039';
    'COSI2040';
    'COSI2041'; % COSI2041: no-show, no session_1
    'COSI2042';
    'COSI2043';
    'COSI2044';
    'COSI2045';
    };
end

if ~exist('prep_eeg','var') || isempty(prep_eeg)
  prep_eeg = true;
end

fprintf('Processing %d subjects',length(subjects));
if prep_eeg
  fprintf(' and aligning EEG data.\n');
else
  fprintf(', behavioral data only.\n');
end

sessions = {'session_0','session_1'};
%sessions = {'session_0'};
%sessions = {'session_1'};

% matlabpool open

for sub = 1:length(subjects)
  for ses = 1:length(sessions)
    fprintf('Getting data for %s, %s...\n',subjects{sub},sessions{ses});
    sesDir = fullfile(dataroot,subjects{sub},sessions{ses});
    
    if ~exist(fullfile(sesDir,'session.log'),'file')
      fprintf('Data for %s does not exist (%s). Moving on.\n',sessions{ses},fullfile(sesDir,'session.log'));
      continue
    end
    
    if prep_eeg
      % find the bad channels for this subject and session
%       sesStruct = filterStruct(infoStruct,'ismember(subject,varargin{1}) & ismember(session,varargin{2})',subjects{sub},sessions{ses});
%       subSesBadChan = sesStruct.badChan;
%       if isempty(sesStruct)
%         error('no subject listing found for this session');
%       end
      %subSesBadChan = [];
      
      nsFile = dir(fullfile(sesDir,'eeg','*.raw'));
      if isempty(nsFile)
        fprintf('Did not find %s raw NS file. Moving on.\n',subjects{sub})
        continue
      else
        [~,nsFile] = fileparts(nsFile.name);
      end
    else
      nsFile = '';
    end
    
    % set the subject events directory
    eventsOutdir_sub = fullfile(sesDir,'events');
    if ~exist(eventsOutdir_sub,'dir')
      mkdir(eventsOutdir_sub);
    end
    
    % set the subject events file
    eventsOutfile_sub = fullfile(eventsOutdir_sub,'events.mat');
    if exist(eventsOutfile_sub,'file')
      fprintf('Events file (%s) already exists! Moving on to next subject...\n',eventsOutfile_sub);
      continue
    else
      %if ~lockFile(eventsOutfile_sub)
      fprintf('Creating events...\n');
      % create the events
      events = cosi2_createEvents(dataroot,subjects{sub},sessions{ses},nsFile);
      fprintf('Saving %s...\n',eventsOutfile_sub);
      % save each subject's events
      saveEvents(events,eventsOutfile_sub);
      
      % release the lockFile
      %releaseFile(eventsOutfile_sub);
    end
    
    %% prep the EEG data
    if prep_eeg
      fprintf('Prepping EEG data...\n');
      % get this subject's session dir
      subEegDir = fullfile(sesDir,'eeg','eeg.noreref');
      pfile = dir(fullfile(subEegDir,[subjects{sub},'*params.txt']));
      
      if ~exist(fullfile(subEegDir,pfile.name),'file')
        curDir = pwd;
        % cd to the session directory since prep_egi_data needs to be there
        cd(sesDir);
        % align the behavioral and EEG data
        prep_egi_data_CU(subjects{sub},sesDir,{fullfile(sesDir,'events/events.mat')},[],'mstime','HCGSN');
        % go back to the previous working directory
        cd(curDir);
      end
      
      % export the events for netstation; saves to the session's events dir
      cosi2_events2ns(dataroot,subjects{sub},sessions{ses});
      
    end % prep_eeg
  end % ses
  fprintf('Done.\n');
end % sub

% matlabpool close
