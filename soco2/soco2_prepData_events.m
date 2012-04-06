function soco2_prepData_events(subjects,prep_eeg)
% soco2_prepData_events(subjects,prep_eeg)
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
%     /Users/username/data/SOCO2/subject/session/events
%
% Assumptions
%   Each subject ran two sessions (session_0, session_1)
%
%   The behavioral data is located in:
%     /Volumes/curranlab/Data/SOCO2/behavioral/subject/session
%

expName = 'SOCO2';

serverDir = fullfile(filesep,'Volumes','curranlab','Data',expName,'behavioral');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data',expName,'behavioral');
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  dataroot = fullfile(getenv('HOME'),'data',expName,'behavioral');
end
%saveDir = dataroot;

if ~exist('subjects','var') || isempty(subjects)
  subjects = {
%     % longer study durations; redo option; source on one side;
%     % (500ms preview, 2000ms, 1125+-125ms ISI)
%     'SOCO2001'; % didn't perform well; very fidgety
%     'SOCO2002';
%     'SOCO2003';
%     'SOCO2004';
%     'SOCO2005'; % no F SC responses
%     'SOCO2006';
%     'SOCO2007';
%     'SOCO2008';
%     'SOCO2009';
%     'SOCO2010';
%     'SOCO2011';
%     'SOCO2012';
%     'SOCO2013';
%     'SOCO2014';
%     'SOCO2015';
%     'SOCO2016';
%     'SOCO2017';
%     'SOCO2018';
%     'SOCO2019';
%     'SOCO2020';
%     'SOCO2021';
%     'SOCO2022';
%     'SOCO2023';
%     'SOCO2024';
%     'SOCO2025';
%     'SOCO2026';
%     'SOCO2027';
        'SOCO2028'; % 2x study pres; 500 prev, 2000 study, 1125 ISI
        'SOCO2029';
        'SOCO2030';
        'SOCO2031';
        'SOCO2032';
        'SOCO2033';
        'SOCO2034';
        'SOCO2035';
        'SOCO2036';
        'SOCO2037';
        'SOCO2038';
        'SOCO2039';
        'SOCO2040';
        'SOCO2041';
        'SOCO2042';
        'SOCO2043';
        'SOCO2044';
        'SOCO2045';
        'SOCO2046';
        'SOCO2047';
        'SOCO2048';
        'SOCO2049';
        'SOCO2050';
    };
end

if ~exist('prep_eeg','var') || isempty(prep_eeg)
  prep_eeg = false;
end

fprintf('Processing %d subjects',length(subjects));
if prep_eeg
  fprintf(' and aligning EEG data.\n');
else
  fprintf(', behavioral data only.\n');
end

sessions = {'session_0'};

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
      events = soco2_createEvents(dataroot,subjects{sub},sessions{ses},nsFile);
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
      soco2_events2ns(dataroot,subjects{sub},sessions{ses});
      
    end % prep_eeg
  end % ses
  fprintf('Done.\n');
end % sub

% matlabpool close
