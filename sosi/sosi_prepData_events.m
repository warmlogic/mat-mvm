function sosi_prepData_events(subjects,prep_eeg)
% sosi_prepData_events(subjects,prep_eeg)
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
%     /Users/username/data/SOSI/subject/session/events
%
% Assumptions
%   Each subject only ran one session (session_0)
%
%   The behavioral data is located in:
%     /Users/username/data/SOSI/subject/session
%

serverDir = '/Volumes/curranlab/Data/SOSI/eeg/behavioral';
serverLocalDir = '/Volumes/RAID/curranlab/Data/SOSI/eeg/behavioral';
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  uroot = getenv('HOME');
  dataroot = fullfile(uroot,'data/SOSI/eeg');
end
saveDir = dataroot;

if nargin == 0
  subjects = {
    'SOSI001';
    'SOSI002';
    'SOSI003';
    'SOSI004';
    'SOSI005';
    'SOSI006';
    'SOSI007';
    'SOSI008';
    'SOSI009';
    'SOSI010';
    'SOSI011';
    'SOSI012';
    'SOSI013';
    'SOSI014';
    'SOSI015';
    'SOSI016';
    'SOSI017';
    'SOSI018';
    'SOSI019'; % reran because the first 19 was bad
    'SOSI020';
    'SOSI021';
    'SOSI022';
    'SOSI023';
    'SOSI024';
    'SOSI025';
    'SOSI026';
    'SOSI027';
    'SOSI028';
    'SOSI029';
    'SOSI030';
    };
%     'SOSI031';
%     'SOSI032';
%     'SOSI033';
%     'SOSI034';
%     'SOSI035';
%     'SOSI036';
%     'SOSI037';
%     'SOSI038';
%     'SOSI039';
%     'SOSI040';
%     'SOSI041';
%     'SOSI042';
%     'SOSI043';
  
  prep_eeg = 1;
end

sessions = {'session_0'};

%matlabpool open

for sub = 1:length(subjects)
  fprintf('Getting data for %s...',subjects{sub});
  for ses = 1:length(sessions)
    fprintf('%s...\n',sessions{ses});
    
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
    eventsOutdir_sub = fullfile(saveDir,subjects{sub},sessions{ses},'events');
    if ~exist(eventsOutdir_sub,'dir')
      mkdir(eventsOutdir_sub);
    end
    
    % set the subject events file
    eventsOutfile_sub = fullfile(eventsOutdir_sub,'events.mat');
    if exist(eventsOutfile_sub,'file')
      %fprintf('%s already exists! Skipping this subject!\n',eventsOutfile_sub);
      %continue
      fprintf('%s already exists! Next subject...\n',eventsOutfile_sub);
      continue
    else
      %if ~lockFile(eventsOutfile_sub)
      fprintf('Creating events...\n');
      % create the events
      events = sosi_createEvents(dataroot,subjects{sub},sessions{ses});
      fprintf('Saving %s...\n',eventsOutfile_sub);
      % save each subject's events
      saveEvents(events,eventsOutfile_sub);
      
      % release the lockFile
      %releaseFile(eventsOutfile_sub);
    end
    
    %% prep the EEG data
    if prep_eeg == 1
      fprintf('Prepping EEG data...\n');
      % get this subject's session dir
      sesDir = fullfile(dataroot,subjects{sub},sessions{ses});
%       % this is where the RAW file is currently stored
%       oldEEGDir = fullfile(sesDir,'eeg');
%       % this is where we want to move the EEG data to
%       newEEGRoot = fullfile(dataroot,subjects{sub},sessions{ses});
%       newEEGDir = fullfile(newEEGRoot,'eeg');
%       if exist(newEEGDir,'dir')
%         sprintf('%s already exists, skipping!\n',newEEGDir);
%         continue
%       else
%         mkdir(newEEGRoot);
%       end
      
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
      sosi_events2ns(dataroot,subjects{sub},sessions{ses});
      
%       % change the location of EEG files
%       if ~exist(newEEGDir,'dir')
%         % do the actual move
%         unix(sprintf('mv %s %s',oldEEGDir,newEEGRoot));
%         % resave the events
%         events = loadEvents(fullfile(sesDir,'events/events.mat'),{oldEEGDir,newEEGDir});
%         saveEvents(events,fullfile(sesDir,'events/events.mat'));
%       else
%         error('Error: %s already exists, not moving any directories!\n',newEEGDir);
%       end

    end % prep_eeg
  end % ses
  fprintf('Done.\n');
end % sub

%matlabpool close

% if do_accuracy
%   if ~exist(fullfile(dataroot,'resp_distrib_head.csv'),'file')
%     fid = fopen(fullfile(dataroot,'resp_distrib_head.csv'),'w');
%     for i = 1:size(headerlines,1)
%       fprintf(fid,'%s',headerlines{i,1});
%       for j = 2:size(headerlines,2)
%         fprintf(fid,',%s',headerlines{i,j});
%       end
%       fprintf(fid,'\n');
%     end
%     fclose(fid);
%   end
%   %xlswrite(fullfile(dataroot,'resp_distrib_head.xls'),headerlines);
%   
%   csvwrite(fullfile(dataroot,'resp_distrib.csv'),resp_data);
% end
