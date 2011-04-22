subjects = {
    'GRUB002',...
    'GRUB003',...
    'GRUB004',...
    'GRUB005',...
    'GRUB006',...
    'GRUB007',...
    'GRUB008',...
    'GRUB009',...
    'GRUB010',...
    'GRUB011',...
    'GRUB012',...
    'GRUB014',...
    'GRUB015',...
    'GRUB016',...
    'GRUB017'...
    'GRUB018'...
    'GRUB019'...
    'GRUB020'...
    'GRUB021'...
    'GRUB022'...
    'GRUB023'...
    'GRUB024'...
    'GRUB025'...
    'GRUB026'...
           };
%    'GRUB001'... % exp was not set up correctly
%    'GRUB013'... % only half the session was recorded

sessions = {'session_0'};

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

eventsOutfile = fullfile(saveDir,'grub_events.mat');

%if ~lockFile(eventsOutfile)
%  error('%s already exists! Not continuing!\n',eventsOutfile);
%end

combinedEvents = [];

for sub = 1:length(subjects)
  fprintf('Getting data for %s...',subjects{sub});
  for ses = 1:length(sessions)
    fprintf('%s...',sessions{ses});
    
    % set the subject events directory
    eventsDir_sub = fullfile(dataroot,subjects{sub},sessions{ses},'events');
    events = loadEvents(fullfile(eventsDir_sub,'events.mat'));
    
    % concatenate every subject's events
    combinedEvents = [combinedEvents;events];
  end
  fprintf('Done.\n');
end

fprintf('Saving %d subjects in %s...',length(subjects),eventsOutfile);
% save all subs
saveEvents(combinedEvents,eventsOutfile);
% release the lockFile
%releaseFile(eventsOutfile);
fprintf('Done.\n');

% % fix event fields
% for sub = 1:length(subjects)
%   fprintf('Getting data for %s...',subjects{sub});
%   for ses = 1:length(sessions)
%     fprintf('%s...',sessions{ses});
    
%     events = grub_fixEvents(dataroot,subjects{sub},sessions{ses});
    
%     saveEvents(events,fullfile(dataroot,subjects{sub},sessions{ses},'events','events.mat'));
    
%   end
%   fprintf('Done.\n');
% end

