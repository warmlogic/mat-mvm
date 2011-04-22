subjects = {
    'SLORK001',...
    'SLORK002',...
    'SLORK003',...
    'SLORK004',...
    'SLORK005',...
    'SLORK006',...
    'SLORK007',...
    'SLORK008',...
           };

sessions = {'session_0'};

uname = getenv('USER');
uroot = getenv('HOME');
dataroot = fullfile(uroot,'data/SLORK/behavioral');
saveDir = dataroot;
% saveDir = fullfile(dataroot,'events');
% if ~exist(saveDir,'dir')
%   mkdir(saveDir);
% end

eventsOutfile = fullfile(saveDir,'slork_events.mat');

if ~lockFile(eventsOutfile)
  error('%s already exists! Not continuing!\n',eventsOutfile);
end

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
releaseFile(eventsOutfile);
fprintf('Done.\n');

% % fix event fields
% for sub = 1:length(subjects)
%   fprintf('Getting data for %s...',subjects{sub});
%   for ses = 1:length(sessions)
%     fprintf('%s...',sessions{ses});
    
%     events = slork_fixEvents(dataroot,subjects{sub},sessions{ses});
    
%     saveEvents(events,fullfile(dataroot,subjects{sub},sessions{ses},'events','events.mat'));
    
%   end
%   fprintf('Done.\n');
% end

