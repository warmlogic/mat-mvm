% sosi_ns_addArtifactInfo

overwriteArtFields = 1;

exper.name = 'SOSI';

exper.subjects = {
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
  'SOSI020';
  'SOSI019';
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
% original SOSI019 was replaced because the first didn't finish

exper.sessions = {'session_0'};

exper.eventValues = {'RCR','RHSC','RHSI'};

% pick the right dirs.dataroot
dirs.serverDir = fullfile('/Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile('/Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile('/data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'data');
if exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
elseif exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
elseif exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
elseif exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
else
  error('Data directory not found.');
end
dirs.dataroot = fullfile(dirs.dataroot,exper.name,'eeg/eppp/-1000_2000');

%% add NS's artifact information to the event structure

evFilters = [];
evFilters.eventValues = exper.eventValues;
% RCR
evFilters.RCR.type = 'TEST_LURE';
evFilters.RCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% RHSC
evFilters.RHSC.type = 'TEST_TARGET';
evFilters.RHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% RHSI
evFilters.RHSI.type = 'TEST_TARGET';
evFilters.RHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},evFilters,129,overwriteArtFields);
  end
end

% %% reject some other events for behavioral reasons (e.g., RT)
% % see eeg_toolbox's filterStruct.m for more info
% 
% badFilters = [];
% badFilters.eventValues = exper.eventValues;
% badFilters.expr = 'src_rt > 4000 | rkn_rt > 4000';
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sessions)
%     ns_rejectEventsBCI(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},badFilters,129,'nsCategory');
%   end
% end
% 
% %% re-add NS's artifact information to the event structure
% 
% evFilters = [];
% evFilters.eventValues = exper.eventValues;
% % RCR
% evFilters.RCR.type = 'TEST_LURE';
% evFilters.RCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % RHSC
% evFilters.RHSC.type = 'TEST_TARGET';
% evFilters.RHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % RHSI
% evFilters.RHSI.type = 'TEST_TARGET';
% evFilters.RHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% 
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sessions)
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},evFilters,129,1);
%   end
% end
