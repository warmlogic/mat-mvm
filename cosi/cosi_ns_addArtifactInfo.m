% cosi_ns_addArtifactInfo

overwriteArtFields = 1;

exper.name = 'COSI';

exper.subjects = {
  'COSI001';
  'COSI002';
  'COSI003';
  'COSI004';
  'COSI005';
  'COSI006';
  'COSI007';
%   'COSI008';
%   'COSI009';
%   'COSI010';
  'COSI011';
  'COSI012';
  'COSI013';
  'COSI014';
  'COSI015';
  'COSI016';
  'COSI017';
  'COSI018';
  'COSI020';
  'COSI019';
%   'COSI021';
  'COSI022';
  'COSI023';
  'COSI024';
  'COSI025';
  'COSI026';
  'COSI027';
  'COSI028';
  'COSI029';
  'COSI030';
  'COSI031';
  };

exper.sessions = {'session_0','session_1'};

exper.eventValues = {'CCR','CHSC','CHSI','SCR','SHSC','SHSI'};

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
dirs.dataroot = fullfile(dirs.dataroot,exper.name,'eeg/nspp/-1000_2000');

%% add NS's artifact information to the event structure

evFilters = [];
evFilters.eventValues = exper.eventValues;
% CCR (color)
evFilters.CCR.type = 'TEST_LURE';
evFilters.CCR.cond = 'color';
evFilters.CCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% CHSC (color)
evFilters.CHSC.type = 'TEST_TARGET';
evFilters.CHSC.cond = 'color';
evFilters.CHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% CHSI (color)
evFilters.CHSI.type = 'TEST_TARGET';
evFilters.CHSI.cond = 'color';
evFilters.CHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

% SCR (side)
evFilters.SCR.type = 'TEST_LURE';
evFilters.SCR.cond = 'side';
evFilters.SCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% SHSC (side)
evFilters.SHSC.type = 'TEST_TARGET';
evFilters.SHSC.cond = 'side';
evFilters.SHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% SHSI (side)
evFilters.SHSI.type = 'TEST_TARGET';
evFilters.SHSI.cond = 'side';
evFilters.SHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},evFilters,129,overwriteArtFields);
  end
end

%% reject some other events for behavioral reasons (e.g., RT)
% see eeg_toolbox's filterStruct.m for more info

badFilters = [];
badFilters.eventValues = exper.eventValues;
badFilters.expr = 'src_rt > 4000 | rkn_rt > 4000';
for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    ns_rejectEventsBCI(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},badFilters,129,'nsCategory');
  end
end

%% add NS's artifact information to the event structure

evFilters = [];
evFilters.eventValues = exper.eventValues;
% CCR (color)
evFilters.CCR.type = 'TEST_LURE';
evFilters.CCR.cond = 'color';
evFilters.CCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% CHSC (color)
evFilters.CHSC.type = 'TEST_TARGET';
evFilters.CHSC.cond = 'color';
evFilters.CHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% CHSI (color)
evFilters.CHSI.type = 'TEST_TARGET';
evFilters.CHSI.cond = 'color';
evFilters.CHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

% SCR (side)
evFilters.SCR.type = 'TEST_LURE';
evFilters.SCR.cond = 'side';
evFilters.SCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% SHSC (side)
evFilters.SHSC.type = 'TEST_TARGET';
evFilters.SHSC.cond = 'side';
evFilters.SHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% SHSI (side)
evFilters.SHSI.type = 'TEST_TARGET';
evFilters.SHSI.cond = 'side';
evFilters.SHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},evFilters,129,1);
  end
end
