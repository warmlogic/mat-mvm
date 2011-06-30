% cosi_ns_addArtifactInfo

overwriteArtFields = 0;

exper.name = 'COSI';

exper.subjects = {
%   'COSI001';
%   'COSI002';
%   'COSI003';
%   'COSI004';
%   'COSI005';
%   'COSI006';
%   'COSI007';
%   'COSI008';
%   'COSI009';
%   'COSI010';
  'COSI011';
  'COSI012';
  'COSI013';
  'COSI014';
  'COSI015';
  'COSI016';
%   'COSI017';
  'COSI018';
%   'COSI020';
%   'COSI019';
%   'COSI021';
%   'COSI022';
%   'COSI023';
%   'COSI024';
%   'COSI025';
%   'COSI026';
%   'COSI027';
%   'COSI028';
%   'COSI029';
%   'COSI030';
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

% add NS's artifact information to the event structure
nsEvFilters = [];
nsEvFilters.eventValues = exper.eventValues;
% CCR (color)
nsEvFilters.CCR.type = 'TEST_LURE';
nsEvFilters.CCR.cond = 'color';
nsEvFilters.CCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% CHSC (color)
nsEvFilters.CHSC.type = 'TEST_TARGET';
nsEvFilters.CHSC.cond = 'color';
nsEvFilters.CHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% CHSI (color)
nsEvFilters.CHSI.type = 'TEST_TARGET';
nsEvFilters.CHSI.cond = 'color';
nsEvFilters.CHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

% SCR (side)
nsEvFilters.SCR.type = 'TEST_LURE';
nsEvFilters.SCR.cond = 'side';
nsEvFilters.SCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% SHSC (side)
nsEvFilters.SHSC.type = 'TEST_TARGET';
nsEvFilters.SHSC.cond = 'side';
nsEvFilters.SHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% SHSI (side)
nsEvFilters.SHSI.type = 'TEST_TARGET';
nsEvFilters.SHSI.cond = 'side';
nsEvFilters.SHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,overwriteArtFields);
  end
end
