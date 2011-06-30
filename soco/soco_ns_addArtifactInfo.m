% soco_ns_addArtifactInfo

overwriteArtFields = 0;

exper.name = 'SOCO';

exper.subjects = {
  'SOCO001';
  'SOCO002';
  'SOCO003';
  'SOCO004';
  'SOCO005';
  'SOCO006';
  'SOCO007';
  'SOCO008';
  'SOCO009';
  'SOCO010';
  'SOCO011';
  'SOCO012';
  'SOCO013';
  'SOCO014';
  'SOCO015';
  'SOCO016';
  'SOCO017';
  'SOCO018';
  'SOCO019';
  'SOCO020';
  'SOCO021';
  'SOCO022';
  'SOCO023';
  'SOCO024';
  'SOCO025';
  'SOCO026';
  'SOCO027';
  'SOCO028';
  'SOCO029';
  'SOCO030';
  };
% SOCO002 ended early by 6(?) trials because of fire alarm

exper.sessions = {'session_0'};

exper.eventValues = {'CR2','HSC2','HSI2','CR6','HSC6','HSI6'};

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
% CR2 (2 colors)
nsEvFilters.CR2.type = 'LURE_PRES';
nsEvFilters.CR2.testType = 'side';
nsEvFilters.CR2.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% HSC2 (2 colors)
nsEvFilters.HSC2.type = 'TARG_PRES';
nsEvFilters.HSC2.testType = 'side';
nsEvFilters.HSC2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% HSI2 (2 colors)
nsEvFilters.HSI2.type = 'TARG_PRES';
nsEvFilters.HSI2.testType = 'side';
nsEvFilters.HSI2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

% CR6 (6 colors)
nsEvFilters.CR6.type = 'LURE_PRES';
nsEvFilters.CR6.testType = 'task';
nsEvFilters.CR6.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% HSC6 (6 colors)
nsEvFilters.HSC6.type = 'TARG_PRES';
nsEvFilters.HSC6.testType = 'task';
nsEvFilters.HSC6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% HSI6 (6 colors)
nsEvFilters.HSI6.type = 'TARG_PRES';
nsEvFilters.HSI6.testType = 'task';
nsEvFilters.HSI6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,overwriteArtFields);
  end
end
