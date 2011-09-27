% soco_ns_addArtifactInfo

overwriteArtFields = 1;

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
dirs.dataroot = fullfile(dirs.dataroot,exper.name,'eeg/eppp/-1000_2000');

%% add NS's artifact information to the event structure
evFilters = [];
evFilters.eventValues = exper.eventValues;
% CR2 (2 colors)
evFilters.CR2.type = 'LURE_PRES';
evFilters.CR2.filters = {'rec_isTarg == 0', 'rec_correct == 1', 'numColors == 2'};
% HSC2 (2 colors)
evFilters.HSC2.type = 'TARG_PRES';
evFilters.HSC2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1', 'numColors == 2'};
% HSI2 (2 colors)
evFilters.HSI2.type = 'TARG_PRES';
evFilters.HSI2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0', 'numColors == 2'};

% CR6 (6 colors)
evFilters.CR6.type = 'LURE_PRES';
evFilters.CR6.filters = {'rec_isTarg == 0', 'rec_correct == 1', 'numColors == 6'};
% HSC6 (6 colors)
evFilters.HSC6.type = 'TARG_PRES';
evFilters.HSC6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1', 'numColors == 6'};
% HSI6 (6 colors)
evFilters.HSI6.type = 'TARG_PRES';
evFilters.HSI6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0', 'numColors == 6'};

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
% %% add NS's artifact information to the event structure
% evFilters = [];
% evFilters.eventValues = exper.eventValues;
% % CR2 (2 colors)
% evFilters.CR2.type = 'LURE_PRES';
% evFilters.CR2.filters = {'rec_isTarg == 0', 'rec_correct == 1', 'numColors == 2'};
% % HSC2 (2 colors)
% evFilters.HSC2.type = 'TARG_PRES';
% evFilters.HSC2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1', 'numColors == 2'};
% % HSI2 (2 colors)
% evFilters.HSI2.type = 'TARG_PRES';
% evFilters.HSI2.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0', 'numColors == 2'};
% 
% % CR6 (6 colors)
% evFilters.CR6.type = 'LURE_PRES';
% evFilters.CR6.filters = {'rec_isTarg == 0', 'rec_correct == 1', 'numColors == 6'};
% % HSC6 (6 colors)
% evFilters.HSC6.type = 'TARG_PRES';
% evFilters.HSC6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1', 'numColors == 6'};
% % HSI6 (6 colors)
% evFilters.HSI6.type = 'TARG_PRES';
% evFilters.HSI6.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0', 'numColors == 6'};
% 
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sessions)
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},evFilters,129,1);
%   end
% end
