% very basic script for looking at TERP behavioral data

expName = 'TERP';

serverDir = fullfile(filesep,'Volumes','curranlab','Data',expName,'beh');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data',expName,'beh');
localDir = fullfile(getenv('HOME'),'data',expName,'beh');
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
elseif exist(localDir,'dir')
  dataroot = localDir;
else
  error('No data directory found.');
end
saveDir = dataroot;

subjects = {
  'TERP001';
  'TERP002';
  };

sessions = {'session_0','session_1'};

for sub = 1:length(subjects)
  fprintf('Getting data for %s (both sessions)...',subjects{sub});
  
  eventsDir_sub = fullfile(saveDir,subjects{sub},'events');
  if ~exist(eventsDir_sub,'dir')
    mkdir(eventsDir_sub);
  end
  
  % set the subject events file
  eventsFile0_sub = fullfile(eventsDir_sub,'events_ses0.mat');
  eventsFile1_sub = fullfile(eventsDir_sub,'events_ses1.mat');
  if exist(eventsFile0_sub,'file') && exist(eventsFile1_sub,'file')
    events_ses0 = loadEvents(eventsFile0_sub);
    events_ses1 = loadEvents(eventsFile1_sub);
  else
    fprintf('%s and %s do not exist! Next subject...\n',eventsOutfile0_sub,eventsOutfile1_sub);
    continue
  end
  
  % analyses TODO: recall during recognition (session_0) and retention
  % (session_1) as a function of session_0 training type (restudy vs quiz)
  
  keyboard
  
  % initial study phase target events, divided by training type (quiz or restudy)
  study_targets_quiz = filterStruct(events_ses0,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(train_type,varargin{2})',{'study'},{'quiz'});
  study_targets_restudy = filterStruct(events_ses0,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(train_type,varargin{2})',{'study'},{'restudy'});
  
  % recognition/recall phase lure events
  recog_lures = filterStruct(events_ses0,'rec_isTarg == 0 & ismember(condition,varargin{1})',{'recogrecall'});
  
  % study event hits (called "old" during recogrecall)
  recog_hit_quiz = filterStruct(study_targets_quiz,'rec_correct == 1');
  recog_hit_restudy = filterStruct(study_targets_restudy,'rec_correct == 1');
  
  % study events correctly recalled during the quiz training phase
  quiz_recall_quiz = filterStruct(study_targets_quiz,'quiz_recall == 1');
  % there is no quivalent for the restudy training phase
  
  % study events correctly recalled during recogrecall, divided by training
  recog_recall_quiz = filterStruct(study_targets_quiz,'rec_recall == 1');
  recog_recall_restudy = filterStruct(study_targets_restudy,'rec_recall == 1');
  
  % study events correctly recalled during retention, divided by training
  reten_recall_quiz = filterStruct(study_targets_quiz,'ret_recall == 1');
  reten_recall_restudy = filterStruct(study_targets_restudy,'ret_recall == 1');
  
  %% debug
  
  evq = filterStruct(events_ses0,'ismember(train_type,varargin{1}) & ismember(condition,varargin{2})',{'quiz'},{'recogrecall'});
  
  evrs = filterStruct(events_ses0,'ismember(train_type,varargin{1}) & ismember(condition,varargin{2})',{'restudy'},{'recogrecall'});
  
end

