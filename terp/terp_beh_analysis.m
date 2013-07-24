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
  'TERP003';
  'TERP004';
  'TERP009';
  'TERP010';
  };

sessions = {'session_0','session_1'};

resultsFile = fullfile(dataroot,sprintf('TERP_results_%dsub.txt',length(subjects)));
%resultsFile = fullfile('~/Desktop',sprintf('TERP_results_%dsub.txt',length(subjects)));
rf = fopen(resultsFile,'w');

% legend
fprintf(rf,'QuTr = quiz-trained\n');
fprintf(rf,'ReTr = restudy-trained\n');

fprintf(rf,'HR = hit rate\n');
fprintf(rf,'MR = miss rate\n');
fprintf(rf,'CRR = correct rejection rate\n');
fprintf(rf,'FAR = false alarm rate\n');

fprintf(rf,'RR = recall rate\n');

fprintf(rf,'Recog = recognition\n');
fprintf(rf,'Reten = retention\n');

fprintf(rf,'\n');

% header
fprintf(rf,'subject\tQuizRecallRate\tQuTr Recog HR\tQuTr Recog MR\tReTr Recog HR\tReTr Recog MR\tRecog CRR\tRecog FAR\tQuTr Recog RR\tReTr Recog RR\tQuTr Reten RR\tReTr Reten RR\n');
  
for sub = 1:length(subjects)
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('Getting data for %s (both sessions)...\n',subjects{sub});
  
  eventsDir_sub = fullfile(saveDir,subjects{sub},'events');
  if ~exist(eventsDir_sub,'dir')
    mkdir(eventsDir_sub);
  end
  
  % set the subject events file
  eventsFile0 = fullfile(eventsDir_sub,'events_ses0.mat');
  eventsFile1 = fullfile(eventsDir_sub,'events_ses1.mat');
  if exist(eventsFile0,'file') && exist(eventsFile1,'file')
    events_ses0 = loadEvents(eventsFile0);
    events_ses1 = loadEvents(eventsFile1);
  else
    fprintf('%s and %s do not exist! Next subject...\n',eventsOutfile0_sub,eventsOutfile1_sub);
    continue
  end
  
  % analyses TODO: recall during recognition (session_0) and retention
  % (session_1) as a function of session_0 training type (restudy vs quiz)
  
  % initial study phase target events, divided by training type (quiz or restudy)
  study_targ_quiz = filterStruct(events_ses0,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(type,varargin{2}) & ismember(train_type,varargin{3})',{'study'},{'STUDY_TARGET'},{'quiz'});
  study_targ_restudy = filterStruct(events_ses0,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(type,varargin{2}) & ismember(train_type,varargin{3})',{'study'},{'STUDY_TARGET'},{'restudy'});
  
%   % quiz phase target events
%   quiz_targ_quiz = filterStruct(events_ses0,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(type,varargin{2}) & ismember(train_type,varargin{3})',{'quiz'},{'STUDY_TARGET'},{'quiz'});

%   % recogrecall phase target events, divided by training type (quiz or restudy)
%   recog_targ_quiz = filterStruct(events_ses0,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(type,varargin{2}) & ismember(train_type,varargin{3})',{'recogrecall'},{'TEST_TARGET'},{'quiz'});
%   recog_targ_restudy = filterStruct(events_ses0,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(type,varargin{2}) & ismember(train_type,varargin{3})',{'recogrecall'},{'TEST_TARGET'},{'restudy'});
  
%   % retention phase target events, divided by training type (quiz or restudy)
%   reten_targ_quiz = filterStruct(events_ses1,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(type,varargin{2}) & ismember(train_type,varargin{3})',{'retention'},{'STUDY_TARGET'},{'quiz'});
%   reten_targ_restudy = filterStruct(events_ses1,'rec_isTarg == 1 & ismember(condition,varargin{1}) & ismember(type,varargin{2}) & ismember(train_type,varargin{3})',{'retention'},{'STUDY_TARGET'},{'restudy'});
  
  % study event hits (old called "old" during recogrecall)
  recog_targ_H_quiz = filterStruct(study_targ_quiz,'rec_correct == 1');
  recog_targ_H_restudy = filterStruct(study_targ_restudy,'rec_correct == 1');
  % study event misses (old called "new" during recogrecall)
  recog_targ_M_quiz = filterStruct(study_targ_quiz,'rec_correct == 0');
  recog_targ_M_quiz_sure = filterStruct(recog_targ_M_quiz,'strcmp(new_resp,varargin{1})',{'SURE'});
  recog_targ_M_quiz_maybe = filterStruct(recog_targ_M_quiz,'strcmp(new_resp,varargin{1})',{'MAYBE'});
  recog_targ_M_restudy = filterStruct(study_targ_restudy,'rec_correct == 0');
  recog_targ_M_restudy_sure = filterStruct(recog_targ_M_restudy,'strcmp(new_resp,varargin{1})',{'SURE'});
  recog_targ_M_restudy_maybe = filterStruct(recog_targ_M_restudy,'strcmp(new_resp,varargin{1})',{'MAYBE'});
  
  % recognition/recall phase lure events
  recog_lure = filterStruct(events_ses0,'rec_isTarg == 0 & ismember(condition,varargin{1}) & ismember(type,varargin{2})',{'recogrecall'},{'TEST_LURE'});
  % recogrecall correct rejections (new called "new" during recogrecall)
  recog_lure_CR = filterStruct(recog_lure,'rec_correct == 1');
  recog_lure_CR_sure = filterStruct(recog_lure_CR,'strcmp(new_resp,varargin{1})',{'SURE'});
  recog_lure_CR_maybe = filterStruct(recog_lure_CR,'strcmp(new_resp,varargin{1})',{'MAYBE'});
  % recogrecall false alarms (new called "old" during recogrecall)
  recog_lure_FA = filterStruct(recog_lure,'rec_correct == 0');
  
  % study events correctly recalled during the quiz training phase
  quiz_recall_quiz = filterStruct(study_targ_quiz,'quiz_recall == 1');
  % study events not recalled during the quiz training phase
  quiz_norecall_quiz = filterStruct(study_targ_quiz,'quiz_recall == 0');
  % there is no equivalent for the restudy training phase
  
  % study events correctly recalled during recogrecall, divided by training
  recog_recall_quiz = filterStruct(recog_targ_H_quiz,'rec_recall == 1');
  recog_recall_restudy = filterStruct(recog_targ_H_quiz,'rec_recall == 1');
  % study events not recalled during recogrecall, divided by training
  recog_norecall_quiz = filterStruct(recog_targ_H_quiz,'rec_recall == 0');
  recog_norecall_restudy = filterStruct(recog_targ_H_quiz,'rec_recall == 0');
  
  % study events correctly recalled during retention, divided by training
  reten_recall_quiz = filterStruct(study_targ_quiz,'ret_recall == 1');
  reten_recall_restudy = filterStruct(study_targ_restudy,'ret_recall == 1');
  % study events not recalled during retention, divided by training
  reten_norecall_quiz = filterStruct(study_targ_quiz,'ret_recall == 0');
  reten_norecall_restudy = filterStruct(study_targ_restudy,'ret_recall == 0');
  
  qt_qrr = (length(quiz_recall_quiz) / length(study_targ_quiz));
  fprintf('QuizTrained quiz recall rate = (%d/%d) = %.3f\n',length(quiz_recall_quiz),length(study_targ_quiz),qt_qrr);
  %fprintf('QuizTrained quiz norecall rate = (%d/%d) = %.3f\n',length(quiz_norecall_quiz),length(study_targ_quiz),(length(quiz_norecall_quiz) / length(study_targ_quiz)));
  
  % during recogrecall
  fprintf('\n');
  qt_recog_hr = (length(recog_targ_H_quiz) / length(study_targ_quiz));
  fprintf('QuizTrained recognition HR = (%d/%d) = %.3f\n',length(recog_targ_H_quiz),length(study_targ_quiz),qt_recog_hr);
  qt_recog_mr = (length(recog_targ_M_quiz) / length(study_targ_quiz));
  fprintf('QuizTrained recognition MR = (%d/%d) = %.3f\n',length(recog_targ_M_quiz),length(study_targ_quiz),qt_recog_mr);
  qt_recog_mr_s = (length(recog_targ_M_quiz_sure) / length(recog_targ_M_quiz));
  fprintf('QuizTrained recognition MR sure = (%d/%d) = %.3f\n',length(recog_targ_M_quiz_sure),length(recog_targ_M_quiz),qt_recog_mr_s);
  qt_recog_mr_m = (length(recog_targ_M_quiz_maybe) / length(recog_targ_M_quiz));
  fprintf('QuizTrained recognition MR maybe = (%d/%d) = %.3f\n',length(recog_targ_M_quiz_maybe),length(recog_targ_M_quiz),qt_recog_mr_m);
  
  rt_recog_hr = (length(recog_targ_H_restudy) / length(study_targ_restudy));
  fprintf('RestudyTrained recognition HR = (%d/%d) = %.3f\n',length(recog_targ_H_restudy),length(study_targ_restudy),rt_recog_hr);
  rt_recog_mr = (length(recog_targ_M_restudy) / length(study_targ_restudy));
  fprintf('RestudyTrained recognition MR = (%d/%d) = %.3f\n',length(recog_targ_M_restudy),length(study_targ_restudy),rt_recog_mr);
  rt_recog_mr_s = (length(recog_targ_M_restudy_sure) / length(recog_targ_M_restudy));
  fprintf('RestudyTrained recognition MR sure = (%d/%d) = %.3f\n',length(recog_targ_M_restudy_sure),length(recog_targ_M_restudy),rt_recog_mr_s);
  rt_recog_mr_m = (length(recog_targ_M_restudy_maybe) / length(recog_targ_M_restudy));
  fprintf('RestudyTrained recognition MR maybe = (%d/%d) = %.3f\n',length(recog_targ_M_restudy_maybe),length(recog_targ_M_restudy),rt_recog_mr_m);
  
  fprintf('\n');
  recog_cr = (length(recog_lure_CR) / length(recog_lure));
  fprintf('Recognition CRR = (%d/%d) = %.3f\n',length(recog_lure_CR),length(recog_lure),recog_cr);
  recog_cr_s = (length(recog_lure_CR_sure) / length(recog_lure_CR));
  fprintf('Recognition CRR sure = (%d/%d) = %.3f\n',length(recog_lure_CR_sure),length(recog_lure_CR),recog_cr_s);
  recog_cr_m = (length(recog_lure_CR_maybe) / length(recog_lure_CR));
  fprintf('Recognition CRR maybe = (%d/%d) = %.3f\n',length(recog_lure_CR_maybe),length(recog_lure_CR),recog_cr_m);
  recog_far = (length(recog_lure_FA) / length(recog_lure));
  fprintf('Recognition FAR = (%d/%d) = %.3f\n',length(recog_lure_FA),length(recog_lure),recog_far);
  
  fprintf('\n');
  qt_recog_rr = (length(recog_recall_quiz) / length(recog_targ_H_quiz));
  fprintf('QuizTrained recognition recall rate = (%d/%d) = %.3f\n',length(recog_recall_quiz),length(recog_targ_H_quiz),qt_recog_rr);
  %fprintf('QuizTrained recognition norecall rate = (%d/%d) = %.3f\n',length(recog_norecall_quiz),length(recog_targ_H_quiz),(length(recog_norecall_quiz) / length(recog_targ_H_quiz)));
  rt_recog_rr = (length(recog_recall_restudy) / length(recog_targ_H_restudy));
  fprintf('RestudyTrained recognition recall rate = (%d/%d) = %.3f\n',length(recog_recall_restudy),length(recog_targ_H_restudy),rt_recog_rr);
  %fprintf('RestudyTrained recognition norecall rate = (%d/%d) = %.3f\n',length(recog_norecall_restudy),length(recog_targ_H_restudy),(length(recog_norecall_restudy) / length(recog_targ_H_restudy)));
  
  fprintf('\n');
  qt_reten_rr = (length(reten_recall_quiz) / length(study_targ_quiz));
  fprintf('QuizTrained retention recall rate = (%d/%d) = %.3f\n',length(reten_recall_quiz),length(study_targ_quiz),qt_reten_rr);
  %fprintf('QuizTrained retention norecall rate = (%d/%d) = %.3f\n',length(reten_norecall_quiz),length(study_targ_quiz),(length(reten_norecall_quiz) / length(study_targ_quiz)));
  rt_reten_rr = (length(reten_recall_restudy) / length(study_targ_restudy));
  fprintf('RestudyTrained retention recall rate = (%d/%d) = %.3f\n',length(reten_recall_restudy),length(study_targ_restudy),rt_reten_rr);
  %fprintf('RestudyTrained retention norecall rate = (%d/%d) = %.3f\n',length(reten_norecall_restudy),length(study_targ_restudy),(length(reten_norecall_restudy) / length(study_targ_restudy)));
  
  fprintf('\n');
  
  fprintf(rf,'%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',...
    subjects{sub},qt_qrr,qt_recog_hr,qt_recog_mr,rt_recog_hr,rt_recog_mr,recog_cr,recog_far,qt_recog_rr,rt_recog_rr,qt_reten_rr,rt_reten_rr);
  
%   %% debug
%   
%   % all quiz events
%   evq = filterStruct(events_ses0,'ismember(train_type,varargin{1}) & ismember(condition,varargin{2})',{'quiz'},{'recogrecall'});
%   
%   % all restudy events
%   evrs = filterStruct(events_ses0,'ismember(train_type,varargin{1}) & ismember(condition,varargin{2})',{'restudy'},{'recogrecall'});
end

fclose(rf);
