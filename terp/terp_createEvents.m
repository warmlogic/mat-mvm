function [events_ses0, events_ses1] = terp_createEvents(dataroot,subject,sessions)
% function [events_ses0, events_ses1] = terp_createEvents(dataroot,subject,sessions)
%
% e.g., [events_ses0, events_ses1] = terp_createEvents('/Volumes/curranlab/Data/TERP/beh/','TERP002',{'session_0','session_1'});
%
% create event struct for TERP
%
% struct fields:
%   subject
%   session
%   mstime
%   msoffset
%   list (#)
%   condition (study, restudy, quiz, recogrecall, retention)
%   trial (#)
%   type
%   train_type (restudy or quiz)
%   item_swa
%   xcoord_swa
%   ycoord_swa
%   item_eng
%   xcoord_eng
%   ycoord_eng
%   serialpos (initial study SP)
%   rec_isTarg (1 or 0)
%   rec_correct (1 or 0)
%   rec_resp
%   rec_resp_rt
%   new_correct
%   new_resp (sure, maybe, -1)
%   new_resp_rt (sure, maybe, -1)
%   quiz_recall (1 or 0)
%   rec_recall (1 or 0)
%   ret_recall (1 or 0)
%   quiz_recall_resp
%   rec_recall_resp
%   ret_recall_resp
%
% [events_ses0, events_ses1] = terp_createEvents('~/data/TERP/beh/','TERP002',{'session_0','session_1'});
%

% event struct recall notation:
%
% for quiz recall, -1 indicates that there was no word recalled; the
% corresponding 'resp' field should be blank (''). 0 indicates that they
% recalled the wrong word, inclduing vocalizations ('<>'). 1 indicates that
% they recalled the correct word.
%
% for recognition recall, -1 indicates that they did not have a chance to
% recall (they incorrectly called the old word "new"), and the resp field
% gets set to blank (''). 0 indicates that they recalled the wrong word
% (including vocalizations). 1 indicates that they recalled the correct
% word.
%
% retention recall should be the same as quiz recall.

% TODO: be able to process only the first session; e.g., for making sure
% the subject is doing a decent job as recalling during quiz and
% recognition

% % debug
% % dataroot = '/Volumes/curranlab/Data/TERP/beh/';
% dataroot = '~/data/TERP/beh/';
% subject = 'TERP001';
% % subject = 'TERP002';
% sessions = {'session_0','session_1'};
% nargout =2;

%% error checking

if nargout ~= 2
  error('%s currently only supports processing both session_0 and session_1. Therefore, output must be something like [events_ses0, events_ses1].',mfilename);
end

if ~iscell(sessions)
  error('input variable ''sessions'' must be a cell, specifically {''session_0'',''session_1''}.');
end

if length(sessions) ~= 2
  error('%s currently only supports processing both session_0 and session_1. Therefore, input variable ''sessions'' must be {''session_0'',''session_1''}.',mfilename);
end

%% set up some constants and initialize some variables
numLists = 12;
maxTestTime = 10000;
quizAnnFile = 'train_eng_';
recogAnnFile = 'recall_targets_eng_';
retenAnnFile = 'reten_eng_';
maxSpeakDur = 3000;

% vocalization symbol in annotation files
vocal = '<>';

%fnames = {'subject', 'session', 'mstime', 'msoffset', 'list', 'condition', 'trial', 'type', 'train_type', 'item_swa', 'xcoord_swa', 'ycoord_swa', 'item_eng', 'xcoord_eng', 'ycoord_eng', 'serialpos', 'rec_isTarg', 'rec_correct', 'rec_resp', 'rec_resp_rt', 'new_correct', 'new_resp', 'new_resp_rt', 'quiz_recall', 'rec_recall', 'ret_recall', 'quiz_recall_resp', 'rec_recall_resp', 'ret_recall_resp'};
% , 'quiz_recall_rt', 'rec_recall_rt', 'ret_recall_ms'

%events = cell2struct(cell(1,length(fnames)),fnames,2);
%events = struct;

% initialize the annotation file start times; half of the quiz start times
% will remain nans because only half of the lists get quizzes
ann_start_time_quiz = nan(1,numLists);
ann_start_time_recog = nan(1,numLists);
ann_start_time_reten = nan(1,numLists);

%% start the event creation

for ses = 1:length(sessions)
  session = sessions{ses};
  sesDir = fullfile(dataroot,subject,session);
  
  logFile = fullfile(sesDir,'session.log');
  if exist(logFile,'file')
    fid = fopen(logFile);
    logData = textscan(fid,'%n%n%s%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
    fclose(fid);
  else
    error('session.log file not found');
  end
  
  sesNum = str2double(strrep(session,'session_',''));
  
  % set all fields here so we can easily concatenate events later
  log = struct('subject',subject,'session',sesNum,'mstime',num2cell(logData{1}),'msoffset',num2cell(logData{2}),...
    'list',[],'condition',[],'trial',[],'type',logData{3},'train_type',[],...
    'item_swa',[], 'xcoord_swa',[], 'ycoord_swa',[],...
    'item_eng',[], 'xcoord_eng',[], 'ycoord_eng',[],...
    'serialpos',[], 'rec_isTarg',[],...
    'rec_correct',[], 'rec_resp','', 'rec_resp_rt',[],...
    'new_correct',[], 'new_resp',[], 'new_resp_rt',[],...
    'quiz_recall',[], 'quiz_recall_resp',[],...
    'rec_recall',[], 'rec_recall_resp',[],...
    'ret_recall',[], 'ret_recall_resp',[]);
  % nb: rec_resp is set to '' so we can filter on it later
  
  % initialize the list number
  listNum = nan;
  
  for i = 1:length(log)
    log(i).mstime = log(i).mstime;
    log(i).msoffset = log(i).msoffset;
    
    % set the list number (confusingly, list is denoted with TRIAL in log)
    if strcmp('TRIAL',logData{3}{i}) % TRIAL = listNum
      listNum = str2double(logData{4}{i});
      % if we're starting a new list, set trialNum to zero
      trialNum_study = 0;
      trialNum_train = 0;
      trialNum_recog = 0;
      trialNum_reten = 0;
    end
    log(i).list = listNum;
    
    switch log(i).type
      
      % collect the annotation start times (REC_START in session.log)
      case {'STUDY_QUIZ_START'}
        ann_start_time_quiz(log(i).list) = logData{1}(i+1);
        
      case {'TEST_BEGIN_KEYHIT'}
        ann_start_time_recog(log(i).list) = logData{1}(i+1);
        
      case {'RETENTION_TEST_START'}
        ann_start_time_reten(log(i).list) = logData{1}(i+1);
        
      case {'STUDY_TARGET','STUDY_BUFFER'}
        % set training type based on the studystudy or studyquiz definition
        if strcmp(logData{4}{i},'studystudy')
          log(i).train_type = 'restudy';
        elseif strcmp(logData{4}{i},'studyquiz')
          log(i).train_type = 'quiz';
        end
        % then set the current condition based on the condition phase number
        % (everyone always does study first)
        if str2double(logData{5}{i}) == 0
          log(i).condition = 'study';
          trialNum_study = trialNum_study + 1;
          log(i).trial = trialNum_study;
        elseif str2double(logData{5}{i}) == 1
          if strcmp(logData{4}{i},'studystudy')
            log(i).condition = 'restudy';
            trialNum_train = trialNum_train + 1;
            log(i).trial = trialNum_train;
          elseif strcmp(logData{4}{i},'studyquiz')
            log(i).condition = 'quiz';
            trialNum_train = trialNum_train + 1;
            log(i).trial = trialNum_train;
          elseif strcmp(logData{4}{i},'retention')
            log(i).condition = 'retention';
            trialNum_reten = trialNum_reten + 1;
            log(i).trial = trialNum_reten;
          end
          
          %log(i).condition = log(i).train_type;
          
        end
        
        % get the english item name
        log(i).item_eng = logData{7}{i};
        % find out where this item was shown on the screen
        log(i).xcoord_eng = str2num(logData{8}{i});
        log(i).ycoord_eng = str2num(logData{9}{i});
        
        % get the swahili item name
        log(i).item_swa = logData{10}{i};
        % find out where this item was shown on the screen
        log(i).xcoord_swa = str2num(logData{11}{i});
        log(i).ycoord_swa = str2num(logData{12}{i});
        
        % mark if it's a target (1) or buffer (0)
        log(i).rec_isTarg = str2double(logData{13}{i});
        
        % get the serial position
        log(i).serialpos = str2double(logData{14}{i});
        
      case {'TEST_TARGET','TEST_LURE'}
        % only do this if it's the swahili presentation
        if str2num(logData{9}{i}) < 0.5 && str2num(logData{9}{i}) > 0
          trialNum_recog = trialNum_recog + 1;
          
          log(i).condition = 'recogrecall';
          
          % get the swahili item name
          log(i).item_swa = logData{5}{i};
          % find out where this item was shown on the screen
          log(i).xcoord_swa = str2num(logData{8}{i});
          log(i).ycoord_swa = str2num(logData{9}{i});
          % mark if it was a target
          log(i).rec_isTarg = str2double(logData{6}{i});
          
          % get the serial position
          log(i).serialpos = str2double(logData{7}{i});
          % get the trial number within this list
          log(i).trial = trialNum_recog;
        else
          % to identify the useless (english) trials that get picked up
          % with this switch/case condition
          log(i).trial = -1;
        end
        
      case {'RECOG_RESP'}
        log(i).condition = log(i-1).condition;
        
        % put in the swahili item name
        log(i).item_swa = log(i-1).item_swa;
        % put in the target status
        log(i).rec_isTarg = log(i-1).rec_isTarg;
        % put in the serial position in which it was shown during study
        log(i).serialpos = log(i-1).serialpos;
        % put in the study location
        log(i).xcoord_swa = log(i-1).xcoord_swa;
        log(i).ycoord_swa = log(i-1).ycoord_swa;
        
        % put in the trial number within this list
        log(i).trial = log(i-1).trial;
        
        % get the response
        log(i).rec_resp = logData{5}{i};
        % was it correct?
        log(i).rec_correct = str2double(logData{8}{i});
        % get the reaction time
        log(i).rec_resp_rt = str2double(logData{6}{i});
        % if they didn't respond
        if log(i).rec_resp_rt == maxTestTime
          log(i).rec_resp = '';
          log(i).rec_correct = -1;
        end
        
        % get the english stimulus info, only if OLD
        if log(i).rec_isTarg && log(i).rec_correct
          log(i).item_eng = logData{5}{i+3};
          log(i).xcoord_eng = str2num(logData{8}{i+3});
          log(i).ycoord_eng = str2num(logData{9}{i+3});
        else
          % items marked NEW don't have english stimulus info logged
          log(i).item_eng = '';
          log(i).xcoord_eng = -1;
          log(i).ycoord_eng = -1;
        end
        
        % if it's an old item and they said "old", there's no new response;
        % if it's a new item and they said "old", there's no new response
        if (log(i).rec_isTarg && log(i).rec_correct) || (~log(i).rec_isTarg && ~log(i).rec_correct)
          log(i).new_correct = -1;
          log(i).new_resp = '';
          log(i).new_resp_rt = -1;
        end
        
        % put this info in test presentation
        log(i-1).rec_correct = log(i).rec_correct;
        log(i-1).rec_resp = log(i).rec_resp;
        log(i-1).rec_resp_rt = log(i).rec_resp_rt;
        
        log(i-1).item_eng = log(i).item_eng;
        log(i-1).xcoord_eng = log(i).xcoord_eng;
        log(i-1).ycoord_eng = log(i).ycoord_eng;
        
        log(i-1).new_correct = log(i).new_correct;
        log(i-1).new_resp = log(i).new_resp;
        log(i-1).new_resp_rt = log(i).new_resp_rt;
        
      case {'NEW_RESP'}
        % bring the RECOG_RESP info forward
        log(i).rec_correct = log(i-1).rec_correct;
        log(i).rec_resp = log(i-1).rec_resp;
        log(i).rec_resp_rt = log(i-1).rec_resp_rt;
        
        log(i).condition = log(i-2).condition;
        
        % put in the item name
        log(i).item_swa = log(i-2).item_swa;
        log(i).item_eng = log(i-2).item_eng;
        % put in the target status
        log(i).rec_isTarg = log(i-2).rec_isTarg;
        % put in the serial position in which it was shown
        log(i).serialpos = log(i-2).serialpos;
        % put in the study location
        log(i).xcoord_swa = log(i-2).xcoord_swa;
        log(i).ycoord_swa = log(i-2).ycoord_swa;
        log(i).xcoord_eng = log(i-2).xcoord_eng;
        log(i).ycoord_eng = log(i-2).ycoord_eng;
        % put in the trial number within this list
        log(i).trial = log(i-2).trial;
        
        % get the response
        log(i).new_resp = logData{5}{i};
        % get the reaction time
        log(i).new_resp_rt = str2double(logData{6}{i});
        % was it correct?
        log(i).new_correct = str2double(logData{8}{i});
        % if they didn't respond
        if log(i).new_resp_rt == maxTestTime
          log(i).new_resp = '';
          log(i).new_correct = -1;
        end
        
        % new responses don't get recognition recall opportunities
        log(i).rec_recall = -1;
        log(i).rec_recall_resp = '';
        
        % put info in recognition response
        log(i-1).new_resp = log(i).new_resp;
        log(i-1).new_resp_rt = log(i).new_resp_rt;
        log(i-1).new_correct = log(i).new_correct;
        log(i-1).rec_recall = log(i).rec_recall;
        log(i-1).rec_recall_resp = log(i).rec_recall_resp;
        % put info in test presentation
        log(i-2).new_resp = log(i).new_resp;
        log(i-2).new_resp_rt = log(i).new_resp_rt;
        log(i-2).new_correct = log(i).new_correct;
        log(i-2).rec_recall = log(i).rec_recall;
        log(i-2).rec_recall_resp = log(i).rec_recall_resp;
        
    end % switch event type
  end % for i = 1:length(log)
  
  % grab only the events we want from our log struct
  events_ses = filterStruct(log,'ismember(type,varargin{1})',{'STUDY_BUFFER','STUDY_TARGET','TEST_TARGET','TEST_LURE','RECOG_RESP','NEW_RESP'});
  % get rid of those trial == -1 events
  events_ses = filterStruct(events_ses,'trial > 0');
  
  % store session_0 and session_1 separately
  if sesNum == 0
    events_ses0 = events_ses;
  elseif sesNum == 1
    events_ses1 = events_ses;
  end
end % sessions

%%  read in annotation file

% go through session_0 events and add annotation info: quiz
% (train_eng_#.ann) and recognition/recall (recall_targets_eng_#.ann)

% go through session_1 events and add annotation info: retention
% (reten_eng_#.ann)

% initialize the counter to track the number of quizzes missing annotation
% files
noquiz = 0;

for ses = 1:length(sessions)
  session = sessions{ses};
  sesDir = fullfile(dataroot,subject,session);
  sesNum = str2double(strrep(session,'session_',''));
  
  if strcmp(session,'session_0')
    theseLists = unique(getStructField(events_ses0,'list'));
  elseif strcmp(session,'session_1')
    theseLists = unique(getStructField(events_ses1,'list'));
  end
  theseLists = theseLists - 1;
  
  for i = 1:length(theseLists)
    if strcmp(session,'session_0')
      % get the quiz annotation file
      annFile = fullfile(sesDir,[quizAnnFile,num2str(theseLists(i)),'.ann']);
      if exist(annFile,'file')
        fprintf('Reading quiz annotation file for list %d (%s%d.ann).\n',i,quizAnnFile,theseLists(i));
        fid = fopen(annFile);
        ann_contents = textscan(fid, '%f%d%s', 'Delimiter', '\t', 'CommentStyle', '#');
        fclose(fid);
        % adjust the word recall onset time
        ann_contents{1} = ann_contents{1} + ann_start_time_quiz(i);
        
        % add the quiz annotation events onto the end of events_ses0
        for w = 1:length(ann_contents{1})
          events_ses0(end+1).subject = subject;
          events_ses0(end).session = sesNum;
          events_ses0(end).mstime = round(ann_contents{1}(w));
          events_ses0(end).msoffset = 0;
          events_ses0(end).list = i;
          events_ses0(end).condition = 'quiz';
          events_ses0(end).type = 'REC_WORD';
          events_ses0(end).quiz_recall_resp = ann_contents{3}{w};
          % add in blank rec_resp so we can filter on this field later
          events_ses0(end).rec_resp = '';
        end
      else
        noquiz = noquiz + 1;
        if noquiz <= (numLists / 2)
          fprintf('No quiz annotation file for list %d (%s%d.ann).\n\tThis is ok, only half of the lists should have them (%d of %d).\n',i,quizAnnFile,theseLists(i),noquiz,numLists);
          %continue
        elseif noquiz > (numLists / 2)
          error('More than half of the lists are missing an annotation files.\n\tThis is not ok, there are too many missing files.');
        end
      end
      
      % get the recognition/recall annotation file
      annFile = fullfile(sesDir,[recogAnnFile,num2str(theseLists(i)),'.ann']);
      if exist(annFile,'file')
        fprintf('Reading recognition/recall annotation file for list %d (%s%d.ann).\n',i,recogAnnFile,theseLists(i));
        fid = fopen(annFile);
        ann_contents = textscan(fid, '%f%d%s', 'Delimiter', '\t', 'CommentStyle', '#');
        fclose(fid);
        % adjust the word recall onset time
        ann_contents{1} = ann_contents{1} + ann_start_time_recog(i);
        
        % add the recogrecall annotation events onto the end of events_ses0
        for w = 1:length(ann_contents{1})
          events_ses0(end+1).subject = subject;
          events_ses0(end).session = sesNum;
          events_ses0(end).mstime = round(ann_contents{1}(w));
          events_ses0(end).msoffset = 0;
          events_ses0(end).list = i;
          events_ses0(end).condition = 'recogrecall';
          events_ses0(end).type = 'REC_WORD';
          events_ses0(end).rec_recall_resp = ann_contents{3}{w};
          % add in blank rec_resp so we can filter on this field later
          events_ses0(end).rec_resp = '';
        end
      else
        error('No recognition/recall annotation file for list %d (%s%d).\n',i,recogAnnFile,theseLists(i));
      end
      
    elseif strcmp(session,'session_1')
      
      % get the retention test annotation file
      annFile = fullfile(sesDir,[retenAnnFile,num2str(theseLists(i)),'.ann']);
      if exist(annFile,'file')
        fprintf('Reading retention annotation file for list %d (%s%d.ann).\n',i,retenAnnFile,theseLists(i));
        fid = fopen(annFile);
        ann_contents = textscan(fid, '%f%d%s', 'Delimiter', '\t', 'CommentStyle', '#');
        fclose(fid);
        % adjust the word recall onset time
        ann_contents{1} = ann_contents{1} + ann_start_time_reten(i);
        
        % now add the retention annotation events onto the end of events_ses1
        for w = 1:length(ann_contents{1})
          events_ses1(end+1).subject = subject;
          events_ses1(end).session = sesNum;
          events_ses1(end).mstime = round(ann_contents{1}(w));
          events_ses1(end).msoffset = 0;
          events_ses1(end).list = i;
          events_ses1(end).condition = 'retention';
          events_ses1(end).type = 'REC_WORD';
          events_ses1(end).ret_recall_resp = ann_contents{3}{w};
        end
      else
        error('No retention annotation file for list %d (%s%d).\n',i,retenAnnFile,theseLists(i));
      end
      
    end % which session
  end % theseLists
  
  % sort events by the mstime field
  % http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/
  if strcmp(session,'session_0')
    ev0fields = fieldnames(events_ses0);
    ev0cell = struct2cell(events_ses0);
    sz = size(ev0cell);
    % Convert to a matrix
    ev0cell = reshape(ev0cell, sz(1), []); % Px(MxN)
    % Make each field a column
    ev0cell = ev0cell'; % (MxN)xP
    % Sort by field 'mstime'
    ev0cell = sortrows(ev0cell, find(ismember(ev0fields,'mstime')));
    % Put back into original cell array format
    ev0cell = reshape(ev0cell', sz);
    % Convert to Struct
    events_ses0 = cell2struct(ev0cell, ev0fields, 1);
  elseif strcmp(session,'session_1')
    ev1fields = fieldnames(events_ses1);
    ev1cell = struct2cell(events_ses1);
    sz = size(ev1cell);
    % Convert to a matrix
    ev1cell = reshape(ev1cell, sz(1), []); % Px(MxN)
    % Make each field a column
    ev1cell = ev1cell'; % (MxN)xP
    % Sort by field 'mstime'
    ev1cell = sortrows(ev1cell, find(ismember(ev0fields,'mstime')));
    % Put back into original cell array format
    ev1cell = reshape(ev1cell', sz);
    % Convert to Struct
    events_ses1 = cell2struct(ev1cell, ev1fields, 1);
  end
end % ses

%% info transfer

% first transfer the study information to the recall events

for i = 1:length(events_ses0)
  % if we find a study event, get the corresponding quiz or recogrecall
  % test event and move the info into the study event
  switch events_ses0(i).type
    case {'REC_WORD'}
      
      switch events_ses0(i).condition
        case {'quiz'}
          % quiz
          %
          % find the swahili item presented just prior to the recalled word
          timeThresh = maxSpeakDur;
          [quizCueEvents,qInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & mstime >= varargin{2} & mstime < varargin{3}',{'STUDY_TARGET'},(events_ses0(i).mstime - timeThresh),events_ses0(i).mstime);
          if ~isempty(quizCueEvents)
            if length(quizCueEvents) == 1
              events_ses0(i).trial = quizCueEvents.trial;
              events_ses0(i).train_type = quizCueEvents.train_type;
              events_ses0(i).item_swa = quizCueEvents.item_swa;
              events_ses0(i).xcoord_swa = quizCueEvents.xcoord_swa;
              events_ses0(i).ycoord_swa = quizCueEvents.ycoord_swa;
              events_ses0(i).item_eng = quizCueEvents.item_eng;
              events_ses0(i).xcoord_eng = quizCueEvents.xcoord_eng;
              events_ses0(i).ycoord_eng = quizCueEvents.ycoord_eng;
              events_ses0(i).serialpos = quizCueEvents.serialpos;
              events_ses0(i).rec_isTarg = quizCueEvents.rec_isTarg;
              
              if strcmp(events_ses0(i).quiz_recall_resp,quizCueEvents.item_eng)
                events_ses0(i).quiz_recall = 1;
              else
                events_ses0(i).quiz_recall = 0;
              end
              
              % put the recall info in the quiz restudy event
              events_ses0(qInd).quiz_recall = events_ses0(i).quiz_recall;
              events_ses0(qInd).quiz_recall_resp = events_ses0(i).quiz_recall_resp;
            else
              % need to figure out a solution if this ever arises
              fprintf('Too many quizCueEvents (n=%d).\n',length(quizCueEvents));
              keyboard
            end
          else
            % TODO: if there was an intrusion or something, make sure this
            % is the right thing to do
            events_ses0(i).trial = -1;
            events_ses0(i).train_type = '';
            events_ses0(i).item_swa = '';
            events_ses0(i).xcoord_swa = -1;
            events_ses0(i).ycoord_swa = -1;
            events_ses0(i).item_eng = '';
            events_ses0(i).xcoord_eng = -1;
            events_ses0(i).ycoord_eng = -1;
            events_ses0(i).serialpos = -1;
            events_ses0(i).rec_isTarg = -1;
            events_ses0(i).rec_correct = -1;
            events_ses0(i).rec_resp = '';
            events_ses0(i).rec_resp_rt = -1;
            events_ses0(i).new_correct = -1;
            events_ses0(i).new_resp = '';
            events_ses0(i).new_resp_rt = -1;
            events_ses0(i).quiz_recall = -1;
            events_ses0(i).rec_recall = -1;
            events_ses0(i).ret_recall = -1;
            %events_ses0(i).quiz_recall_resp = '';
            events_ses0(i).rec_recall_resp = '';
            events_ses0(i).ret_recall_resp = '';
          end
          
        case {'recogrecall'}
          % recogrecall
          %
          % find the swahili item presented just prior to the recalled word
          % (with a little extra time built in)
          timeThresh = maxSpeakDur + 2000;
          [recCueEvents,rInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & mstime >= varargin{2} & mstime < varargin{3} & ismember(rec_resp,varargin{4})',{'RECOG_RESP'},(events_ses0(i).mstime - timeThresh),events_ses0(i).mstime,{'OLD'});
          if ~isempty(recCueEvents)
            if length(recCueEvents) > 1
              % if there's more than one we need to figure out which was the
              % most recent and make sure it got an OLD response
              times = getStructField(recCueEvents,'mstime');
              [y,minInd] = min(events_ses0(i).mstime - times);
              if (recCueEvents(minInd).rec_isTarg && recCueEvents(minInd).rec_correct) || (~recCueEvents(minInd).rec_isTarg && ~recCueEvents(minInd).rec_correct)
                recCueEvents = recCueEvents(minInd);
                % and get rid of the indices for the other events
                these_rInd = find(rInd);
                notMin = ones(size(these_rInd));
                notMin(minInd) = 0;
                rInd(these_rInd(find(these_rInd .* notMin))) = 0;
              else
                fprintf('The most recent event did not get an OLD response. WTF no way this is the match for the recalled word.'\n);
                keyboard
              end
            end
            if length(recCueEvents) == 1
              events_ses0(i).trial = recCueEvents.trial;
              %events_ses0(i).train_type = recCueEvents.train_type;
              events_ses0(i).item_swa = recCueEvents.item_swa;
              events_ses0(i).xcoord_swa = recCueEvents.xcoord_swa;
              events_ses0(i).ycoord_swa = recCueEvents.ycoord_swa;
              events_ses0(i).item_eng = recCueEvents.item_eng;
              events_ses0(i).xcoord_eng = recCueEvents.xcoord_eng;
              events_ses0(i).ycoord_eng = recCueEvents.ycoord_eng;
              events_ses0(i).serialpos = recCueEvents.serialpos;
              events_ses0(i).rec_isTarg = recCueEvents.rec_isTarg;
              events_ses0(i).rec_correct = recCueEvents.rec_correct;
              events_ses0(i).rec_resp = recCueEvents.rec_resp;
              events_ses0(i).rec_resp_rt = recCueEvents.rec_resp_rt;
              events_ses0(i).new_correct = recCueEvents.new_correct;
              events_ses0(i).new_resp = recCueEvents.new_resp;
              events_ses0(i).new_resp_rt = recCueEvents.new_resp_rt;
              
              if strcmp(events_ses0(i).rec_recall_resp,recCueEvents.item_eng)
                events_ses0(i).rec_recall = 1;
              else
                events_ses0(i).rec_recall = 0;
              end
              
              % put the recall info in the recognition event
              events_ses0(rInd).rec_recall = events_ses0(i).rec_recall;
              events_ses0(rInd).rec_recall_resp = events_ses0(i).rec_recall_resp;
              % make sure we get the right test target event
              if strcmp(events_ses0(rInd).item_swa,events_ses0(find(rInd)-1).item_swa)
                events_ses0(find(rInd)-1).rec_recall = events_ses0(i).rec_recall;
                events_ses0(find(rInd)-1).rec_recall_resp = events_ses0(i).rec_recall_resp;
              else
                timeThresh = maxSpeakDur + 3000;
                [testTargEvents,tInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & mstime >= varargin{2} & mstime < varargin{3} & ismember(rec_resp,varargin{4})',{'TEST_TARGET'},(events_ses0(i).mstime - timeThresh),events_ses0(i).mstime,{'OLD'});
                if length(testTargEvents) == 1
                  events_ses0(tInd).rec_recall = events_ses0(i).rec_recall;
                  events_ses0(tInd).rec_recall_resp = events_ses0(i).rec_recall_resp;
                else
                  keyboard
                end
              end
            else
              fprintf('Too many recCueEvents (n=%d).\n',length(recCueEvents));
              keyboard
            end
          else
            % TODO: if there was an intrusion or something, make sure this
            % is the right thing to do
            events_ses0(i).trial = -1;
            events_ses0(i).train_type = '';
            events_ses0(i).item_swa = '';
            events_ses0(i).xcoord_swa = -1;
            events_ses0(i).ycoord_swa = -1;
            events_ses0(i).item_eng = '';
            events_ses0(i).xcoord_eng = -1;
            events_ses0(i).ycoord_eng = -1;
            events_ses0(i).serialpos = -1;
            events_ses0(i).rec_isTarg = -1;
            events_ses0(i).rec_correct = -1;
            events_ses0(i).rec_resp = '';
            events_ses0(i).rec_resp_rt = -1;
            events_ses0(i).new_correct = -1;
            events_ses0(i).new_resp = '';
            events_ses0(i).new_resp_rt = -1;
            events_ses0(i).quiz_recall = -1;
            events_ses0(i).rec_recall = -1;
            events_ses0(i).ret_recall = -1;
            events_ses0(i).quiz_recall_resp = '';
            %events_ses0(i).rec_recall_resp = '';
            events_ses0(i).ret_recall_resp = '';
          end
      end
  end
end

%% add in the english translation for targets if we need it

% because of how the experiment currently logs data, item_eng is blank when
% targets are labeled "new" during recogrecall. we want that data, so get
% it.

for i = 1:length(events_ses0)
  switch events_ses0(i).type
    case {'TEST_TARGET','RECOG_RESP','NEW_RESP'}
      if events_ses0(i).rec_isTarg && isempty(events_ses0(i).item_eng)
        [studyEvent] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(condition,varargin{3}) & list == varargin{4}',{'STUDY_TARGET'},{events_ses0(i).item_swa},{'study'},events_ses0(i).list);
        events_ses0(i).item_eng = studyEvent.item_eng;
      end
  end
end

%% add in retention test accuracy to events_ses1 events

for i = 1:length(events_ses1)
  % if we find a study event, get the corresponding test event
  % and move the test info into the study event
  switch events_ses1(i).type
    case {'STUDY_TARGET'}
      % initialize all test items with blank/no recall
      events_ses1(i).ret_recall = 0;
      events_ses1(i).ret_recall_resp = '';
      
    case {'REC_WORD'}
      % retention test
      %
      % find the swahili item presented just prior to the recalled word
      timeThresh = maxSpeakDur;
      [retenCueEvents,rInd] = filterStruct(events_ses1,'ismember(type,varargin{1}) & mstime >= varargin{2} & mstime < varargin{3}',{'STUDY_TARGET'},(events_ses1(i).mstime - timeThresh),events_ses1(i).mstime);
      if ~isempty(retenCueEvents)
        if length(retenCueEvents) == 1
          events_ses1(i).trial = retenCueEvents.trial;
          events_ses1(i).train_type = retenCueEvents.train_type;
          events_ses1(i).item_swa = retenCueEvents.item_swa;
          events_ses1(i).xcoord_swa = retenCueEvents.xcoord_swa;
          events_ses1(i).ycoord_swa = retenCueEvents.ycoord_swa;
          events_ses1(i).item_eng = retenCueEvents.item_eng;
          events_ses1(i).xcoord_eng = retenCueEvents.xcoord_eng;
          events_ses1(i).ycoord_eng = retenCueEvents.ycoord_eng;
          events_ses1(i).serialpos = retenCueEvents.serialpos;
          events_ses1(i).rec_isTarg = retenCueEvents.rec_isTarg;
          
          if strcmp(events_ses1(i).ret_recall_resp,retenCueEvents.item_eng)
            events_ses1(i).ret_recall = 1;
          else
            events_ses1(i).ret_recall = 0;
          end
          
          % put the recall info in the retention event
          events_ses1(rInd).ret_recall = events_ses1(i).ret_recall;
          events_ses1(rInd).ret_recall_resp = events_ses1(i).ret_recall_resp;
        else
          fprintf('Too many retenCueEvents (n=%d).\n',length(retenCueEvents));
          keyboard
        end
      else
        % I think this is an (extra-list?) intrusion. Or maybe a very rare
        % instance where they make an ELI for a word that will later be on
        % a study list. Whatever the case, it doesn't have this info.
        events_ses1(i).trial = -1;
        events_ses1(i).train_type = '';
        events_ses1(i).item_swa = '';
        events_ses1(i).xcoord_swa = -1;
        events_ses1(i).ycoord_swa = -1;
        events_ses1(i).item_eng = '';
        events_ses1(i).xcoord_eng = -1;
        events_ses1(i).ycoord_eng = -1;
        events_ses1(i).serialpos = -1;
        events_ses1(i).rec_isTarg = -1;
        
        % more info?
        events_ses1(i).rec_correct = -1;
        events_ses1(i).rec_resp = '';
        events_ses1(i).rec_resp_rt = -1;
        events_ses1(i).new_correct = -1;
        events_ses1(i).new_resp = '';
        events_ses1(i).new_resp_rt = -1;
        events_ses1(i).quiz_recall = -1;
        events_ses1(i).rec_recall = -1;
        events_ses1(i).ret_recall = -1;
        events_ses1(i).quiz_recall_resp = '';
        events_ses1(i).rec_recall_resp = '';
        % % don't reset the ret_recall_resp field
        % events_ses1(i).ret_recall_resp = '';
      end
  end
end

%% more info transfer

% put the info from recognition/recall test into the study/train events,
% and the info from study/train into the recognition/recall test events
for i = 1:length(events_ses0)
  % if we find a study event, get the corresponding test event
  % and move the test info into the study event
  switch events_ses0(i).type
    case {'STUDY_TARGET'}
      [testEvent,tInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3})',{'TEST_TARGET'},{events_ses0(i).item_swa},{events_ses0(i).item_eng});
      if length(testEvent) > 1
        error('Found multiple events when searching for study target %s!',events_ses0(i).item_swa);
      elseif isempty(testEvent)
        % this would only happen when they didn't finish a session
        warning('Subject %s possibly did not finish session_0',subject);
        keyboard
        
        % study presentation - test info
        events_ses0(i).rec_correct = -1;
        events_ses0(i).rec_resp = '';
        events_ses0(i).rec_resp_rt = -1;
        events_ses0(i).new_correct = -1;
        events_ses0(i).new_resp = '';
        events_ses0(i).new_resp_rt = -1;
        
        events_ses0(i).rec_recall = -1;
        events_ses0(i).rec_recall_resp = '';
        
      elseif length(testEvent) == 1
        % put it in the study presentation
        events_ses0(i).rec_correct = testEvent.rec_correct;
        events_ses0(i).rec_resp = testEvent.rec_resp;
        events_ses0(i).rec_resp_rt = testEvent.rec_resp_rt;
        events_ses0(i).new_correct = testEvent.new_correct;
        events_ses0(i).new_resp = testEvent.new_resp;
        events_ses0(i).new_resp_rt = testEvent.new_resp_rt;
        
        if ~isempty(testEvent.rec_recall)
          events_ses0(i).rec_recall = testEvent.rec_recall;
          events_ses0(i).rec_recall_resp = testEvent.rec_recall_resp;
        elseif isempty(testEvent.rec_recall) && testEvent.rec_isTarg && testEvent.rec_correct
          % if it's empty but was correctly called old then there was a
          % chance for recall
          events_ses0(i).rec_recall = 0;
          events_ses0(i).rec_recall_resp = '';
          events_ses0(tInd).rec_recall = 0;
          events_ses0(tInd).rec_recall_resp = '';
        elseif isempty(testEvent.rec_recall) && testEvent.rec_isTarg && ~testEvent.rec_correct
          % if it's empty but was correctly called new then there was no
          % chance for recall, so it's -1
          events_ses0(i).rec_recall = -1;
          events_ses0(i).rec_recall_resp = '';
          events_ses0(tInd).rec_recall = -1;
          events_ses0(tInd).rec_recall_resp = '';
        else
          events_ses0(i).rec_recall = -1;
          events_ses0(i).rec_recall_resp = '';
          events_ses0(tInd).rec_recall = -1;
          events_ses0(tInd).rec_recall_resp = '';
        end
      end
      
      % add quiz info to study event
      [quizEvent,qInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & ismember(condition,varargin{4})',{'STUDY_TARGET'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},{'quiz'});
      if length(quizEvent) > 1
        error('Found multiple events when searching for study target %s!',events_ses0(i).item_swa);
      elseif isempty(quizEvent)
        % % this might only happen when they didn't finish a session
        % warning('Subject %s possibly did not finish session_0',subject);
        % keyboard
        
        % this could also happen when the training was a restudy
        
        events_ses0(i).quiz_recall = -1;
        events_ses0(i).quiz_recall_resp = '';
        
      elseif length(quizEvent) == 1
        if ~isempty(quizEvent.quiz_recall)
          events_ses0(i).quiz_recall = quizEvent.quiz_recall;
          events_ses0(i).quiz_recall_resp = quizEvent.quiz_recall_resp;
        elseif isempty(quizEvent.quiz_recall)
          events_ses0(i).quiz_recall = 0;
          events_ses0(i).quiz_recall_resp = '';
          events_ses0(qInd).quiz_recall = 0;
          events_ses0(qInd).quiz_recall_resp = '';
        else
          events_ses0(i).quiz_recall = -1;
          events_ses0(i).quiz_recall_resp = '';
          events_ses0(qInd).quiz_recall = -1;
          events_ses0(qInd).quiz_recall_resp = '';
        end
      end
      
      % add retention info to study event
      [retenCueEvent] = filterStruct(events_ses1,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & ismember(condition,varargin{4})',{'STUDY_TARGET'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},{'retention'});
      if length(retenCueEvent) > 1
        error('Found multiple events when searching for study target %s!',events_ses1(i).item_swa);
      elseif isempty(retenCueEvent)
        % % this might only happen when they didn't finish a session
        % warning('Subject %s possibly did not finish session_0',subject);
        % keyboard
        
        % this could also happen when the training was a restudy
        
        events_ses1(i).ret_recall = -1;
        events_ses1(i).ret_recall_resp = '';
        
      elseif length(retenCueEvent) == 1
        [quizEvent,qInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & ismember(condition,varargin{4})',{'STUDY_TARGET'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},{'quiz','restudy'});
        if length(quizEvent) ~= 1
          keyboard
        end
        if ~isempty(retenCueEvent.ret_recall)
          % study event
          events_ses0(i).ret_recall = retenCueEvent.ret_recall;
          events_ses0(i).ret_recall_resp = retenCueEvent.ret_recall_resp;
          % quiz event
          events_ses0(qInd).ret_recall = retenCueEvent.ret_recall;
          events_ses0(qInd).ret_recall_resp = retenCueEvent.ret_recall_resp;
        elseif isempty(retenCueEvent.ret_recall)
          events_ses0(i).ret_recall = 0;
          events_ses0(i).ret_recall_resp = '';
          events_ses0(qInd).ret_recall = 0;
          events_ses0(qInd).ret_recall_resp = '';
        else
          events_ses0(i).ret_recall = -1;
          events_ses0(i).ret_recall_resp = '';
          events_ses0(qInd).ret_recall = -1;
          events_ses0(qInd).ret_recall_resp = '';
        end
      end
      
      % if it was a buffer, put in -1s
    case {'STUDY_BUFFER'}
      % study presentation - test info
      events_ses0(i).rec_correct = -1;
      events_ses0(i).rec_resp = '';
      events_ses0(i).rec_resp_rt = -1;
      events_ses0(i).new_correct = -1;
      events_ses0(i).new_resp = '';
      events_ses0(i).new_resp_rt = -1;
      events_ses0(i).quiz_recall = -1;
      events_ses0(i).rec_recall = -1;
      events_ses0(i).ret_recall = -1;
      events_ses0(i).quiz_recall_resp = '';
      events_ses0(i).rec_recall_resp = '';
      events_ses0(i).ret_recall_resp = '';
      
    case {'TEST_TARGET'}
      % if we find a test event, get the corresponding study event
      % and move the study info into the test event
      [studyEvent] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & ismember(condition,varargin{4})',{'STUDY_TARGET'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},{'study'});
      if isempty(studyEvent)
        keyboard
      elseif length(studyEvent) > 1
        error('Found multiple events when searching for study target %s!',events_ses0(i).item_swa);
      elseif length(studyEvent) == 1
        % there will always be more than 1 studyEvent because of the
        % study-training pairings, but we narrow it down to only 1 event
        % with condition==study (also we use both english and swahili to
        % make sure we have the right event)
        
        events_ses0(i).train_type = studyEvent.train_type;
        events_ses0(i).quiz_recall = studyEvent.quiz_recall;
        events_ses0(i).quiz_recall_resp = studyEvent.quiz_recall_resp;
        events_ses0(i).ret_recall = studyEvent.ret_recall;
        events_ses0(i).ret_recall_resp = studyEvent.ret_recall_resp;
        
        % add it to the RECOG_RESP
        [rrEvents,rrInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & rec_isTarg == 1 & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & list == varargin{4}',{'RECOG_RESP'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},events_ses0(i).list);
        if length(rrEvents) == 1
          events_ses0(rrInd).train_type = studyEvent.train_type;
          events_ses0(rrInd).quiz_recall = studyEvent.quiz_recall;
          events_ses0(rrInd).quiz_recall_resp = studyEvent.quiz_recall_resp;
          % TODO: do we need rec_recall? yes I think so
          events_ses0(rrInd).rec_recall = studyEvent.rec_recall;
          events_ses0(rrInd).rec_recall_resp = studyEvent.rec_recall_resp;
          events_ses0(rrInd).ret_recall = studyEvent.ret_recall;
          events_ses0(rrInd).ret_recall_resp = studyEvent.ret_recall_resp;
        else
          keyboard
        end
        
        if (~events_ses0(i).rec_isTarg && events_ses0(i).rec_correct) || (events_ses0(i).rec_isTarg && ~events_ses0(i).rec_correct)
          % if it was a "new" response
          [newRespEvents,nrInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & list == varargin{3}',{'NEW_RESP'},{events_ses0(i).item_swa},events_ses0(i).list);
          if length(newRespEvents) == 1
            events_ses0(nrInd).train_type = studyEvent.train_type;
            events_ses0(nrInd).quiz_recall = studyEvent.quiz_recall;
            events_ses0(nrInd).quiz_recall_resp = studyEvent.quiz_recall_resp;
            events_ses0(nrInd).ret_recall = studyEvent.ret_recall;
            events_ses0(nrInd).ret_recall_resp = studyEvent.ret_recall_resp;
          else
            keyboard
          end
        end
      end
      
    case {'TEST_LURE'}
      % if we find a lure test event, put in -1s
      events_ses0(i).train_type = '';
      events_ses0(i).quiz_recall = -1;
      events_ses0(i).quiz_recall_resp = '';
      events_ses0(i).ret_recall = -1;
      events_ses0(i).ret_recall_resp = '';
      events_ses0(i).rec_recall = -1;
      events_ses0(i).rec_recall_resp = '';
      
      % put this info in the RECOG_RESP
      [rRespEvents,rrInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & list == varargin{4}',{'RECOG_RESP'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},events_ses0(i).list);
      if length(rRespEvents) == 1
        events_ses0(rrInd).train_type = events_ses0(i).train_type;
        events_ses0(rrInd).quiz_recall = events_ses0(i).quiz_recall;
        events_ses0(rrInd).quiz_recall_resp = events_ses0(i).quiz_recall_resp;
        events_ses0(rrInd).rec_recall = events_ses0(i).ret_recall;
        events_ses0(rrInd).rec_recall_resp = events_ses0(i).ret_recall_resp;
        events_ses0(rrInd).ret_recall = events_ses0(i).ret_recall;
        events_ses0(rrInd).ret_recall_resp = events_ses0(i).ret_recall_resp;
      else
        keyboard
      end
      
      if events_ses0(i).rec_correct == 1
        % if it's a new item and they got it correct, then there will be a
        % third log event ("NEW_RESP"); if they got it wrong (called it
        % "old") then they will be prompted for a vocal response and this
        % won't get accessed
        [newRespEvents,nrInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & list == varargin{4}',{'NEW_RESP'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},events_ses0(i).list);
        
        if length(newRespEvents) == 1
          events_ses0(nrInd).train_type = events_ses0(i).train_type;
          events_ses0(nrInd).quiz_recall = events_ses0(i).quiz_recall;
          events_ses0(nrInd).quiz_recall_resp = events_ses0(i).quiz_recall_resp;
          events_ses0(nrInd).rec_recall = events_ses0(i).ret_recall;
          events_ses0(nrInd).rec_recall_resp = events_ses0(i).ret_recall_resp;
          events_ses0(nrInd).ret_recall = events_ses0(i).ret_recall;
          events_ses0(nrInd).ret_recall_resp = events_ses0(i).ret_recall_resp;
        else
          keyboard
        end
      end
      
    case {'REC_WORD'}
      switch events_ses0(i).condition
        case {'quiz'}
          timeThresh = maxSpeakDur;
          [testTargEvents,tInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & mstime >= varargin{2} & mstime < varargin{3}',{'STUDY_TARGET'},(events_ses0(i).mstime - timeThresh),events_ses0(i).mstime);
        case {'recogrecall'}
          %timeThresh = maxSpeakDur + 5000;
          %[testTargEvents,tInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & mstime >= varargin{2} & mstime < varargin{3} & ismember(item_swa,varargin{4}) & list = varargin{5}',{'TEST_TARGET','TEST_LURE'},{events_ses0(i).item_swa},events_ses0(i).list);
          [testTargEvents,tInd] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & list == varargin{4}',{'TEST_TARGET','TEST_LURE'},{events_ses0(i).item_swa},{events_ses0(i).item_eng},events_ses0(i).list);
      end
      
      if length(testTargEvents) == 1
        events_ses0(i).train_type = events_ses0(tInd).train_type;
        events_ses0(i).rec_correct = events_ses0(tInd).rec_correct;
        events_ses0(i).rec_resp = events_ses0(tInd).rec_resp;
        events_ses0(i).rec_resp_rt = events_ses0(tInd).rec_resp_rt;
        events_ses0(i).new_correct = events_ses0(tInd).new_correct;
        events_ses0(i).new_resp = events_ses0(tInd).new_resp;
        events_ses0(i).new_resp_rt = events_ses0(tInd).new_resp_rt;
        events_ses0(i).quiz_recall = events_ses0(tInd).quiz_recall;
        events_ses0(i).quiz_recall_resp = events_ses0(tInd).quiz_recall_resp;
        events_ses0(i).rec_recall = events_ses0(tInd).rec_recall;
        events_ses0(i).rec_recall_resp = events_ses0(tInd).rec_recall_resp;
        events_ses0(i).ret_recall = events_ses0(tInd).ret_recall;
        events_ses0(i).ret_recall_resp = events_ses0(tInd).ret_recall_resp;
      elseif isempty(testTargEvents) && strcmp(events_ses0(i).rec_recall_resp,vocal)
        % TODO: just some kind of intrusion or random vocalization?
        fprintf('Might have found an intrusion or vocalization at %s %s list %d, event %d (not associated with a list?).\n',events_ses0(i).subject,events_ses0(i).condition,events_ses0(i).list,i);
      else
        keyboard
      end
  end
end

%% transfer info regarding whether word was recalled during retention

for i = 1:length(events_ses1)
  % put the study and recogrecall info in the retention test events
  switch events_ses1(i).type
    case {'STUDY_TARGET','REC_WORD'}
      [studyEvent] = filterStruct(events_ses0,'ismember(type,varargin{1}) & ismember(item_swa,varargin{2}) & ismember(item_eng,varargin{3}) & ismember(condition,varargin{4})',{'STUDY_TARGET'},{events_ses1(i).item_swa},{events_ses1(i).item_eng},{'study'});
      if ~isempty(studyEvent)
        events_ses1(i).train_type = studyEvent.train_type;
        events_ses1(i).rec_correct = studyEvent.rec_correct;
        events_ses1(i).rec_resp = studyEvent.rec_resp;
        events_ses1(i).rec_resp_rt = studyEvent.rec_resp_rt;
        events_ses1(i).new_correct = studyEvent.new_correct;
        events_ses1(i).new_resp = studyEvent.new_resp;
        events_ses1(i).new_resp_rt = studyEvent.new_resp_rt;
        events_ses1(i).quiz_recall = studyEvent.quiz_recall;
        events_ses1(i).quiz_recall_resp = studyEvent.quiz_recall_resp;
        events_ses1(i).rec_recall = studyEvent.rec_recall;
        events_ses1(i).rec_recall_resp = studyEvent.rec_recall_resp;
      end
  end
end

%% put fields in an orderly manner
fnames = {'subject', 'session', 'mstime', 'msoffset', 'list', 'condition', 'trial', 'type', 'train_type', 'item_swa', 'xcoord_swa', 'ycoord_swa', 'item_eng', 'xcoord_eng', 'ycoord_eng', 'serialpos', 'rec_isTarg', 'rec_correct', 'rec_resp', 'rec_resp_rt', 'new_correct', 'new_resp', 'new_resp_rt', 'quiz_recall', 'quiz_recall_resp', 'rec_recall', 'rec_recall_resp', 'ret_recall', 'ret_recall_resp'};
events_ses0 = orderfields(events_ses0,fnames);
events_ses1 = orderfields(events_ses1,fnames);

%% matlab junk

%#ok<*ST2NM>

%% more debugging stuff

% % debug
% fn0 = fieldnames(events_ses0);
% for i = 1:length(fn0)
%   fprintf('%s: %d\n',fn0{i},length(getStructField(events_ses0,fn0{i})));
%   if length(getStructField(events_ses0,fn0{i})) ~= length(events_ses0)
%     fprintf('ERROR: %s has a length of %d\n',fn0{i},length(getStructField(events_ses0,fn0{i})));
%   end
% end
% fn1 = fieldnames(events_ses1);
% for i = 1:length(fn1)
%   fprintf('%s: %d\n',fn1{i},length(getStructField(events_ses1,fn1{i})));
%   if length(getStructField(events_ses1,fn1{i})) ~= length(events_ses1)
%     fprintf('ERROR: %s has a length of %d\n',fn1{i},length(getStructField(events_ses1,fn1{i})));
%   end
% end
%
% % debug: check whether we can grab all the fields and look for empties
% fn0 = fieldnames(events_ses0);
% for i = 1:length(fn0)
%   thisev0 = getStructField(events_ses0,fn0{i});
%   for j = 1:length(thisev0)
%     if iscell(thisev0(j))
%       tester = thisev0{j};
%     else
%       tester = thisev0(j);
%     end
%     if ~ischar(tester) && isempty(tester)
%       fprintf('uh oh, empty found! events_ses0(%d).%s\n',j,fn0{i});
%       keyboard
%     end
%   end
% end
%
% fn1 = fieldnames(events_ses1);
% for i = 1:length(fn1)
%   thisev1 = getStructField(events_ses1,fn1{i});
%   for j = 1:length(thisev1)
%     if iscell(thisev1(j))
%       tester = thisev1{j};
%     else
%       tester = thisev1(j);
%     end
%     if ~ischar(tester) && isempty(tester)
%       fprintf('uh oh, empty found! events_ses1(%d).%s\n',j,fn1{i});
%       keyboard
%     end
%   end
% end
