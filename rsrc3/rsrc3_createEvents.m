function [study_events,test_events] = rsrc3_createEvents(dataroot,subject,session)
% function [study_events,test_events] = rsrc3_events(dataroot,subject,session)
%
% create event struct for RSRC
%
% event struct fields:
%   subject
%   session
%   mstime
%   msoffset
%   list
%   trial
%   type
%   item
%   xcoord
%   ycoord
%   serialpos
%   study_resp
%   study_correct
%   study_rt
%   study_tooslow
%   rec_isTarg
%   rec_resp
%   rec_conf
%   rec_correct
%   rec_rt
%   rec_tooslow
%   src_isTarg
%   src_resp
%   src_correct
%   src_rt
%   src_tooslow
%
% For source response, correct is calculated relative to a "right
% side" (R) response (i.e., "hit" if R response for R source, and
% "false alarm" if R response for L source). This could be switched
% around to be relative to L.
%

% dataroot = '/Volumes/curranlab/Data/RSRC3/eeg/behavioral';
% subject = 'RSRC3001';
% session = 1;

if session == 1
  sessions = {'session_0','session_1'};
elseif session == 2
  sessions = {'session_2','session_3'};
end

if str2double(subject(end)) > 0 && str2double(subject(end)) <= 5
  % confidence numbers are [1 2 3 4 5 6]: RS, RO, DF, MF, MU, DU (old to new)
  confRange = [1 2 3 4 5 6];
  confRange_rec_str = {'REMEMBER_SOURCE','REMEMBER_OTHER','DEFINITELY_FAMILIAR','MAYBE_FAMILIAR','MAYBE_UNFAMILIAR','DEFINITELY_UNFAMILIAR'};
elseif str2double(subject(end)) == 0 || str2double(subject(end)) > 5
  % confidence numbers are [6 5 4 3 2 1]: DU, MU, MF, DF, RO, RS (new to old)
  confRange = [6 5 4 3 2 1];
  confRange_rec_str = {'DEFINITELY_UNFAMILIAR','MAYBE_UNFAMILIAR','MAYBE_FAMILIAR','DEFINITELY_FAMILIAR','REMEMBER_OTHER','REMEMBER_SOURCE'};
end

% constants
maxTestTime = 10000;

for ses = 1:length(sessions)
  
  sesDir = fullfile(dataroot,subject,sessions{ses});
  
  logfile = fullfile(sesDir,'session.log');
  
  fid = fopen(logfile);
  logdata = textscan(fid,'%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
  fclose(fid);
  
  sesNum = str2double(sessions{ses}(end));
  log = struct('subject',subject,'session',sesNum,'mstime',logdata{1},'msoffset',logdata{2},'trial',[],'type',logdata{3});
  
  listNum = NaN;
  
  for i = 1:length(log)
    % convert to numbers
    log(i).mstime = str2num(log(i).mstime);
    log(i).msoffset = str2num(log(i).msoffset);
    
    % set the trial number
    if strcmp('TRIAL',logdata{3}{i})
      listNum = str2num(logdata{4}{i});
      % if we're starting a new list, set trialNum to zero
      trialNum = 0;
    end
    log(i).list = listNum;
    
    switch log(i).type
      
      case 'STUDY_BUFFER'
        % get the image name
        log(i).item = logdata{5}{i};
        % mark that it's a buffer (and not a target)
        log(i).rec_isTarg = str2num(logdata{6}{i});
        % find out where this item was shown on the screen
        log(i).xcoord = str2num(logdata{7}{i});
        log(i).ycoord = str2num(logdata{8}{i});
        % get the serial position
        log(i).serialpos = str2num(logdata{9}{i});
        % get the trial number within this list
        log(i).trial = log(i).serialpos;
        % put in -1 where we won't have info
        log(i).rec_resp = '';
        log(i).rec_conf = -1;
        log(i).rec_correct = -1;
        log(i).rec_rt = -1;
        log(i).rec_tooslow = -1;
        log(i).src_isTarg = -1;
        log(i).src_resp = '';
        log(i).src_correct = -1;
        log(i).src_rt = -1;
        log(i).src_tooslow = -1;
        
      case 'STUDY_TARGET'
        % get the image name
        log(i).item = logdata{5}{i};
        % mark that it's a target (and not a buffer)
        log(i).rec_isTarg = str2num(logdata{6}{i});
        % find out where this item was shown on the screen
        log(i).xcoord = str2num(logdata{7}{i});
        log(i).ycoord = str2num(logdata{8}{i});
        % get the serial position
        log(i).serialpos = str2num(logdata{9}{i});
        % get the trial number within this list
        log(i).trial = log(i).serialpos;
        % we'll mark the source target status later because this is not
        % dependent on study target status
        
        % put in -1 where we won't have info
        log(i).rec_resp = '';
        log(i).rec_conf = -1;
        log(i).rec_correct = -1;
        log(i).rec_rt = -1;
        log(i).rec_tooslow = -1;
        log(i).src_isTarg = -1;
        log(i).src_resp = '';
        log(i).src_correct = -1;
        log(i).src_rt = -1;
        log(i).src_tooslow = -1;
        
      case 'STUDY_RESP'
        % get the study response and reaction time; find out if they were
        % correct or were too slow
        log(i).study_resp = logdata{5}{i};
        log(i).study_correct = str2num(logdata{8}{i});
        log(i).study_rt = str2num(logdata{6}{i});
        log(i).study_tooslow = str2num(logdata{9}{i});
        
        % put it in the study presentation
        log(i-1).study_resp = log(i).study_resp;
        log(i-1).study_correct = log(i).study_correct;
        log(i-1).study_rt = log(i).study_rt;
        log(i-1).study_tooslow = log(i).study_tooslow;
        
        % get it from the study presentation
        log(i).item = log(i-1).item;
        log(i).rec_isTarg = log(i-1).rec_isTarg;
        log(i).xcoord = log(i-1).xcoord;
        log(i).ycoord = log(i-1).ycoord;
        log(i).serialpos = log(i-1).serialpos;
        log(i).trial = log(i-1).trial;
        
        % put in -1 where we won't have info
        log(i).rec_resp = '';
        log(i).rec_conf = -1;
        log(i).rec_correct = -1;
        log(i).rec_rt = -1;
        log(i).rec_tooslow = -1;
        log(i).src_isTarg = -1;
        log(i).src_resp = '';
        log(i).src_correct = -1;
        log(i).src_rt = -1;
        log(i).src_tooslow = -1;
        
      case {'TARG_PRES','LURE_PRES'}
        trialNum = trialNum + 1;
        % get the image name
        log(i).item = logdata{5}{i};
        % mark that it was a target
        log(i).rec_isTarg = str2num(logdata{6}{i});
        % find out where this item was shown on the screen during study
        log(i).xcoord = str2num(logdata{7}{i});
        log(i).ycoord = str2num(logdata{8}{i});
        % get the study serial position
        log(i).serialpos = str2num(logdata{9}{i});
        % get the trial number within this list
        log(i).trial = trialNum;
        
        
      case 'RK_RESP'
        % put in the item name
        log(i).item = log(i-1).item;
        % put in the target status
        log(i).rec_isTarg = log(i-1).rec_isTarg;
        % put in the study coordinates
        log(i).xcoord = log(i-1).xcoord;
        log(i).ycoord = log(i-1).ycoord;
        % put in the serial position in which it was shown during study
        log(i).serialpos = log(i-1).serialpos;
        % and the test list trial number
        log(i).trial = log(i-1).trial;
        
        % get the response
        log(i).rec_resp = logdata{5}{i};
        % put it in the test presentation
        log(i-1).rec_resp = log(i).rec_resp;
        
        % get the confidence level response and convert it to a number; sure
        % old to sure new
        if strcmp(logdata{5}{i},confRange_rec_str{1})
          log(i).rec_conf = confRange(1);
        elseif strcmp(logdata{5}{i},confRange_rec_str{2})
          log(i).rec_conf = confRange(2);
        elseif strcmp(logdata{5}{i},confRange_rec_str{3})
          log(i).rec_conf = confRange(3);
        elseif strcmp(logdata{5}{i},confRange_rec_str{4})
          log(i).rec_conf = confRange(4);
        elseif strcmp(logdata{5}{i},confRange_rec_str{5})
          log(i).rec_conf = confRange(5);
        elseif strcmp(logdata{5}{i},confRange_rec_str{6})
          log(i).rec_conf = confRange(6);
        end
        % put info in test presentation
        log(i-1).rec_conf = log(i).rec_conf;
        
        % was it correct?
        log(i).rec_correct = str2num(logdata{8}{i});
        % put info in test presentation
        log(i-1).rec_correct = log(i).rec_correct;
        
        % get the reaction time
        log(i).rec_rt = str2num(logdata{6}{i});
        % put info in test presentation
        log(i-1).rec_rt = log(i).rec_rt;
        if log(i).rec_rt == maxTestTime
          log(i).rec_resp = '';
          log(i).rec_correct = -1;
        end
        
        % were they slow to respond?
        log(i).rec_tooslow = str2num(logdata{9}{i});
        % put info in test presentation
        log(i-1).rec_tooslow = log(i).rec_tooslow;
        
      case 'SOURCE_RESP'
        % bring the RK_RESP info forward
        log(i).rec_isTarg = log(i-1).rec_isTarg;
        log(i).rec_resp = log(i-1).rec_resp;
        log(i).rec_conf = log(i-1).rec_conf;
        log(i).rec_correct = log(i-1).rec_correct;
        log(i).rec_rt = log(i-1).rec_rt;
        log(i).rec_tooslow = log(i-1).rec_tooslow;
        
        % put in the item name
        log(i).item = log(i-2).item;
        
        % put in the coordinates
        log(i).xcoord = log(i-2).xcoord;
        log(i).ycoord = log(i-2).ycoord;
        
        % put in the serial position in which it was shown
        log(i).serialpos = log(i-2).serialpos;
        % and the test list trial number
        log(i).trial = log(i-2).trial;
        
        % get the response
        log(i).src_resp = logdata{5}{i};
        % put info in recognition response
        log(i-1).src_resp = log(i).src_resp;
        % put it in the test presentation
        log(i-2).src_resp = log(i).src_resp;
        
        % put in the target status where all R sources are targets and all L
        % sources are lures; -1 if otherwise; this is described at the top
        % of the function
        if log(i).xcoord > 0.5 && log(i).xcoord <= 1 % if it's a R source then
          % it's a target
          log(i).src_isTarg = 1;
          % put info in recognition response
          log(i-1).src_isTarg = log(i).src_isTarg;
          % put info in test presentation
          log(i-2).src_isTarg = log(i).src_isTarg;
        elseif log(i).xcoord < 0.5 && log(i).xcoord >= 0 % if it's a L source
          % then it's a lure
          log(i).src_isTarg = 0;
          % put info in recognition response
          log(i-1).src_isTarg = log(i).src_isTarg;
          % put info in test presentation
          log(i-2).src_isTarg = log(i).src_isTarg;
        else % otherwise it was a real lure, and we mark that as -1
          log(i).src_isTarg = -1;
          % put info in recognition response
          log(i-1).src_isTarg = log(i).src_isTarg;
          % put info in test presentation
          log(i-2).src_isTarg = log(i).src_isTarg;
        end
        
        % mark whether they got the source correct
        % otherwise mark it as -1;
        if log(i).src_isTarg == 1 || log(i).src_isTarg == 0
          log(i).src_correct = str2num(logdata{8}{i});
        else
          log(i).src_correct = -1;
        end
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct;
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct;
        
        % get the reaction time
        log(i).src_rt = str2num(logdata{6}{i});
        % put info in recognition response
        log(i-1).src_rt = log(i).src_rt;
        % put info in test presentation
        log(i-2).src_rt = log(i).src_rt;
        
        % were they slow to respond?
        log(i).src_tooslow = str2num(logdata{9}{i});
        % put info in recognition presentation
        log(i-1).src_tooslow = log(i).src_tooslow;
        % put info in test presentation
        log(i-2).src_tooslow = log(i).src_tooslow;
        
        % put in a -1 where we won't have info
        log(i).study_resp = '';
        log(i).study_correct = -1;
        log(i).study_rt = -1;
        log(i).study_tooslow = -1;
        log(i-1).study_resp = '';
        log(i-1).study_correct = -1;
        log(i-1).study_rt = -1;
        log(i-1).study_tooslow = -1;
        log(i-2).study_resp = '';
        log(i-2).study_correct = -1;
        log(i-2).study_rt = -1;
        log(i-2).study_tooslow = -1;
    end % switch event type
  end % for i = 1:length(log)
  if ses == 1
    study_log = filterStruct(log,'ismember(type,varargin{1})',{'STUDY_BUFFER','STUDY_TARGET','STUDY_RESP'});
    study_events = orderfields(study_log,{'subject','session','mstime','msoffset','list','trial','type','item','xcoord','ycoord','serialpos','study_resp','study_correct','study_rt','study_tooslow','rec_isTarg','rec_resp','rec_conf','rec_correct','rec_rt','rec_tooslow','src_isTarg','src_resp','src_correct','src_rt','src_tooslow'});
  elseif ses == 2
    test_log = filterStruct(log,'ismember(type,varargin{1})',{'TARG_PRES','LURE_PRES','RK_RESP','SOURCE_RESP'});
    test_events = orderfields(test_log,{'subject','session','mstime','msoffset','list','trial','type','item','xcoord','ycoord','serialpos','study_resp','study_correct','study_rt','study_tooslow','rec_isTarg','rec_resp','rec_conf','rec_correct','rec_rt','rec_tooslow','src_isTarg','src_resp','src_correct','src_rt','src_tooslow'});
  end
end

% % collapse the repeated lists
% for i = 1:length(study_log)
%   switch study_log(i).type
%     case {'STUDY_BUFFER'}
%       % find the two events
%       studyEvent = filterStruct(study_log,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_BUFFER'},{study_log(i).item});
%       study_log_new = studyEvent(1);
%       for j = 2:length(studyEvent)
%         study_log_new.mstime = cat(2,study_log_new.mstime,study_log_new(j).mstime);
%         study_log_new.list = cat(2,study_log_new.list,study_log_new(j).list);
%         study_log_new.trial = cat(2,study_log_new.trial,study_log_new(j).trial);
%         study_log_new.serialpos = cat(2,study_log_new.serialpos,study_log_new(j).serialpos);
%         study_log_new.study_resp = cat(2,study_log_new.study_resp,study_log_new(j).study_resp);
%         study_log_new.study_correct = cat(2,study_log_new.study_correct,study_log_new(j).study_correct);
%         study_log_new.study_rt = cat(2,study_log_new.study_rt,study_log_new(j).study_rt);
%         study_log_new.study_tooslow = cat(2,study_log_new.study_tooslow,study_log_new(j).study_tooslow);
%       end
%     case {'STUDY_TARGET'}
%       studyEvent = filterStruct(study_log,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_TARGET'},{study_log(i).item});
%     case {'STUDY_RESP'}
%       studyEvent = filterStruct(study_log,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_RESP'},{study_log(i).item});
%   end
% end

% % combine study and test sessions
% combined_log = [study_log;test_log];
% 
% % grab only the events we want from our log struct
% %combined_log = filterStruct(combined_log,'ismember(type,varargin{1})',{'STUDY_BUFFER','STUDY_TARGET','STUDY_RESP','TARG_PRES','LURE_PRES','RK_RESP','SOURCE_RESP'});
% % put fields in an orderly manner
% events =
% orderfields(combined_log,{'subject','session','mstime','msoffset','list','trial','type','item','xcoord','ycoord','serialpos','study_resp','study_correct','study_rt','study_tooslow','rec_isTarg','rec_resp','rec_conf','rec_correct','rec_rt','rec_tooslow','src_isTarg','src_resp','src_correct','src_rt','src_tooslow'});

for i = 1:length(study_events)
  switch study_events(i).type
    case {'STUDY_TARGET'}
      % put the info from test into the study events
      testEvent = filterStruct(test_events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'TARG_PRES'},{study_events(i).item});
      if length(testEvent) > 1
        error('Found multiple test events when searching for study target %s!',study_events(i).item);
      else
        % put it in the study response
        study_events(i).rec_isTarg = testEvent.rec_isTarg;
        study_events(i).rec_resp = testEvent.rec_resp;
        study_events(i).rec_conf = testEvent.rec_conf;
        study_events(i).rec_correct = testEvent.rec_correct;
        study_events(i).rec_rt = testEvent.rec_rt;
        study_events(i).rec_tooslow = testEvent.rec_tooslow;
        study_events(i).src_isTarg = testEvent.src_isTarg;
        study_events(i).src_resp = testEvent.src_resp;
        study_events(i).src_correct = testEvent.src_correct;
        study_events(i).src_rt = testEvent.src_rt;
        study_events(i).src_tooslow = testEvent.src_tooslow;
        % and put it in the study presentation
        study_events(i+1).rec_isTarg = testEvent.rec_isTarg;
        study_events(i+1).rec_resp = testEvent.rec_resp;
        study_events(i+1).rec_conf = testEvent.rec_conf;
        study_events(i+1).rec_correct = testEvent.rec_correct;
        study_events(i+1).rec_rt = testEvent.rec_rt;
        study_events(i+1).rec_tooslow = testEvent.rec_tooslow;
        study_events(i+1).src_isTarg = testEvent.src_isTarg;
        study_events(i+1).src_resp = testEvent.src_resp;
        study_events(i+1).src_correct = testEvent.src_correct;
        study_events(i+1).src_rt = testEvent.src_rt;
        study_events(i+1).src_tooslow = testEvent.src_tooslow;
      end
  end
end

for i = 1:length(test_events)
  switch test_events(i).type
    case {'TARG_PRES'}
      % put the info from study into the test events
      studyEvent = filterStruct(study_events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_TARGET'},{test_events(i).item});
      if length(studyEvent) > 2 % this is 2 because the lists repeat
        error('Found multiple study events when searching for test target %s!',test_events(i).item);
      else
        % put it in the study presentation
        test_events(i).serialpos = getStructField(studyEvent,'serialpos');
        test_events(i).study_resp = getStructField(studyEvent,'study_resp');
        test_events(i).study_correct = getStructField(studyEvent,'study_correct');
        test_events(i).study_rt = getStructField(studyEvent,'study_rt');
        test_events(i).study_tooslow = getStructField(studyEvent,'study_tooslow');
        % and the source response
        test_events(i+1).serialpos = getStructField(studyEvent,'serialpos');
        test_events(i+1).study_resp = getStructField(studyEvent,'study_resp');
        test_events(i+1).study_correct = getStructField(studyEvent,'study_correct');
        test_events(i+1).study_rt = getStructField(studyEvent,'study_rt');
        test_events(i+1).study_tooslow = getStructField(studyEvent,'study_tooslow');
        % and the recognition response
        test_events(i+2).serialpos = getStructField(studyEvent,'serialpos');
        test_events(i+2).study_resp = getStructField(studyEvent,'study_resp');
        test_events(i+2).study_correct = getStructField(studyEvent,'study_correct');
        test_events(i+2).study_rt = getStructField(studyEvent,'study_rt');
        test_events(i+2).study_tooslow = getStructField(studyEvent,'study_tooslow');
      end
  end
end
