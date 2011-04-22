function events = sosi_createEvents(dataroot,subject,session)
% function events = sosi_createEvents(dataroot,subject,session)
%
% create event struct for SOSI
%
% struct fields:
%   subject
%   session
%   mstime
%   msoffset
%   list
%   trial
%   type
%   item
%   serialpos
%   xcoord
%   ycoord
%   rec_isTarg
%   rec_correct: didn't say "new" to an old item, but can get source wrong
%   src_isTarg
%   src_resp
%   src_rt
%   src_correct
%   rkn_resp
%   rkn_rt
%   rkn_correct
%
% For source response, correct is calculated relative to a "SIZE" response
% (i.e., "hit" if SIZE response for SIZE source, and "false alarm" if SIZE
% response for LIFE source). This could be switched around to be relative to
% LIFE.
%
% For source response, correct is calculated relative to a "right
% side" (R) response (i.e., "hit" if R response for R source, and
% "false alarm" if R response for L source). This could be switched
% around to be relative to L.
%

%dataroot = '/Volumes/curranlab/Data/SOSI/eeg/behavioral';
%subject = 'SOSI001';
%session = 'session_0';

sesDir = fullfile(dataroot,subject,session);

logfile = fullfile(sesDir,'session.log');

fid = fopen(logfile);
logdata = textscan(fid,'%n%n%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);

%[l1 l2 l3 l4 l5 l6 l7 l8] = textread(logfile,'%d%d%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);

sesNum = str2double(strrep(session,'session_',''));
log = struct('subject',subject,'session',sesNum,'mstime',num2cell(logdata{1}),'msoffset',num2cell(logdata{2}),'trial',[],'type',logdata{3});

% constants
%numLists = 4;
maxTestTime = 30000;

% initialize
listNum = NaN;

for i = 1:length(log)
  log(i).mstime = log(i).mstime;
  log(i).msoffset = log(i).msoffset;
  
  % set the list number and test type
  if strcmp('TRIAL',logdata{3}{i}) % TRIAL = listNum
    listNum = str2double(logdata{4}{i});
    % if we're starting a new list, set trialNum to zero
    trialNum = 0;
  end
  log(i).list = listNum;
  
  switch log(i).type
    
    case {'STUDY_TARGET','STUDY_BUFFER'}
      % get the image name
      log(i).item = logdata{5}{i};
      % mark if it's a target (1) or buffer (0)
      log(i).rec_isTarg = str2double(logdata{6}{i});
      % find out where this item was shown on the screen
      log(i).xcoord = str2num(logdata{7}{i});
      log(i).ycoord = str2num(logdata{8}{i});
      % get the serial position
      log(i).serialpos = str2double(logdata{9}{i});
      % get the trial number within this list
      log(i).trial = log(i).serialpos;
      
    case {'TEST_TARGET','TEST_LURE','TARG_PRES','LURE_PRES'}
      trialNum = trialNum + 1;
      % get the image name
      log(i).item = logdata{5}{i};
      % mark if it was a target
      log(i).rec_isTarg = str2double(logdata{6}{i});
      % get the serial position
      log(i).serialpos = str2double(logdata{7}{i});
      % find out where this item was shown on the screen
      %log(i).xcoord = str2num(logdata{8}{i});
      %log(i).ycoord = str2num(logdata{9}{i});
      % get the trial number within this list
      log(i).trial = trialNum;
      
    case 'SOURCE_RESP'
      % put in the item name
      log(i).item = log(i-1).item;
      % put in the target status
      log(i).rec_isTarg = log(i-1).rec_isTarg;
      % put in the serial position in which it was shown during study
      log(i).serialpos = log(i-1).serialpos;
      % put in the study location
      %log(i).xcoord = log(i-1).xcoord;
      %log(i).ycoord = log(i-1).ycoord;
      % put in the trial number within this list
      log(i).trial = log(i-1).trial;

      % get the response
      log(i).src_resp = logdata{5}{i};
      % was it correct?
      log(i).src_correct = str2double(logdata{8}{i});
      % get the reaction time
      log(i).src_rt = str2double(logdata{6}{i});
      if log(i).src_rt == maxTestTime
        log(i).src_resp = '';
        log(i).src_correct = -1;
      end
      
      % determine if they correctly recognized it as either old
      % (independent of source accuracy) or new
      if log(i).rec_isTarg == 1 && ~strcmp(log(i).src_resp,'NEW')
        log(i).rec_correct = 1;
      elseif log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,'NEW')
        log(i).rec_correct = 0;
      elseif log(i).rec_isTarg == 0 && strcmp(log(i).src_resp,'NEW')
        log(i).rec_correct = 1;
      elseif log(i).rec_isTarg == 0 && ~strcmp(log(i).src_resp,'NEW')
        log(i).rec_correct = 0;
      else
        log(i).rec_correct = -1;
      end
      
      % put info in test presentation
      log(i-1).rec_correct = log(i).rec_correct;
      log(i-1).src_resp = log(i).src_resp;
      log(i-1).src_rt = log(i).src_rt;
      log(i-1).src_correct = log(i).src_correct;
      
    case {'RK_RESP','NEW_RESP'}
      % bring the SOURCE_RESP info forward
      %log(i).src_isTarg = log(i-1).src_isTarg;
      log(i).src_resp = log(i-1).src_resp;
      log(i).src_rt = log(i-1).src_rt;
      log(i).src_correct = log(i-1).src_correct;
      log(i).rec_correct = log(i-1).rec_correct;
      
      % put in the item name
      log(i).item = log(i-2).item;
      % put in the target status
      log(i).rec_isTarg = log(i-2).rec_isTarg;
      % put in the serial position in which it was shown
      log(i).serialpos = log(i-2).serialpos;
      % put in the study location
      %log(i).xcoord = log(i-2).xcoord;
      %log(i).ycoord = log(i-2).ycoord;
      % put in the trial number within this list
      log(i).trial = log(i-2).trial;
      
      % get the response
      log(i).rkn_resp = logdata{5}{i};
      % get the reaction time
      log(i).rkn_rt = str2double(logdata{6}{i});
      % was it correct?
      log(i).rkn_correct = str2double(logdata{8}{i});
      if log(i).rkn_rt == maxTestTime
        log(i).rkn_resp = '';
        log(i).rkn_correct = -1;
      end
      
      % put info in recognition response
      log(i-1).rkn_resp = log(i).rkn_resp;
      log(i-1).rkn_rt = log(i).rkn_rt;
      log(i-1).rkn_correct = log(i).rkn_correct;
      % put info in test presentation
      log(i-2).rkn_resp = log(i).rkn_resp;
      log(i-2).rkn_rt = log(i).rkn_rt;
      log(i-2).rkn_correct = log(i).rkn_correct;
      
  end % switch event type
end % for i = 1:length(log)

% grab only the events we want from our log struct
events = filterStruct(log,'ismember(type,varargin{1})',{'STUDY_BUFFER','STUDY_TARGET','STUDY_RESP','TARG_PRES','LURE_PRES','SOURCE_RESP','RK_RESP','NEW_RESP','TEST_TARGET','TEST_LURE'});

% put the info from test into the study events, and the info from
% study into the test events
for i = 1:length(events)
  % if we find a study event, get the corresponding test event
  % and move the test info into the study event
  switch events(i).type
    case {'STUDY_TARGET'}
      testEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'TEST_TARGET'},{events(i).item});
      if length(testEvent) > 1
        error('Found multiple events when searching for study target %s!',events(i).item);
      elseif isempty(testEvent)
        
        % when would this happen? only when they didn't finish a session
        %keyboard
        
        % study presentation - test info
        events(i).rec_correct = -1;
        events(i).src_isTarg = -1;
        events(i).src_resp = -1;
        events(i).src_rt = -1;
        events(i).src_correct = -1;
        events(i).rkn_resp = -1;
        events(i).rkn_rt = -1;
        events(i).rkn_correct = -1;
      elseif length(testEvent) == 1
        % put it in the study presentation
        events(i).rec_correct = testEvent.rec_correct;
        events(i).src_resp = testEvent.src_resp;
        events(i).src_correct = testEvent.src_correct;
        events(i).src_rt = testEvent.src_rt;
        events(i).rkn_resp = testEvent.rkn_resp;
        events(i).rkn_correct = testEvent.rkn_correct;
        events(i).rkn_rt = testEvent.rkn_rt;
        
        % put in the target status where all RIGHT SIDE sources are targets
        % and all LEFT SIDE  sources are lures; -1 if otherwise; this is
        % described at the top of the function
        if events(i).rec_isTarg == 1 && (events(i).xcoord > 0.5 && events(i).xcoord <= 1)
          % if it was on right then it's a target
          events(i).src_isTarg = 1;
        elseif events(i).rec_isTarg == 1 && (events(i).xcoord < 0.5 && events(i).xcoord >= 0)
          % if it was on left then it's a lure
          events(i).src_isTarg = 0;
        end
      end
      
      % if it was a buffer, put in -1s
    case {'STUDY_BUFFER'}
      % study presentation - test info
      events(i).rec_correct = -1;
      events(i).src_isTarg = -1;
      events(i).src_resp = '';
      events(i).src_rt = -1;
      events(i).src_correct = -1;
      events(i).rkn_resp = '';
      events(i).rkn_rt = -1;
      events(i).rkn_correct = -1;
      
      % if we find a test event, get the corresponding study event
      % and move the study info into the test event
    case {'TEST_TARGET','TARG_PRES'}
      studyEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_TARGET'},{events(i).item});
      if length(studyEvent) > 1
        error('Found multiple events when searching for study target %s!',events(i).item);
      elseif isempty(studyEvent)
        keyboard
      elseif length(studyEvent) == 1
        events(i).src_isTarg = studyEvent.src_isTarg;
        events(i+1).src_isTarg = studyEvent.src_isTarg;
        events(i+2).src_isTarg = studyEvent.src_isTarg;
        
        events(i).xcoord = studyEvent.xcoord;
        events(i+1).xcoord = studyEvent.xcoord;
        events(i+2).xcoord = studyEvent.xcoord;
        events(i).ycoord = studyEvent.ycoord;
        events(i+1).ycoord = studyEvent.ycoord;
        events(i+2).ycoord = studyEvent.ycoord;
      end
      
      % if we find a lure test event, put in -1s
    case {'TEST_LURE','LURE_PRES'}
      % test presentation
      events(i).src_isTarg = -1;
      % recognition response
      events(i+1).src_isTarg = -1;
      % source response
      events(i+2).src_isTarg = -1;
      
      events(i).xcoord = -1;
      events(i+1).xcoord = -1;
      events(i+2).xcoord = -1;
      events(i).ycoord = -1;
      events(i+1).ycoord = -1;
      events(i+2).ycoord = -1;
  end
end

% put fields in an orderly manner
events = orderfields(events,{'subject','session','mstime','msoffset','list','trial','type','item','serialpos','xcoord','ycoord','rec_isTarg','rec_correct','src_isTarg','src_resp','src_rt','src_correct','rkn_resp','rkn_rt','rkn_correct'});

% % debug
% fn = fieldnames(events);
% for i = 1:length(fn)
%   fprintf('%s: %d\n',fn{i},length(getStructField(events,fn{i})));
%   if length(getStructField(events,fn{i})) ~= length(events)
%     fprintf('ERROR: %s has a length of %d\n',fn{i},length(getStructField(events,fn{i})));
%   end
% end

