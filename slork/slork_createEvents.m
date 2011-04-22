function events = slork_createEvents(dataroot,subject,session)
% function events = slork_createEvents(dataroot,subject,session)
%
% create event struct for SLORK
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
%   pool
%   pool_correct
%   xcoord
%   ycoord
%   serialpos
%   study_cond
%   study_resp
%   study_rt
%   rec_isTarg
%   rec_correct: didn't say "new" to an old item, but can get source wrong
%   testType
%   src_isTarg
%   src_resp
%   src_rt
%   src_correct
%   rkn_resp
%   rkn_rt
%   rkn_correct: always == 1
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

%dataroot = '/Users/matt/data/SLORK/behavioral';
%subject = 'SLORK003';
%session = 'session_0';

sesDir = fullfile(dataroot,subject,session);

logfile = fullfile(sesDir,'session.log');

fid = fopen(logfile);
logdata = textscan(fid,'%n%n%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);

%[l1 l2 l3 l4 l5 l6 l7 l8] = textread(logfile,'%d%d%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);

sesNum = str2num(strrep(session,'session_',''));
log = struct('subject',subject,'session',sesNum,'mstime',num2cell(logdata{1}),'msoffset',num2cell(logdata{2}),'trial',[],'type',logdata{3});

% constants
numLists = 4;
maxStudyTime = 20000;
maxTestTime = 30000;

% initialize
listNum = NaN;
testTypes = cell(1,numLists);
testType = NaN;

% get the test types
for i = 1:length(log)
  if strcmp('TESTBLOCKS',logdata{3}{i})
    testTypes{1} = logdata{4}{i};
    testTypes{2} = logdata{5}{i};
    testTypes{3} = logdata{6}{i};
    testTypes{4} = logdata{7}{i};
    break
  end
end

for i = 1:length(log)
  log(i).mstime = log(i).mstime;
  log(i).msoffset = log(i).msoffset;
  
  % set the list number and test type
  if strcmp('TRIAL',logdata{3}{i}) % TRIAL = listNum
    listNum = str2num(logdata{4}{i});
    testType = testTypes{listNum};
    % if we're starting a new list, set trialNum to zero
    trialNum = 0;
  end
  log(i).list = listNum;
  log(i).testType = testType;
  
  switch log(i).type
    
    case {'STUDY_TARGET','STUDY_BUFFER'}
      % get the image name
      log(i).item = logdata{5}{i};
      % mark if it's a target (1) or buffer (0)
      log(i).rec_isTarg = str2num(logdata{6}{i});
      % find out where this item was shown on the screen
      log(i).xcoord = str2num(logdata{7}{i});
      log(i).ycoord = str2num(logdata{8}{i});
      % get the serial position
      log(i).serialpos = str2num(logdata{9}{i});
      % get the trial number within this list
      log(i).trial = log(i).serialpos;
      % get the study condition (judgment type)
      log(i).study_cond = logdata{10}{i};
      
    case 'STUDY_RESP'
      log(i).study_resp = logdata{5}{i};
      if strcmp(log(i).study_resp,'NONE')
        log(i).study_resp = '';
      end
      log(i).study_rt = str2num(logdata{6}{i});
      % % if they hit the max time, they didn't answer; bug in exp (3Mar2010)
      % if log(i).study_rt == maxStudyTime
      %   log(i).study_resp = '';
      % end
      % get info from presentation
      log(i).item = log(i-1).item;
      log(i).rec_isTarg = log(i-1).rec_isTarg;
      log(i).xcoord = log(i-1).xcoord;
      log(i).ycoord = log(i-1).ycoord;
      log(i).serialpos = log(i-1).serialpos;
      log(i).trial = log(i-1).trial;
      log(i).study_cond = log(i-1).study_cond;
      % put info in the presentation event
      log(i-1).study_resp = log(i).study_resp;
      log(i-1).study_rt = log(i).study_rt;
      
      %    case 'TEST_START'
      %     testType = logdata{4}{i};
      
    case {'TARG_PRES','LURE_PRES'}
      trialNum = trialNum + 1;
      % get the image name
      log(i).item = logdata{5}{i};
      % mark if it was a target
      log(i).rec_isTarg = str2num(logdata{6}{i});
      % find out where this item was shown on the screen
      log(i).xcoord = str2num(logdata{7}{i});
      log(i).ycoord = str2num(logdata{8}{i});
      % get the serial position
      log(i).serialpos = str2num(logdata{9}{i});
      % get the trial number within this list
      log(i).trial = trialNum;
      
    case 'SOURCE_RESP'
      % put in the item name
      log(i).item = log(i-1).item;
      % put in the coordinates
      log(i).xcoord = log(i-1).xcoord;
      log(i).ycoord = log(i-1).ycoord;
      % put in the serial position in which it was shown during study
      log(i).serialpos = log(i-1).serialpos;
      % put in the target status
      log(i).rec_isTarg = log(i-1).rec_isTarg;
      % put in the trial number within this list
      log(i).trial = log(i-1).trial;

      % get the response
      log(i).src_resp = logdata{5}{i};
      % was it correct?
      log(i).src_correct = str2num(logdata{8}{i});
      % get the reaction time
      log(i).src_rt = str2num(logdata{6}{i});
      if log(i).src_rt == maxTestTime
        log(i).src_resp = '';
        log(i).src_correct = -1;
      end
      
      % determine if they correctly recognized it as either old
      % (independent of source accuracy) or new
      switch log(i).testType
        case {'task'}
          if log(i).rec_isTarg == 1 && (strcmp(log(i).src_resp,'SIZE') || strcmp(log(i).src_resp,'LIFE'))
            log(i).rec_correct = 1;
          elseif log(i).rec_isTarg == 0 && (strcmp(log(i).src_resp,'SIZE') || strcmp(log(i).src_resp,'LIFE'))
            log(i).rec_correct = 0;
          elseif log(i).rec_isTarg == 0 && strcmp(log(i).src_resp,'NEW')
            log(i).rec_correct = 1;
          elseif log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,'NEW')
            log(i).rec_correct = 0;
          else
            log(i).rec_correct = -1;
          end
        case {'side'}
          if log(i).rec_isTarg == 1 && (strcmp(log(i).src_resp,'LEFT') || strcmp(log(i).src_resp,'RIGHT'))
            log(i).rec_correct = 1;
          elseif log(i).rec_isTarg == 0 && (strcmp(log(i).src_resp,'LEFT') || strcmp(log(i).src_resp,'RIGHT'))
            log(i).rec_correct = 0;
          elseif log(i).rec_isTarg == 0 && strcmp(log(i).src_resp,'NEW')
            log(i).rec_correct = 1;
          elseif log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,'NEW')
            log(i).rec_correct = 0;
          else
            log(i).rec_correct = -1;
          end
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
      % put in the coordinates
      log(i).xcoord = log(i-2).xcoord;
      log(i).ycoord = log(i-2).ycoord;
      % put in the serial position in which it was shown
      log(i).serialpos = log(i-2).serialpos;
      % put in the target status
      log(i).rec_isTarg = log(i-2).rec_isTarg;
      % put in the trial number within this list
      log(i).trial = log(i-2).trial;
      
      % get the reaction time
      log(i).rkn_rt = str2num(logdata{6}{i});
      % was it correct?
      log(i).rkn_correct = str2num(logdata{8}{i});
      % get the response
      log(i).rkn_resp = logdata{5}{i};
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

% search the video log for the names of the images, the times they
% were presented, and grab the pool name from the path
vidlogfile = fullfile(sesDir,'video.vidlog');
fid = fopen(vidlogfile);
vidlogdata = textscan(fid,'%n%n%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);
image_ms = [];
image_pool = {};
image_name = {};
count = 0;
for i = 1:length(vidlogdata{7})
  if strcmp(vidlogdata{7}(i),'IMAGE') && ~isempty(strfind(cell2mat(vidlogdata{8}(i)),'object'))
    count = count + 1;
    % save the ms presentation time for this image
    image_ms = [image_ms;vidlogdata{1}(i)];
    % save the image name
    [pathstr,image_name{count}] = fileparts(cell2mat(vidlogdata{8}(i)));
    % save the image pool
    [pathstr,image_pool{count}] = fileparts(pathstr);
  end
end

% grab only the events we want from our log struct
events = filterStruct(log,'ismember(type,varargin{1})',{'STUDY_BUFFER','STUDY_TARGET','STUDY_RESP','TARG_PRES','LURE_PRES','SOURCE_RESP','RK_RESP','NEW_RESP'});

% get rid of practice images
expStartInd = find(image_ms == events(1).mstime);
image_ms = image_ms(expStartInd:end);
image_pool = image_pool(expStartInd:end);
image_name = image_name(expStartInd:end);

for i = 1:length(image_ms)
  for j = 1:length(events)
    
    % insert the name of the image pool into the field 'pool'
    if strcmp(events(j).item,image_name(i))
      events(j).pool = image_pool{i};
    end
    
    switch events(j).type
      case {'STUDY_TARGET','STUDY_BUFFER'}
        
        % find out if they got the pool correct (by our classification)
        separator = strfind(events(j).pool,'_');
        lifepool = events(j).pool(1:separator-1);
        sizepool = events(j).pool(separator+1:end);
        
        if strcmp(lifepool,'living') && strcmp(events(j).study_resp,'LIVING')
          events(j).pool_correct = 1;
        elseif strcmp(lifepool,'nonliving') && strcmp(events(j).study_resp,'NONLIVING')
          events(j).pool_correct = 1;
        elseif strcmp(sizepool,'small') && strcmp(events(j).study_resp,'SMALLER')
          events(j).pool_correct = 1;
        elseif strcmp(sizepool,'large') && strcmp(events(j).study_resp,'BIGGER')
          events(j).pool_correct = 1;
        elseif isempty(events(j).study_resp)
          events(j).pool_correct = -1;
        else
          events(j).pool_correct = 0;
        end
    end
  end
end

% put the info from test into the study events, and the info from
% study into the test events
for i = 1:length(events)
  % if we find a study event, get the corresponding test event
  % and move the test info into the study event
  switch events(i).type
    case {'STUDY_TARGET'}
      testEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'TARG_PRES'},{events(i).item});
      if length(testEvent) > 1
        error('Found multiple events when searching for study target %s!',events(i).item);
      else
        % put it in the study presentation
        events(i).rec_correct = testEvent.rec_correct;
        events(i).src_resp = testEvent.src_resp;
        events(i).src_correct = testEvent.src_correct;
        events(i).src_rt = testEvent.src_rt;
        events(i).rkn_resp = testEvent.rkn_resp;
        events(i).rkn_correct = testEvent.rkn_correct;
        events(i).rkn_rt = testEvent.rkn_rt;
        % and the study response
        events(i+1).rec_correct = testEvent.rec_correct;
        events(i+1).src_resp = testEvent.src_resp;
        events(i+1).src_correct = testEvent.src_correct;
        events(i+1).src_rt = testEvent.src_rt;
        events(i+1).rkn_resp = testEvent.rkn_resp;
        events(i+1).rkn_correct = testEvent.rkn_correct;
        events(i+1).rkn_rt = testEvent.rkn_rt;
        
        % study response for if the study classification was correct
        events(i+1).pool_correct = events(i).pool_correct;
        
        % put in the target status where all SIZE and RIGHT SIDE sources are
        % targets and all LIFE and LEFT SIDE sources are lures; -1 if
        % otherwise; this is described at the top of the function
        if strcmp(events(i).testType,'task')
          if events(i).rec_isTarg == 1 && strcmp(events(i).study_cond,'SIZE')
            % if it's Size source then it's a target
            events(i).src_isTarg = 1;
            events(i+1).src_isTarg = 1;
          elseif events(i).rec_isTarg == 1 && strcmp(events(i).study_cond,'LIFE')
            % if it's Life source then it's a lure
            events(i).src_isTarg = 0;
            events(i+1).src_isTarg = 0;
          end
        elseif strcmp(events(i).testType,'side')
          if events(i).rec_isTarg == 1 && (events(i).xcoord > 0.5 && events(i).xcoord <= 1)
            % if it's Right source then it's a target
            events(i).src_isTarg = 1;
            events(i+1).src_isTarg = 1;
          elseif events(i).rec_isTarg == 1 && (events(i).xcoord < 0.5 && events(i).xcoord >= 0)
            % if it's Left source then it's a lure
            events(i).src_isTarg = 0;
            events(i+1).src_isTarg = 0;
          end
        end
      end
      
      % if it was a buffer, put in -1s
    case {'STUDY_BUFFER'}
      % study presentation
      events(i).rec_correct = -1;
      events(i).src_isTarg = -1;
      events(i).src_resp = -1;
      events(i).src_rt = -1;
      events(i).src_correct = -1;
      events(i).rkn_resp = -1;
      events(i).rkn_rt = -1;
      events(i).rkn_correct = -1;
      % study response
      events(i+1).rec_correct = -1;
      events(i+1).src_isTarg = -1;
      events(i+1).src_resp = -1;
      events(i+1).src_rt = -1;
      events(i+1).src_correct = -1;
      events(i+1).rkn_resp = -1;
      events(i+1).rkn_rt = -1;
      events(i+1).rkn_correct = -1;
      
      % study response for if the study classification was correct
      events(i+1).pool_correct = events(i).pool_correct;
      
      % if we find a test event, get the corresponding study event
      % and move the study info into the test event
    case {'TARG_PRES'}
      studyEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_TARGET'},{events(i).item});
      if length(studyEvent) > 1
        error('Found multiple events when searching for study target %s!',events(i).item);
      else
        % put it in the study presentation
        events(i).study_resp = studyEvent.study_resp;
        events(i).study_rt = studyEvent.study_rt;
        events(i).study_cond = studyEvent.study_cond;
        events(i).pool_correct = studyEvent.pool_correct;
        events(i).src_isTarg = studyEvent.src_isTarg;
        % and the source response
        events(i+1).study_resp = studyEvent.study_resp;
        events(i+1).study_rt = studyEvent.study_rt;
        events(i+1).study_cond = studyEvent.study_cond;
        events(i+1).pool_correct = studyEvent.pool_correct;
        events(i+1).src_isTarg = studyEvent.src_isTarg;
        % and the remem/know/new response
        events(i+2).study_resp = studyEvent.study_resp;
        events(i+2).study_rt = studyEvent.study_rt;
        events(i+2).study_cond = studyEvent.study_cond;
        events(i+2).pool_correct = studyEvent.pool_correct;
        events(i+2).src_isTarg = studyEvent.src_isTarg;
      end
      
      % if we find a lure test event, put in -1s
    case {'LURE_PRES'}
      % test presentation
      events(i).study_cond = '';
      events(i).study_resp = '';
      events(i).study_rt = -1;
      events(i).pool_correct = -1;
      events(i).src_isTarg = -1;
      % recognition response
      events(i+1).study_cond = '';
      events(i+1).study_resp = '';
      events(i+1).study_rt = -1;
      events(i+1).pool_correct = -1;
      events(i+1).src_isTarg = -1;
      % source response
      events(i+2).study_cond = '';
      events(i+2).study_resp = '';
      events(i+2).study_rt = -1;
      events(i+2).pool_correct = -1;
      events(i+2).src_isTarg = -1;
  end
end

% put fields in an orderly manner
events = orderfields(events,{'subject','session','mstime','msoffset','list','trial','type','item','pool','pool_correct','xcoord','ycoord','serialpos','study_cond','study_resp','study_rt','rec_isTarg','rec_correct','testType','src_isTarg','src_resp','src_rt','src_correct','rkn_resp','rkn_rt','rkn_correct'});
