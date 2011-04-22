function events = grub_createEvents(dataroot,subject,session)
% function events = grub_events(dataroot,subject,session)
%
% create event struct for GRUB
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
%   rec_resp
%   rec_rt
%   rec_correct
%   src_isTarg
%   src_resp
%   src_rt
%   src_correct
%
% For source responses, accuracy is calculated such that the size source
% is the target distribution and the life source is the lure distribution.
% Here, the designation of the target distribution is arbitrary; the same
% results would be obtained if the life source was the target distribution
% and the size source was the lure distribution.
%
% Hit if "SIZE" response to SIZE source
% Miss if "LIFE" response to SIZE source
% CR if "LIFE" response to LIFE source
% FA if "SIZE" response to LIFE source
% 

%dataroot = '/Volumes/curranlab/Data/GRUB/eeg/behavioral';
%subject = 'GRUB002';
%session = 'session_0';

%confKeys = ['Z','/']; % keys, sure old to sure new, or
                                      % sure L to sure R
%confRange_src = [1 2 3 4 5 6]; % sure L to sure R

sesDir = fullfile(dataroot,subject,session);

logfile = fullfile(sesDir,'session.log');

fid = fopen(logfile);
logdata = textscan(fid,'%s%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);

%[l1 l2 l3 l4 l5 l6 l7 l8] = textread(logfile,'%d%d%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);

% constants
maxStudyTime = 20000;
maxTestTime = 20000;
expCompImageDir = '/home/curranlab/experiments/pyGrubEtal08/images/';

sesNum = str2num(strrep(session,'session_',''));
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
    
   case {'STUDY_TARGET','STUDY_BUFFER'}
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
    % get the study condition
    log(i).study_cond = logdata{10}{i};
    
   case 'STUDY_RESP'
    log(i).study_resp = logdata{5}{i};
    log(i).study_rt = str2num(logdata{6}{i});
    % if they hit the max time, they didn't answer; bug in exp (3Mar2010)
    if log(i).study_rt == maxStudyTime
      log(i).study_resp = '';
    end
    % get info from presentation
    log(i).item = log(i-1).item;
    log(i).rec_isTarg = log(i-1).rec_isTarg;
    log(i).xcoord = log(i-1).xcoord;
    log(i).ycoord = log(i-1).ycoord;
    log(i).serialpos = log(i-1).serialpos;
    log(i).trial = log(i-1).trial;
    log(i).study_cond = log(i-1).study_cond;
    % put info in the presentation event
    %log(i-1).study_cond = log(i).study_cond;
    log(i-1).study_resp = log(i).study_resp;
    log(i-1).study_rt = log(i).study_rt;
    
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
    
   case 'RECOG_RESP'
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
    log(i).rec_resp = logdata{5}{i};
    % get the reaction time
    log(i).rec_rt = str2num(logdata{6}{i});
    if log(i).rec_rt == maxTestTime
      log(i).rec_resp = '';
    end
    % was it correct?
    log(i).rec_correct = str2num(logdata{8}{i});
    
    % put info in test presentation
    log(i-1).rec_resp = log(i).rec_resp;
    log(i-1).rec_rt = log(i).rec_rt;
    log(i-1).rec_correct = log(i).rec_correct;
    
%     if log(i).rec_isTarg == 1 & log(i).rec_correct == 0
%         log(i).src_correct = 0;
%         log(i-1).src_correct = 0;
%     end
%     if log(i).rec_isTarg == 0 & log(i).rec_correct == 1
%         log(i).src_correct = 1;
%         log(i-1).src_correct = 1;
%     end
    
   case 'SOURCE_RESP'
    % bring the RECOG_RESP info forward
    log(i).rec_isTarg = log(i-1).rec_isTarg;
    log(i).rec_resp = log(i-1).rec_resp;
    log(i).rec_rt = log(i-1).rec_rt;
    log(i).rec_correct = log(i-1).rec_correct;
    
    % put in the item name
    log(i).item = log(i-2).item;
    % put in the coordinates
    log(i).xcoord = log(i-2).xcoord;
    log(i).ycoord = log(i-2).ycoord;
    % put in the serial position in which it was shown
    log(i).serialpos = log(i-2).serialpos;
    % put in the trial number within this list
    log(i).trial = log(i-2).trial;
    
    % get the response
    log(i).src_resp = logdata{5}{i};
    % get the reaction time
    log(i).src_rt = str2num(logdata{6}{i});
    % if they hit the max time, they didn't answer; bug in exp (3Mar2010)
    if log(i).src_rt == maxTestTime
      log(i).src_resp = '';
      log(i).src_isTarg = -1;
    end
    % was it correct?
    log(i).src_correct = str2num(logdata{8}{i});
    
%     % put in the target status where all SIZE sources are targets and all
%     % LIFE sources are lures; -1 if otherwise; this is described at the top
%     % of the function
%     if log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,'SIZE') && log(i).src_correct == 1
%       % if it's Size source & they got it right then it's a target (= size)
%       log(i).src_isTarg = 1;
%     elseif log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,'SIZE') && log(i).src_correct == 0
%       % if it's Size source & they got it wrong then it's a lure (= life)
%       log(i).src_isTarg = 0;
%     elseif log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,'ANIMACY') && log(i).src_correct == 1
%       % if it's Life source & they got it right then it's a lure (= life)
%       log(i).src_isTarg = 0;
%     elseif log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,'ANIMACY') && log(i).src_correct == 0
%       % if it's Life source & they got it wrong then it's a target (= size)
%       log(i).src_isTarg = 1;
%     elseif log(i).rec_isTarg == 0
%       % otherwise it was a real lure, and we mark that as -1
%       log(i).src_isTarg = -1;
%     end
    
    % put info in recognition response
    log(i-1).src_resp = log(i).src_resp;
    log(i-1).src_rt = log(i).src_rt;
    log(i-1).src_correct = log(i).src_correct;
    % put info in test presentation
    log(i-2).src_resp = log(i).src_resp;
    log(i-2).src_rt = log(i).src_rt;
    log(i-2).src_correct = log(i).src_correct;
    
  end % switch event type
end % for i = 1:length(log)

% search the video log for the names of the images, the times they
% were presented, and grab the pool name from the path
vidlogfile = fullfile(sesDir,'video.vidlog');
fid = fopen(vidlogfile);
vidlogdata = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);
image_ms = [];
image_pool = {};
image_name = {};
count = 0;
for i = 1:length(vidlogdata{7})
  if strcmp(vidlogdata{7}(i),'IMAGE')
    count = count + 1;
    image_ms = [image_ms;str2num(cell2mat(vidlogdata{1}(i)))];
    pool_imagename = cell2mat(strrep(vidlogdata{8}(i),expCompImageDir,''));
    image_pool{count} = pool_imagename(1:strfind(pool_imagename,'/')-1);
    image_name{count} = pool_imagename(strfind(pool_imagename,'/')+1:end-4);
  end
end

% grab only the events we want from our log struct
events = filterStruct(log,'ismember(type,varargin{1})',{'STUDY_TARGET','STUDY_RESP','TARG_PRES','LURE_PRES','RECOG_RESP','SOURCE_RESP'});

% get rid of practice images
expStartInd = find(image_ms == events(1).mstime);
image_ms = image_ms(expStartInd:end);
image_pool = image_pool(expStartInd:end);
image_name = image_name(expStartInd:end);

for i = 1:length(image_ms)
  for j = 1:length(events)
    % specify which Snake (large or small) and insert the name of the
    % image pool into the field 'pool'; special case for the snake
    % image because we had one large and one small, both with the
    % name 'Snake'; changed as of 9 July 2009
    snakePool = [];
    if strcmp(events(j).item,'Snake')
      snakePool = cell2mat(image_pool(find(events(j).mstime == image_ms)));
      if ~isempty(snakePool)
        % put in the image pool
        events(j).pool = snakePool;
        % rename the item with large or small
        events(j).item = strcat([events(j).item,' ',strrep(snakePool,'living_','')]);
      else
        % if it is a response, we won't find the image_ms,
        % so grab it from the previous event
        events(j).pool = events(j-1).pool;
        events(j).item = events(j-1).item;
      end

    % insert the name of the image pool into the field 'pool'
    elseif strcmp(events(j).item,image_name(i))
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
        events(i).rec_isTarg = testEvent.rec_isTarg;
        events(i).rec_resp = testEvent.rec_resp;
        events(i).rec_rt = testEvent.rec_rt;
        events(i).rec_correct = testEvent.rec_correct;
        %events(i).src_isTarg = testEvent.src_isTarg;
        events(i).src_resp = testEvent.src_resp;
        events(i).src_rt = testEvent.src_rt;
        events(i).src_correct = testEvent.src_correct;
        % and the study response
        events(i+1).rec_isTarg = testEvent.rec_isTarg;
        events(i+1).rec_resp = testEvent.rec_resp;
        events(i+1).rec_rt = testEvent.rec_rt;
        events(i+1).rec_correct = testEvent.rec_correct;
        %events(i+1).src_isTarg = testEvent.src_isTarg;
        events(i+1).src_resp = testEvent.src_resp;
        events(i+1).src_rt = testEvent.src_rt;
        events(i+1).src_correct = testEvent.src_correct;
      end
      
      % study response for if the study classification was correct
      events(i+1).pool_correct = events(i).pool_correct;
      
      % put in the target status where all SIZE sources are targets and all
      % LIFE sources are lures; -1 if otherwise; this is described at the top
      % of the function
      if events(i).rec_isTarg == 1 && strcmp(events(i).study_cond,'SIZE')
        % if it's Size source then it's a target
        events(i).src_isTarg = 1;
        events(i+1).src_isTarg = 1;
      elseif events(i).rec_isTarg == 1 && strcmp(events(i).study_cond,'LIFE')
        % if it's Animacy source then it's a lure
        events(i).src_isTarg = 0;
        events(i+1).src_isTarg = 0;
      elseif events(i).rec_isTarg == 0
        % otherwise it was a real lure, and we mark that as -1
        events(i).src_isTarg = -1;
        events(i+1).src_isTarg = -1;
      end
      
      % if there was no src_resp, put in -1s
      if strcmp(events(i).rec_resp,'NEW')
        % study presentation
        events(i).src_isTarg = -1;
        events(i).src_resp = '';
        events(i).src_rt = -1;
        events(i).src_correct = -1;
        % study response
        events(i+1).src_isTarg = -1;
        events(i+1).src_resp = '';
        events(i+1).src_rt = -1;
        events(i+1).src_correct = -1;
      end
      
      % if we find a test event, get the corresponding study event
      % and move the study info into the test event
    case {'TARG_PRES'}
      studyEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_TARGET'},{events(i).item});
      if length(studyEvent) > 1
        error('Found multiple events when searching for study target %s!',events(i).item);
      else
        % put it in the study presentation
        events(i).study_cond = studyEvent.study_cond;
        events(i).study_resp = studyEvent.study_resp;
        events(i).study_rt = studyEvent.study_rt;
        events(i).pool_correct = studyEvent.pool_correct;
        events(i).src_isTarg = studyEvent.src_isTarg;
        % and the recognition response
        events(i+1).study_cond = studyEvent.study_cond;
        events(i+1).study_resp = studyEvent.study_resp;
        events(i+1).study_rt = studyEvent.study_rt;
        events(i+1).pool_correct = studyEvent.pool_correct;
        events(i+1).src_isTarg = studyEvent.src_isTarg;
        % and the source response if they said "OLD"
        if strcmp(events(i).rec_resp,'OLD')
          events(i+2).study_cond = studyEvent.study_cond;
          events(i+2).study_resp = studyEvent.study_resp;
          events(i+2).study_rt = studyEvent.study_rt;
          events(i+2).pool_correct = studyEvent.pool_correct;
          events(i+2).src_isTarg = studyEvent.src_isTarg;
        elseif strcmp(events(i).rec_resp,'NEW')
          % if there was no src_resp, put in -1s
          %
          % test presentation
          events(i).src_isTarg = -1;
          events(i).src_resp = '';
          events(i).src_rt = -1;
          events(i).src_correct = -1;
          % recognition response
          events(i+1).src_isTarg = -1;
          events(i+1).src_resp = '';
          events(i+1).src_rt = -1;
          events(i+1).src_correct = -1;
        end
      end
      
      % if we find a lure test event, put in -1s
    case {'LURE_PRES'}
      if strcmp(events(i).rec_resp,'OLD')
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
      elseif strcmp(events(i).rec_resp,'NEW')
        % test presentation
        events(i).study_cond = '';
        events(i).study_resp = '';
        events(i).study_rt = -1;
        events(i).pool_correct = -1;
        events(i).src_isTarg = -1;
        events(i).src_resp = '';
        events(i).src_rt = -1;
        events(i).src_correct = -1;
        % recognition response
        events(i+1).study_cond = '';
        events(i+1).study_resp = '';
        events(i+1).study_rt = -1;
        events(i+1).pool_correct = -1;
        events(i+1).src_isTarg = -1;
        events(i+1).src_resp = '';
        events(i+1).src_rt = -1;
        events(i+1).src_correct = -1;
      end
  end
end

% put fields in an orderly manner
events = orderfields(events,{'subject','session','mstime','msoffset','list','trial','type','item','pool','pool_correct','xcoord','ycoord','serialpos','study_cond','study_resp','study_rt','rec_isTarg','rec_resp','rec_rt','rec_correct','src_isTarg','src_resp','src_rt','src_correct'});
