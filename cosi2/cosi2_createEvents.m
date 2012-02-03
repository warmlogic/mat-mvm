function events = cosi2_createEvents(dataroot,subject,session,nsFile)
% function events = cosi2_createEvents(dataroot,subject,session,nsFile)
%
% create event struct for COSI2
%
% struct fields:
%   subject
%   session
%   nsFile (only if processing EEG)
%   mstime
%   msoffset
%   list
%   numColors
%   trial
%   type
%   cond
%   item
%   serialpos
%   study_loc_x
%   study_loc_y
%   study_color
%   study_color_x
%   lure_color
%   lure_color_x
%   rec_isTarg
%   rec_correct: didn't say "new" to an old item, but can get source wrong
%   src_isTarg
%   src_resp
%   src_rt
%   src_correct
%   rkn_resp
%   rkn_rt
%   rkn_correct
%   redo
%
% Color:
% For source response, correct is calculated relative to a " Blue" (B)
% response (i.e., "hit" if B response for B source, and "false alarm" if B
% response for Yellow source; CR if Y to Y source, M if Y to B source).
% This could be switched around to be relative to Y.
%
% Side:
% For source response, correct is calculated relative to a "right
% side" (R) response (i.e., "hit" if R response for R source, and
% "false alarm" if R response for L source; CR if L to L source, M if L to
% R). This could be switched around to be relative to L.
%
%
% 'redo' codes: 0=not redone. 1=changed to correct. 2=changed to incorrect.
% 3=same.
%

%dataroot = '/Volumes/curranlab/Data/COSI2/eeg/behavioral';
%subject = 'COSI2001';
%session = 'session_0';

if nargin < 4
  nsFile = '';
end

sesDir = fullfile(dataroot,subject,session);

logfile = fullfile(sesDir,'session.log');

fid = fopen(logfile);
logdata = textscan(fid,'%n%n%s%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);

% constants
NUMLISTS = 4;
MAXTESTTIME = 30000;
REDO_STR = 'REDO';
NEW_STR = 'NEW';

% this only works with two colors
SRC_TARG_COLOR = 'Blue';
SRC_LURE_COLOR = 'Yellow';

% find the places where they did a redo response; will exclude these
redoInd = strcmp(logdata{5},REDO_STR);
if ~isempty(redoInd)
  % find the redone source responses; we want to keep these and mark them
  redoToKeep = double(redoInd);
  % 1. do want to keep the response following REDO (must get set first!!)
  redoToKeep(find(redoInd == 1) + 1) = 1;
  % 2. don't want to keep the actual redo response
  redoToKeep(redoInd == 1) = -1;
  % 3. don't want to keep the source response that came before it; by being
  % set last this also gets rid of any multiple redos for a single trial
  redoToKeep(find(redoInd == 1) - 1) = -1;
  
  % this is where the final redo source responses occurred
  redoToKeepInd = find(redoToKeep == 1);
  
  % find out whether the final redo source responses were accurate
  redoRespAcc = str2double(logdata{8}(redoToKeepInd));
  
  % get indices of actual redo responses to see if there were multiple
  % redos for a single trial
  redoRespInd = find(redoInd == 1);
  if length(redoRespInd) ~= length(redoRespAcc)
    multiRedos = diff(redoRespInd) == 2;
    % add a first 0 value because the first redo response will always be
    % the one we want. then find the other locations where multiple redos
    % were not contiguous (i.e., diff == 2)
    redoRespInd = redoRespInd([0; multiRedos] == 0);
  end
  
  % find out whether the first source responses were accurate
  firstRespAcc = str2double(logdata{8}(redoRespInd - 1));
  
  % Add more info: 1=change to correct, 2=change to incorrect, 3=same.
  for i = 1:length(redoRespAcc)
    if firstRespAcc(i) == 0 && redoRespAcc(i) == 1
      redoToKeep(redoToKeepInd(i)) = 1;
    elseif firstRespAcc(i) == 1 && redoRespAcc(i) == 0
      redoToKeep(redoToKeepInd(i)) = 2;
    elseif firstRespAcc(i) == 1 && redoRespAcc(i) == 1
      redoToKeep(redoToKeepInd(i)) = 3;
    elseif firstRespAcc(i) == 0 && redoRespAcc(i) == 0
      redoToKeep(redoToKeepInd(i)) = 3;
    end
  end
  
  % remove the redo data
  newlogdata = cell(size(logdata));
  for i = 1:length(logdata)
    newlogdata{i} = logdata{i}(redoToKeep ~= -1);
  end
  %oldlogdata = logdata;
  logdata = newlogdata;
  redoToKeep = redoToKeep(redoToKeep ~= -1);
  % put the redo info in the item presentation
  redoToKeepInd = find(redoToKeep ~= 0);
  for i = 1:length(redoToKeepInd)
    redoToKeep(redoToKeepInd(i) - 1) = redoToKeep(redoToKeepInd(i));
  end
else
  % otherwise if they made no redo responses then it gets set to all 0s
  redoToKeep = redoInd;
end

sesNum = str2double(strrep(session,'session_',''));
if isempty(nsFile)
  log = struct('subject',subject,'session',sesNum,'mstime',num2cell(logdata{1}),'msoffset',num2cell(logdata{2}),'trial',[],'type',logdata{3},'redo',num2cell(redoToKeep));
else
  log = struct('subject',subject,'session',sesNum,'nsFile',nsFile,'mstime',num2cell(logdata{1}),'msoffset',num2cell(logdata{2}),'trial',[],'type',logdata{3},'redo',num2cell(redoToKeep));
end

% initialize
listNum = NaN;
numColors = cell(1,NUMLISTS);
numColor = NaN;

% get the test types
for i = 1:length(log)
  if strcmp('COLORS',logdata{3}{i})
    numColors{1} = str2double(logdata{4}{i});
    numColors{2} = str2double(logdata{5}{i});
    numColors{3} = str2double(logdata{6}{i});
    numColors{4} = str2double(logdata{7}{i});
    break
  end
end

for i = 1:length(log)
  log(i).mstime = log(i).mstime;
  log(i).msoffset = log(i).msoffset;
  
  % set the list number and test type
  if strcmp('TRIAL',logdata{3}{i}) % TRIAL = listNum
    listNum = str2double(logdata{4}{i});
    numColor = numColors{listNum};
    % if we're starting a new list, set trialNum to zero
    trialNum = 0;
  end
  log(i).list = listNum;
  if strcmp(logdata{4}{i},'color')
    log(i).numColors = numColor;
  elseif strcmp(logdata{4}{i},'side')
    log(i).numColors = -1;
  end
  
  switch log(i).type
    
    case {'STUDY_TARGET','STUDY_BUFFER'}
      % get the condition
      log(i).cond = logdata{4}{i};
      % get the image name
      log(i).item = logdata{6}{i};
      % mark if it's a target (1) or buffer (0)
      log(i).rec_isTarg = str2double(logdata{7}{i});
      % find out where this item was shown on the screen
      log(i).study_loc_x = str2num(logdata{8}{i});
      log(i).study_loc_y = str2num(logdata{9}{i});
      % get the serial position
      log(i).serialpos = str2double(logdata{10}{i});
      % get the trial number within this list
      log(i).trial = log(i).serialpos;
      if strcmp(log(i).cond,'color')
        % get the study color
        log(i).study_color = logdata{11}{i};
        %log(i).study_color_rgb = logdata{12}{i};
      elseif strcmp(log(i).cond,'side')
        log(i).study_color = 'none';
        %log(i).study_color_rgb = 'none';
      end
      
    case {'TEST_TARGET','TEST_LURE'}
      trialNum = trialNum + 1;
      % get the condition
      log(i).cond = logdata{4}{i};
      % get the image name
      log(i).item = logdata{6}{i};
      % mark if it was a target
      log(i).rec_isTarg = str2double(logdata{7}{i});
      % get the serial position
      log(i).serialpos = str2double(logdata{8}{i});
      if strcmp(log(i).cond,'side')
        % find out where this item was shown on the screen
        log(i).study_loc_x = str2num(logdata{8}{i});
        log(i).study_loc_y = str2num(logdata{9}{i});
        log(i).study_color = 'none';
        log(i).study_color_x = -1;
        log(i).lure_color = 'none';
        log(i).lure_color_x = -1;
      elseif strcmp(log(i).cond,'color')
        log(i).study_loc_x = 0.500;
        log(i).study_loc_y = 0.500;
        log(i).study_color = logdata{9}{i};
        log(i).study_color_x = str2num(logdata{10}{i});
        log(i).lure_color = logdata{11}{i};
        log(i).lure_color_x = str2num(logdata{12}{i});
      end
      % get the trial number within this list
      log(i).trial = trialNum;
      
    case 'SOURCE_RESP'
      % get the condition
      log(i).cond = log(i-1).cond;
      % put in the item name
      log(i).item = log(i-1).item;
      % put in the target status
      log(i).rec_isTarg = log(i-1).rec_isTarg;
      % put in the serial position in which it was shown during study
      log(i).serialpos = log(i-1).serialpos;
      % put in the study locs and study/lure colors
      log(i).study_loc_x = log(i-1).study_loc_x;
      log(i).study_loc_y = log(i-1).study_loc_y;
      log(i).study_color = log(i-1).study_color;
      log(i).study_color_x = log(i-1).study_color_x;
      log(i).lure_color = log(i-1).lure_color;
      log(i).lure_color_x = log(i-1).lure_color_x;
      log(i).numColors = log(i-1).numColors;
      % put in the trial number within this list
      log(i).trial = log(i-1).trial;
      
      % get the response
      log(i).src_resp = logdata{5}{i};
      % was it correct?
      log(i).src_correct = str2double(logdata{8}{i});
      % get the reaction time
      log(i).src_rt = str2double(logdata{6}{i});
      if log(i).src_rt == MAXTESTTIME
        log(i).src_resp = '';
        log(i).src_correct = -1;
      end
      
      % determine if they correctly recognized it as either old
      % (independent of source accuracy) or new
      if log(i).rec_isTarg == 1 && ~strcmp(log(i).src_resp,NEW_STR)
        % hit
        log(i).rec_correct = 1;
      elseif log(i).rec_isTarg == 1 && strcmp(log(i).src_resp,NEW_STR)
        % miss
        log(i).rec_correct = 0;
      elseif log(i).rec_isTarg == 0 && strcmp(log(i).src_resp,NEW_STR)
        % correct rejection
        log(i).rec_correct = 1;
      elseif log(i).rec_isTarg == 0 && ~strcmp(log(i).src_resp,NEW_STR)
        % false alarm
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
      
      % get the condition
      log(i).cond = log(i-2).cond;
      % put in the item name
      log(i).item = log(i-2).item;
      % put in the target status
      log(i).rec_isTarg = log(i-2).rec_isTarg;
      % put in the serial position in which it was shown during study
      log(i).serialpos = log(i-2).serialpos;
      % put in the study locs and study/lure colors
      log(i).study_loc_x = log(i-2).study_loc_x;
      log(i).study_loc_y = log(i-2).study_loc_y;
      log(i).study_color = log(i-2).study_color;
      log(i).study_color_x = log(i-2).study_color_x;
      log(i).lure_color = log(i-2).lure_color;
      log(i).lure_color_x = log(i-2).lure_color_x;
      log(i).numColors = log(i-1).numColors;
      % put in the trial number within this list
      log(i).trial = log(i-2).trial;
      
      % get the reaction time
      log(i).rkn_rt = str2double(logdata{6}{i});
      % was it correct?
      log(i).rkn_correct = str2double(logdata{8}{i});
      % get the response
      log(i).rkn_resp = logdata{5}{i};
      if log(i).rkn_rt == MAXTESTTIME
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
events = filterStruct(log,'ismember(type,varargin{1})',{'STUDY_BUFFER','STUDY_TARGET','TEST_TARGET','TEST_LURE','SOURCE_RESP','RK_RESP','NEW_RESP'});

% put the info from test into the study events, and the info from
% study into the test events
for i = 1:length(events)
  switch events(i).type
    
    case {'STUDY_BUFFER'}
      % if it was a buffer, put in -1s
      
      % study presentation - study info
      events(i).study_color = 'none';
      events(i).study_color_x = -1;
      %events(i).study_color_rgb = -1;
      events(i).lure_color = 'none';
      events(i).lure_color_x = -1;
      % study presentation - test info
      events(i).rec_correct = -1;
      events(i).src_isTarg = -1;
      events(i).src_resp = '';
      events(i).src_rt = -1;
      events(i).src_correct = -1;
      events(i).rkn_resp = '';
      events(i).rkn_rt = -1;
      events(i).rkn_correct = -1;
      
    case {'STUDY_TARGET'}
      % if we find a study event, get the corresponding test event
      % and move the test info into the study event
      
      testEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'TEST_TARGET'},{events(i).item});
      if length(testEvent) > 1
        error('Found multiple events when searching for study target %s!',events(i).item);
      elseif isempty(testEvent)
        % study presentation - study info
        events(i).study_color = 'none';
        events(i).study_color_x = -1;
        events(i).lure_color = 'none';
        events(i).lure_color_x = -1;
        % study presentation - test info
        events(i).rec_correct = -1;
        events(i).src_isTarg = -1;
        events(i).src_resp = '';
        events(i).src_rt = -1;
        events(i).src_correct = -1;
        events(i).rkn_resp = '';
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
        
        if isempty(events(i).study_color)
          events(i).study_color = testEvent.study_color;
        end
        events(i).study_color_x = testEvent.study_color_x;
        events(i).lure_color = testEvent.lure_color;
        events(i).lure_color_x = testEvent.lure_color_x;
        
        if strcmp(events(i).cond,'color')
          % put in the target status where all RIGHT SIDE TARGET COLOR
          % sources are targets and all LEFT SIDE TARGET COLOR sources are
          % lures; -1 if otherwise; this is described at the top of the
          % function
          if events(i).rec_isTarg == 1 && strcmp(events(i).study_color,SRC_TARG_COLOR)
            % if source color was SRC_TARG_COLOR then it's a target
            events(i).src_isTarg = 1;
          elseif events(i).rec_isTarg == 1 && strcmp(events(i).study_color,SRC_LURE_COLOR)
            % if source color was SRC_LURE_COLOR then it's a lure
            events(i).src_isTarg = 0;
          end
        elseif strcmp(events(i).cond,'side')
          % put in the target status where all RIGHT SIDE sources are targets
          % and all LEFT SIDE  sources are lures; -1 if otherwise; this is
          % described at the top of the function
          if events(i).rec_isTarg == 1 && (events(i).study_loc_x > 0.5 && events(i).study_loc_x <= 1)
            % if it was on right then it's a target
            events(i).src_isTarg = 1;
          elseif events(i).rec_isTarg == 1 && (events(i).study_loc_x < 0.5 && events(i).study_loc_x >= 0)
            % if it was on left then it's a lure
            events(i).src_isTarg = 0;
          end
        end
      end
      
    case {'TEST_TARGET'}
      % if we find a test event, get the corresponding study event
      % and move the study info into the test event
      
      studyEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'STUDY_TARGET'},{events(i).item});
      if length(studyEvent) > 1
        error('Found multiple events when searching for study target %s!',events(i).item);
      elseif isempty(studyEvent)
        keyboard
      elseif length(studyEvent) == 1
        events(i).src_isTarg = studyEvent.src_isTarg;
        events(i+1).src_isTarg = studyEvent.src_isTarg;
        events(i+2).src_isTarg = studyEvent.src_isTarg;
        
        events(i).study_loc_x = studyEvent.study_loc_x;
        events(i+1).study_loc_x = studyEvent.study_loc_x;
        events(i+2).study_loc_x = studyEvent.study_loc_x;
        events(i).study_loc_y = studyEvent.study_loc_y;
        events(i+1).study_loc_y = studyEvent.study_loc_y;
        events(i+2).study_loc_y = studyEvent.study_loc_y;
      end
      
    case {'TEST_LURE'}
      % if we find a lure test event, put in -1s
      
      % test presentation
      events(i).src_isTarg = -1;
      % recognition response
      events(i+1).src_isTarg = -1;
      % source response
      events(i+2).src_isTarg = -1;
      
      events(i).study_loc_x = -1;
      events(i+1).study_loc_x = -1;
      events(i+2).study_loc_x = -1;
      events(i).study_loc_y = -1;
      events(i+1).study_loc_y = -1;
      events(i+2).study_loc_y = -1;
  end
end

% put fields in an orderly manner
if isempty(nsFile)
  events = orderfields(events,{'subject','session','mstime','msoffset','list','numColors','trial','type','cond','item','serialpos','study_loc_x','study_loc_y','study_color','study_color_x','lure_color','lure_color_x','rec_isTarg','rec_correct','src_isTarg','src_resp','src_rt','src_correct','rkn_resp','rkn_rt','rkn_correct','redo'});
else
  events = orderfields(events,{'subject','session','nsFile','mstime','msoffset','list','numColors','trial','type','cond','item','serialpos','study_loc_x','study_loc_y','study_color','study_color_x','lure_color','lure_color_x','rec_isTarg','rec_correct','src_isTarg','src_resp','src_rt','src_correct','rkn_resp','rkn_rt','rkn_correct','redo'});
end
