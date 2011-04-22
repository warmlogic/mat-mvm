function events = rsrc2_createEvents(dataroot,subject,session)
% function events = rsrc2_events(dataroot,subject,session)
%
% create event struct for RSRC
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
%   xcoord
%   ycoord
%   serialpos
%   rec_isTarg
%   rec_correct
%   rec_resp
%   rec_conf
%   rec_rt
%   src_isTarg
%   src_correct
%   src_resp
%   src_conf
%   src_rt
%
% For source response, correct is calculated relative to a "right
% side" (R) response (i.e., "hit" if R response for R source, and
% "false alarm" if R response for L source). This could be switched
% around to be relative to L.
%

%dataroot = '/Users/matt/data/RSRC/behavioral';
%subject = 'RSRC001';
%session = 'session_0';

% confidence numbers are [1 2 3 4 5 6], sure old to sure new, same for
% sure L to sure R
confRange = [1 2 3 4 5 6];
confRange_rec_str = {'OLD_SURE','OLD_PROBABLY','OLD_GUESS','NEW_GUESS','NEW_PROBABLY','NEW_SURE'};
confRange_src_str = {'LEFT_SURE','LEFT_PROBABLY','LEFT_GUESS','RIGHT_GUESS','RIGHT_PROBABLY','RIGHT_SURE'};

% confKeys = ['Z','X','C',',','.','/']; % keys, sure old/new to sure
%                                       % new/old, and sure L to sure R
% if str2num(subject(end-2:end)) < 21
%   confRange_rec = [6 5 4 3 2 1]; % sure old to sure new
% else
%   confRange_rec = [1 2 3 4 5 6]; % sure new to sure old
% end
% confRange_src = [1 2 3 4 5 6]; % sure L to sure R

sesDir = fullfile(dataroot,subject,session);

logfile = fullfile(sesDir,'session.log');

fid = fopen(logfile);
logdata = textscan(fid,'%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);

sesNum = str2num(strrep(session,'session_',''));
log = struct('subject',subject,'session',sesNum,'mstime',logdata{1},'msoffset',logdata{2},'trial',[],'type',logdata{3});

% constants
maxTestTime = 30000;

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
    log(i).rec_resp = -1;
    log(i).rec_conf = -1;
    log(i).rec_rt = -1;
    log(i).rec_correct = -1;
    % mark that it's a buffer (and not a target)
    log(i).src_isTarg = str2num(logdata{6}{i});
    % put in -1 where we won't have info
    log(i).src_resp = -1;
    log(i).src_conf = -1;
    log(i).src_rt = -1;
    log(i).src_correct = -1;
    
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
    % we'll mark the source target status later
    
   case {'TARG_PRES','LURE_PRES'}
    trialNum = trialNum + 1;
    % get the image name
    log(i).item = logdata{5}{i};
    % mark that it was a target
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
    % put in the serial position in which it was shown
    log(i).serialpos = log(i-1).serialpos;
    % and the test list trial number
    log(i).trial = log(i-1).trial;
    % put in the target status
    log(i).rec_isTarg = log(i-1).rec_isTarg;
    
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
    if log(i).rec_rt == maxTestTime
      log(i).rec_resp = '';
      log(i).rec_correct = -1;
    end
    
    % put info in test presentation
    log(i-1).rec_rt = log(i).rec_rt;
    
   case 'SOURCE_RESP'
    % bring the RECOG_RESP info forward
    log(i).rec_isTarg = log(i-1).rec_isTarg;
    log(i).rec_resp = log(i-1).rec_resp;
    log(i).rec_conf = log(i-1).rec_conf;
    log(i).rec_rt = log(i-1).rec_rt;
    log(i).rec_correct = log(i-1).rec_correct;
    
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
    % put it in the test presentation
    log(i-2).src_resp = log(i).src_resp;
    log(i-1).src_resp = log(i).src_resp;
    
    % get the confidence level response and convert it to a number; sure
    % old to sure new
    if strcmp(logdata{5}{i},confRange_src_str{1})
      log(i).src_conf = confRange(1);
    elseif strcmp(logdata{5}{i},confRange_src_str{2})
      log(i).src_conf = confRange(2);
    elseif strcmp(logdata{5}{i},confRange_src_str{3})
      log(i).src_conf = confRange(3);
    elseif strcmp(logdata{5}{i},confRange_src_str{4})
      log(i).src_conf = confRange(4);
    elseif strcmp(logdata{5}{i},confRange_src_str{5})
      log(i).src_conf = confRange(5);
    elseif strcmp(logdata{5}{i},confRange_src_str{6})
      log(i).src_conf = confRange(6);
    end
    % put info in test presentation
    log(i-2).src_conf = log(i).src_conf;
    log(i-1).src_conf = log(i).src_conf;
    
    % get the reaction time
    log(i).src_rt = str2num(logdata{6}{i});
    % put info in recognition response
    log(i-1).src_rt = log(i).src_rt;
    % put info in test presentation
    log(i-2).src_rt = log(i).src_rt;
    
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
    
    % mark whether they got the source correct; this relative to R method
    % is described at the top of the function.
    %
    % if it was an R source (isTarg = 1) marked as 'R', it is correct;
    % if it was an R source (isTarg = 1) marked as 'L', it is wrong;
    % otherwise mark it as -1;
    if log(i).src_isTarg == 1 % right source
      if ismember(log(i).src_conf,confRange(4:6)) % right conf
        log(i).src_correct = 1; % hit
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct; % hit
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct; % hit
      elseif ismember(log(i).src_conf,confRange(1:3)) % left conf
        log(i).src_correct = 0; % miss
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct; % miss
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct; % miss
      end
    elseif log(i).src_isTarg == 0 % left source
      if ismember(log(i).src_conf,confRange(4:6)) % right conf
        log(i).src_correct = 0; % false alarm
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct; % false alarm
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct; % false alarm
      elseif ismember(log(i).src_conf,confRange(1:3)) % left conf
        log(i).src_correct = 1; % correct rejection
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct; % correct rejection
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct; % correct rejection
      end
    elseif log(i).src_isTarg == -1 % real lure
      log(i).src_correct = -1; % real lure
      % put info in recognition response
      log(i-1).src_correct = log(i).src_correct; % real lure
      % put info in test presentation
      log(i-2).src_correct = log(i).src_correct; % real lure
    end
    
  end % switch event type
end % for i = 1:length(log)

% grab only the events we want from our log struct
events = filterStruct(log,'ismember(type,varargin{1})',{'STUDY_BUFFER','STUDY_TARGET','TARG_PRES','LURE_PRES','RECOG_RESP','SOURCE_RESP'});
% put fields in an orderly manner
events = orderfields(events,{'subject','session','mstime','msoffset','list','trial','type','item','xcoord','ycoord','serialpos','rec_isTarg','rec_correct','rec_resp','rec_conf','rec_rt','src_isTarg','src_correct','src_resp','src_conf','src_rt'});

% put the info from test into the study events
for i = 1:length(events)
  if strcmp(events(i).type,'STUDY_TARGET')
    testEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'TARG_PRES'},{events(i).item});
    if length(testEvent) > 1
      error('Found multiple events when searching for study target %s!',events(i).item);
    else
      events(i).rec_isTarg = testEvent.rec_isTarg;
      events(i).rec_resp = testEvent.rec_resp;
      events(i).rec_conf = testEvent.rec_conf;
      events(i).rec_rt = testEvent.rec_rt;
      events(i).rec_correct = testEvent.rec_correct;
      events(i).src_isTarg = testEvent.src_isTarg;
      events(i).src_resp = testEvent.src_resp;
      events(i).src_conf = testEvent.src_conf;
      events(i).src_rt = testEvent.src_rt;
      events(i).src_correct = testEvent.src_correct;
    end
  end
end
