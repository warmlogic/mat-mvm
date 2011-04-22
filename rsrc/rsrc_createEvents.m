function events = rsrc_createEvents(dataroot,subject,session)
% function events = rsrc_events(dataroot,subject,session)
%
% create event struct for RSRC
%
% struct fields:
%   subject
%   session
%   mstime
%   msoffset
%   trial
%   type
%   item
%   xcoord
%   ycoord
%   serialpos
%   rec_isTarg
%   rec_conf
%   rec_rt
%   rec_correct
%   src_isTarg
%   src_conf
%   src_rt
%   src_correct
%
% For source response, correct is calculated relative to a "right
% side" (R) response (i.e., "hit" if R response for R source, and
% "false alarm" if R response for L source). This could be switched
% around to be relative to L.
%

%dataroot = '/Users/matt/data/RSRC/behavioral';
%subject = 'RSRC001';
%session = 'session_0';

confKeys = ['Z','X','C',',','.','/']; % keys, sure old to sure new, or
                                      % sure L to sure R
if str2num(subject(end-2:end)) < 21
  confRange_rec = [6 5 4 3 2 1]; % sure old to sure new
else
  confRange_rec = [1 2 3 4 5 6]; % sure new to sure old
end
confRange_src = [1 2 3 4 5 6]; % sure L to sure R

sesDir = fullfile(dataroot,subject,session);

logfile = fullfile(sesDir,'session.log');

fid = fopen(logfile);
logdata = textscan(fid,'%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
fclose(fid);

%[l1 l2 l3 l4 l5 l6 l7 l8] = textread(logfile,'%d%d%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);

sesNum = str2num(strrep(session,'session_',''));
log = struct('subject',subject,'session',sesNum,'mstime',logdata{1},'msoffset',logdata{2},'trial',[],'type',logdata{3});

trialNum = NaN;

for i = 1:length(log)
  % convert to numbers
  log(i).mstime = str2num(log(i).mstime);
  log(i).msoffset = str2num(log(i).msoffset);
  
  % set the trial number
  if strcmp('TRIAL',logdata{3}{i})
    trialNum = str2num(logdata{4}{i});
  end
  log(i).trial = trialNum;
  
  switch log(i).type
    
   case 'STUDY_BUFFER'
    % get the image name
    log(i).item = logdata{5}{i};
    % find out where this item was shown on the screen
    log(i).xcoord = str2num(logdata{7}{i});
    log(i).ycoord = str2num(logdata{8}{i});
    % get the serial position
    log(i).serialpos = str2num(logdata{9}{i});
    % mark that it's a buffer (and not a target)
    log(i).rec_isTarg = str2num(logdata{6}{i});
    % put in -1 where we won't have info
    log(i).rec_conf = -1;
    log(i).rec_rt = -1;
    log(i).rec_correct = -1;
    % mark that it's a buffer (and not a target)
    log(i).src_isTarg = str2num(logdata{6}{i});
    % put in -1 where we won't have info
    log(i).src_conf = -1;
    log(i).src_rt = -1;
    log(i).src_correct = -1;
    
   case 'STUDY_TARGET'
    % get the image name
    log(i).item = logdata{5}{i};
    % find out where this item was shown on the screen
    log(i).xcoord = str2num(logdata{7}{i});
    log(i).ycoord = str2num(logdata{8}{i});
    % get the serial position
    log(i).serialpos = str2num(logdata{9}{i});
    % mark that it's a target (and not a buffer)
    log(i).rec_isTarg = str2num(logdata{6}{i});
    % we'll mark the source target status later
    
   case 'TARG_PRES'
    % get the image name
    log(i).item = logdata{5}{i};
    % find out where this item was shown on the screen
    log(i).xcoord = str2num(logdata{7}{i});
    log(i).ycoord = str2num(logdata{8}{i});
    % get the serial position
    log(i).serialpos = str2num(logdata{9}{i});
    % mark that it was a target
    log(i).rec_isTarg = str2num(logdata{6}{i});
    
   case 'LURE_PRES'
    % get the image name
    log(i).item = logdata{5}{i};
    % find out where this item was shown on the screen
    log(i).xcoord = str2num(logdata{7}{i});
    log(i).ycoord = str2num(logdata{8}{i});
    % get the serial position
    log(i).serialpos = str2num(logdata{9}{i});
    % mark that it was not a target
    log(i).rec_isTarg = str2num(logdata{6}{i});
    
   case 'RECOG_RESP'
    % put in the item name
    log(i).item = log(i-1).item;
    
    % put in the coordinates
    log(i).xcoord = log(i-1).xcoord;
    log(i).ycoord = log(i-1).ycoord;
    
    % put in the serial position in which it was shown
    log(i).serialpos = log(i-1).serialpos;
    
    % put in the target status
    log(i).rec_isTarg = log(i-1).rec_isTarg;
    
    % get the confidence level response and convert it to a number; sure
    % old to sure new
    rec_conf = NaN;
    if strcmp(logdata{4}{i},confKeys(1))
      rec_conf = confRange_rec(1);
    elseif strcmp(logdata{4}{i},confKeys(2))
      rec_conf = confRange_rec(2);
    elseif strcmp(logdata{4}{i},confKeys(3))
      rec_conf = confRange_rec(3);
    elseif strcmp(logdata{4}{i},confKeys(4))
      rec_conf = confRange_rec(4);
    elseif strcmp(logdata{4}{i},confKeys(5))
      rec_conf = confRange_rec(5);
    elseif strcmp(logdata{4}{i},confKeys(6))
      rec_conf = confRange_rec(6);
    end
    log(i).rec_conf = rec_conf;
    % put info in test presentation
    log(i-1).rec_conf = rec_conf;
    
    % get the reaction time
    log(i).rec_rt = str2num(logdata{5}{i});
    % put info in test presentation
    log(i-1).rec_rt = str2num(logdata{5}{i});
    
    if str2num(subject(end-2:end)) < 21
      % if it was an old item and they marked it as 'old', mark it correct;
      % if it was an old item and they marked it as 'new', mark it wrong
      if log(i).rec_isTarg == 1
        if ismember(log(i).rec_conf,confRange_rec(1:3)) % old resp
          log(i).rec_correct = 1; % hit
          % put info in test presentation
          log(i-1).rec_correct = 1; % hit
        elseif ismember(log(i).rec_conf,confRange_rec(4:6)) % new resp
          log(i).rec_correct = 0; % miss
          % put info in test presentation
          log(i-1).rec_correct = 0; % miss
        end
      elseif log(i).rec_isTarg == 0
        if ismember(log(i).rec_conf,confRange_rec(1:3)) % old resp
          log(i).rec_correct = 0; % false alarm
          % put info in test presentation
          log(i-1).rec_correct = 0; % false alarm
        elseif ismember(log(i).rec_conf,confRange_rec(4:6)) % new resp
          log(i).rec_correct = 1; % correct rejection
          % put info in test presentation
          log(i-1).rec_correct = 1; % correct rejection
        end
      end
    else % for the newer eeg subjects
      % if it was an old item and they marked it as 'old', mark it correct;
      % if it was an old item and they marked it as 'new', mark it wrong
      if log(i).rec_isTarg == 1
        if ismember(log(i).rec_conf,confRange_rec(4:6)) % old resp
          log(i).rec_correct = 1; % hit
          % put info in test presentation
          log(i-1).rec_correct = 1; % hit
        elseif ismember(log(i).rec_conf,confRange_rec(1:3)) % new resp
          log(i).rec_correct = 0; % miss
          % put info in test presentation
          log(i-1).rec_correct = 0; % miss
        end
      elseif log(i).rec_isTarg == 0
        if ismember(log(i).rec_conf,confRange_rec(4:6)) % old resp
          log(i).rec_correct = 0; % false alarm
          % put info in test presentation
          log(i-1).rec_correct = 0; % false alarm
        elseif ismember(log(i).rec_conf,confRange_rec(1:3)) % new resp
          log(i).rec_correct = 1; % correct rejection
          % put info in test presentation
          log(i-1).rec_correct = 1; % correct rejection
        end
      end
    end
    
   case 'SOURCE_RESP'
    % bring the RECOG_RESP info forward
    log(i).rec_isTarg = log(i-1).rec_isTarg;
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
    
    % get the confidence level response and convert it to a number; sure
    % left to sure right
    src_conf = NaN;
    if strcmp(logdata{4}{i},confKeys(1))
      src_conf = confRange_src(1);
    elseif strcmp(logdata{4}{i},confKeys(2))
      src_conf = confRange_src(2);
    elseif strcmp(logdata{4}{i},confKeys(3))
      src_conf = confRange_src(3);
    elseif strcmp(logdata{4}{i},confKeys(4))
      src_conf = confRange_src(4);
    elseif strcmp(logdata{4}{i},confKeys(5))
      src_conf = confRange_src(5);
    elseif strcmp(logdata{4}{i},confKeys(6))
      src_conf = confRange_src(6);
    end
    log(i).src_conf = src_conf;
    % put info in recognition response
    log(i-1).src_conf = log(i).src_conf;
    % put info in test presentation
    log(i-2).src_conf = log(i).src_conf;
    
    % get the reaction time
    log(i).src_rt = str2num(logdata{5}{i});
    % put info in recognition response
    log(i-1).src_rt = log(i).src_rt;
    % put info in test presentation
    log(i-2).src_rt = log(i).src_rt;
    
%     % old method:
%     % mark whether they got the source correct
%     if log(i).isTarg % targ
%       if log(i).xcoord > 0.5 & log(i).xcoord <= 1 % right source
%         if ismember(log(i).conf,confRange(1:3)) % left conf
%           log(i).src_correct = 0;
%         elseif ismember(log(i).conf,confRange(4:6)) % right conf
%           log(i).src_correct = 1;
%         end
%       elseif log(i).xcoord < 0.5 & log(i).xcoord >= 0 % left source
%         if ismember(log(i).conf,confRange(1:3)) % left conf
%           log(i).src_correct = 1;
%         elseif ismember(log(i).conf,confRange(4:6)) % right conf
%           log(i).src_correct = 0;
%         end
%       end
%     elseif log(i).isTarg == 0 % lure
%       log(i).src_correct = -1;
%     end
    
%     % I think this is incorrect:
%     % mark whether they got the source correct - this is relative to R
%     % response as described at top of function
%     if log(i).isTarg % targ
%       if log(i).xcoord > 0.5 & log(i).xcoord <= 1 % right source
%         if ismember(log(i).conf,confRange(1:3)) % left conf
%           log(i).src_correct = -1;
%         elseif ismember(log(i).conf,confRange(4:6)) % right conf
%           log(i).src_correct = 1;
%         end
%       elseif log(i).xcoord < 0.5 & log(i).xcoord >= 0 % left source
%         if ismember(log(i).conf,confRange(1:3)) % left conf
%           log(i).src_correct = -1;
%         elseif ismember(log(i).conf,confRange(4:6)) % right conf
%           log(i).src_correct = 0;
%         end
%       end
%     elseif log(i).isTarg == 0 % lure
%       log(i).src_correct = -1;
%     end
    
    % mark whether they got the source correct; this relative to R method
    % is described at the top of the function.
    %
    % if it was an R source (isTarg = 1) marked as 'R', it is correct;
    % if it was an R source (isTarg = 1) marked as 'L', it is wrong;
    % otherwise mark it as -1;
    if log(i).src_isTarg == 1 % right source
      if ismember(log(i).src_conf,confRange_src(4:6)) % right conf
        log(i).src_correct = 1; % hit
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct; % hit
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct; % hit
      elseif ismember(log(i).src_conf,confRange_src(1:3)) % left conf
        log(i).src_correct = 0; % miss
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct; % miss
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct; % miss
      end
    elseif log(i).src_isTarg == 0 % left source
      if ismember(log(i).src_conf,confRange_src(4:6)) % right conf
        log(i).src_correct = 0; % false alarm
        % put info in recognition response
        log(i-1).src_correct = log(i).src_correct; % false alarm
        % put info in test presentation
        log(i-2).src_correct = log(i).src_correct; % false alarm
      elseif ismember(log(i).src_conf,confRange_src(1:3)) % left conf
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
events = orderfields(events,{'subject','session','mstime','msoffset','trial','type','item','xcoord','ycoord','serialpos','rec_isTarg','rec_conf','rec_rt','rec_correct','src_isTarg','src_conf','src_rt','src_correct'});

% put the info from test into the study events
for i = 1:length(events)
  if strcmp(events(i).type,'STUDY_TARGET')
    testEvent = filterStruct(events,'ismember(type,varargin{1}) & ismember(item,varargin{2})',{'TARG_PRES'},{events(i).item});
    if length(testEvent) > 1
      error('Found multiple events when searching for study target %s!',events(i).item);
    else
      events(i).rec_isTarg = testEvent.rec_isTarg;
      events(i).rec_conf = testEvent.rec_conf;
      events(i).rec_rt = testEvent.rec_rt;
      events(i).rec_correct = testEvent.rec_correct;
      events(i).src_isTarg = testEvent.src_isTarg;
      events(i).src_conf = testEvent.src_conf;
      events(i).src_rt = testEvent.src_rt;
      events(i).src_correct = testEvent.src_correct;
    end
  end
end
