function [exper,ana] = mm_threshSubs_multiSes(exper,ana,thresh,commonGoodChan,printDirection,eventValues)
%MM_THRESHSUBS_MULTISES Mark subjects as bad when they have fewer than X events
%
%   [exper,ana] = mm_threshSubs_multiSes(exper,ana,thresh,commonGoodChan,printDirection)
%
%   If one event in ana.eventValues has fewer than 'thresh' trials, then
%   the subject is marked bad for all events values for that session
%
% Input:
%  exper  = exper struct (manually set bad subjects in exper.badBehSub={};
%           the outer cell contains a cell for each session)
%  ana    = ana struct containing eventValues of interest, for each session
%  thresh = the event count threshold below (<) which a subject will be
%           counted as having a too low trial count
%  commonGoodChan = true/false, whether to find the common good channels
%                   across the good subjects
%  printDirection = how to print out trial counts ('horiz' or 'vert')
%
% Output:
%  exper.badSub                  = overall, who's bad for any reason
%  exper.nTrials.(sesStr).lowNum = subjects with a low number of trials
%  ana.elecCommon                = the common good channels across good
%                                  subjects (only if commonGoodChan=true)
%

%% setup

% whether to exclude a bad subject in any session from all sessions,
% regardless of whether they were good in the others
collapseSessions = false;


if nargin < 6
  eventValues = [];
  if nargin < 5
    printDirection = [];
    if nargin < 4
      commonGoodChan = [];
      if nargin < 3
        thresh = [];
      end
    end
  end
end

if ~exist('eventValues','var') || isempty(eventValues)
  eventValues = ana.eventValues;
end

if ~exist('printDirection','var') || isempty(printDirection)
  printDirection = 'horiz';
end

if ~exist('commonGoodChan','var') || isempty(commonGoodChan)
  commonGoodChan = false;
end

% make sure thresh is set; save it in exper struct
if ~exist('thresh','var') || isempty(thresh)
  error('Must specify ''thresh'' for the event count threshold');
end

if ~isfield(exper,'badBehSub')
  exper.badBehSub = cell(1,length(exper.sesStr));
end

for ses = 1:length(exper.sesStr)
  exper.nTrials.(exper.sesStr{ses}).thresh = thresh;
  
  if isempty(exper.badBehSub{ses})
    exper.badBehSub{ses} = {};
  else
    exper.badBehSub{ses} = sort(exper.badBehSub{ses}(ismember(exper.badBehSub{ses},exper.subjects)));
  end
  
  % initialize to store who is below threshold
  exper.nTrials.(exper.sesStr{ses}).lowNum = false(length(exper.subjects),1);
  
  % figure out who is below threshold
  for sub = 1:length(exper.subjects)
    %for ses = 1:length(exper.sessions)
    for evVal = 1:length(eventValues{ses})
      for es = 1:length(eventValues{ses}{evVal})
        % compare nTrials to threshold
        if exper.nTrials.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es})(sub) < exper.nTrials.(exper.sesStr{ses}).thresh
          exper.nTrials.(exper.sesStr{ses}).lowNum(sub) = true;
        end
      end
    end
    %end
  end
  
end

%% print some info

tabLimit = 7;

tabchar_sub = '';
if length(exper.subjects{1}) > tabLimit
  for i = 1:floor(length(exper.subjects{1}) / tabLimit)
    tabchar_sub = cat(2,tabchar_sub,'\t');
  end
end

fprintf('\n');

% print out the trial counts for each subject
for ses = 1:length(exper.sessions)
  fprintf('========================\n');
  fprintf('Session: %s\n',exper.sesStr{ses});
  fprintf('========================\n');
  for evVal = 1:length(eventValues{ses})
    if strcmp(printDirection,'horiz')
    fprintf('Subject%s%s\n',sprintf(tabchar_sub),sprintf(repmat('\t%s',1,length(eventValues{ses}{evVal})),eventValues{ses}{evVal}{:}));
    elseif strcmp(printDirection,'vert')
      fprintf('Subject\n');
    end
    for sub = 1:length(exper.subjects)
      subStr = exper.subjects{sub};
      for es = 1:length(eventValues{ses}{evVal})
        if strcmp(printDirection,'horiz')
          tabchar_ev = '';
          if length(eventValues{ses}{evVal}{es}) > tabLimit
            for i = 1:floor(length(eventValues{ses}{evVal}{es}) / tabLimit)
              tabchar_ev = cat(2,tabchar_ev,'\t');
            end
          end
          
          subStr = cat(2,subStr,sprintf('\t%d%s',exper.nTrials.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es})(sub),sprintf(tabchar_ev)));
        elseif strcmp(printDirection,'vert')
          if es == 1
            subStr = sprintf('%s\n',subStr);
          end
          subStr = cat(2,subStr,sprintf('\t%s: %d\n',eventValues{ses}{evVal}{es},exper.nTrials.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es})(sub)));
        end
      end
      fprintf('%s\n',subStr);
    end
  end
  
  % display who we're rejecting because of low trial counts
  if sum(exper.nTrials.(exper.sesStr{ses}).lowNum) > 0
    fprintf('\nSubjects with low trial counts:\n');
    for i = 1:length(exper.subjects)
      if exper.nTrials.(exper.sesStr{ses}).lowNum(i)
        fprintf('%s\n',exper.subjects{i});
      end
    end
  else
    fprintf('\nNo subjects with low trial counts.\n');
  end
  
  % print pre-defined subjects who we want to exclude
  if ~isempty(exper.badBehSub{ses})
    fprintf('\nSubjects with bad behavior:\n');
    for i = 1:length(exper.badBehSub{ses})
      fprintf('%s\n',exper.badBehSub{ses}{i});
    end
  else
    fprintf('\nNo subjects with bad behavior.\n');
  end
end

% mark bad subjects and sessions
exper.badSub = zeros(length(exper.subjects),length(exper.sessions));
for ses = 1:length(exper.sessions)
  exper.badSub(:,ses) = sign(sum([ismember(exper.subjects,exper.badBehSub{ses}) exper.nTrials.(exper.sesStr{ses}).lowNum(:)],2));
end
exper.badSub = logical(exper.badSub);

if collapseSessions
  exper.badSub = logical(sum(exper.badSub,2));
  exper.badSub = repmat(exper.badSub,1,length(exper.sessions));
end

if sum(~exper.badSub) > 0
  % print a summary
  for ses = 1:length(exper.sessions)
    fprintf('\nSession %s:\n\tNumber of events included in EEG analyses (%d subjects; threshold: %d events):\n',exper.sesStr{ses},sum(~exper.badSub(:,ses)),thresh);
    for evVal = 1:length(eventValues{ses})
      for es = 1:length(eventValues{ses}{evVal})
        meanTrials = mean(exper.nTrials.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es})(~exper.badSub(:,ses)));
        semTrials = std(exper.nTrials.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es})(~exper.badSub(:,ses)),0,1)/sqrt(sum(~exper.badSub(:,ses)));
        stdTrials = std(exper.nTrials.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es})(~exper.badSub(:,ses)),0,1);
        fprintf('%s:\tM=%.3f,\tSEM=%.3f,\tSD=%.3f\n',eventValues{ses}{evVal}{es},meanTrials,semTrials,stdTrials);
      end
    end
  end
else
  fprintf('\n');
  warning('No good subjects remain when threshold=%d trials!',thresh);
end

% find the channels that are common across subjects; start with full set
if commonGoodChan
  cfg_com = [];
  cfg_com.channel = ana.elec.label;
  fprintf('Finding the common good channels for good subjects');
  for sub = 1:length(exper.subjects)
    fprintf('.');
    for ses = 1:length(exper.sessions)
      if ~exper.badSub(sub,ses)
        % this subject's channels
        subChan = eval(sprintf('{''all'' ''-Fid*''%s};',sprintf(repmat(' ''-E%d''',1,length(exper.badChan{ses}{sub})),exper.badChan{ses}{sub})));
        % whittle down to common channels
        cfg_com.channel = ft_channelselection(subChan,cfg_com.channel);
      end
    end
  end
  fprintf('Done.\n');
  % save the common good channels
  ana.elecCommon = cfg_com.channel;
else
  % return all normal channels
  ana.elecCommon = ft_channelselection({'all' '-Fid*'},ana.elec.label);
end

end % function
