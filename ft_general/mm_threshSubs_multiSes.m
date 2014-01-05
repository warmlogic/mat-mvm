function [exper,ana] = mm_threshSubs_multiSes(exper,ana,thresh,commonGoodChan)
%MM_THRESHSUBS_MULTISES Mark subjects as bad when they have fewer than X events
%
%   [exper,ana] = mm_threshSubs_multiSes(exper,ana,thresh,commonGoodChan)
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
%
% Output:
%  exper.badSub                  = overall, who's bad for any reason
%  exper.nTrials.(sesStr).lowNum = subjects with a low number of trials
%  ana.elecCommon                = the common good channels across good
%                                  subjects (only if commonGoodChan=true)
%

%% setup

% make sure thresh is set; save it in exper struct
if nargin < 4
  commonGoodChan = false;
  if nargin < 3
    error('Must specify ''thresh'' for the event count threshold');
  end
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
  
  % get the events from ana; those are the ones we care about for analysis
  %eventValues = ana.eventValues{ses};
  
  % initialize to store who is below threshold
  exper.nTrials.(exper.sesStr{ses}).lowNum = false(length(exper.subjects),length(exper.sessions));
  
  % figure out who is below threshold
  for sub = 1:length(exper.subjects)
    %for ses = 1:length(exper.sessions)
    for evVal = 1:length(ana.eventValues{ses})
      % compare nTrials to threshold
      if exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal})(sub) < exper.nTrials.(exper.sesStr{ses}).thresh
        exper.nTrials.(exper.sesStr{ses}).lowNum(sub) = true;
      end
    end
    %end
  end
  
end

%% print some info

if length(exper.subjects{1}) > 7
  tabchar_sub = '\t';
else
  tabchar_sub = '';
end

fprintf('\n');

% print out the trial counts for each subject
for ses = 1:length(exper.sessions)
  fprintf('========================\n');
  fprintf('Session: %s\n',exper.sesStr{ses});
  fprintf('========================\n');
  fprintf('Subject%s%s\n',sprintf(tabchar_sub),sprintf(repmat('\t%s',1,length(ana.eventValues{ses})),ana.eventValues{ses}{:}));
  for sub = 1:length(exper.subjects)
    subStr = exper.subjects{sub};
    for evVal = 1:length(ana.eventValues{ses})
      if length(ana.eventValues{ses}{evVal}) > 7
        tabchar_ev = '\t';
      else
        tabchar_ev = '';
      end
      subStr = cat(2,subStr,sprintf('\t%d%s',exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal})(sub),sprintf(tabchar_ev)));
    end
    fprintf('%s\n',subStr);
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

if sum(~exper.badSub) > 0
  % print a summary
  fprintf('\nNumber of events included in EEG analyses (%d subjects; threshold: %d events):\n',sum(~exper.badSub),thresh);
  for ses = 1:length(exper.sessions)
    for evVal = 1:length(ana.eventValues{ses})
      meanTrials = mean(exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal})(~exper.badSub));
      semTrials = std(exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal})(~exper.badSub),0,1)/sqrt(sum(~exper.badSub));
      stdTrials = std(exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal})(~exper.badSub),0,1);
      fprintf('%s:\tM=%.3f,\tSEM=%.3f,\tSD=%.3f\n',ana.eventValues{ses}{evVal},meanTrials,semTrials,stdTrials);
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
