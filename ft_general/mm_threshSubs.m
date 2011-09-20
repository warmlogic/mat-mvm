function [exper] = mm_threshSubs(exper,ana,thresh)
%MM_THRESHSUBS Mark subjects as bad when they have fewer than X events
%
%   [exper] = mm_threshSubs(exper,ana,thresh)
%
%   If one event in ana.eventValues has fewer than 'thresh' trials, then
%   the subject is marked bad for all events values for that session
%
% Input:
%  exper  = exper struct (manually set bad subjects in exper.badBehSub={})
%  ana    = ana struct containing eventValues of interest
%  thresh = the event count threshold below
%
% Output:
%  exper.badSub         = overall, who's bad for any reason
%  exper.nTrials.lowNum = subjects with a low number of trials

%% setup

% make sure thresh is set; save it in exper struct
if nargin < 3
  error('Must specify ''thresh'' for the event count threshold');
else
  exper.nTrials.thresh = thresh;
end

if ~isfield(exper,'badBehSub')
  exper.badBehSub = {};
else
  exper.badBehSub = sort(exper.badBehSub(ismember(exper.badBehSub,exper.subjects)));
end

% get the events from ana; those are the ones we care about for analysis
events = cat(2,ana.eventValues{:});

% initialize to store who is below threshold
exper.nTrials.lowNum = false(length(exper.subjects),length(exper.sessions));

% figure out who is below threshold
for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    for evVal = 1:length(events)
      % compare nTrials to threshold
      if exper.nTrials.(events{evVal})(sub,ses) < exper.nTrials.thresh
        exper.nTrials.lowNum(sub,ses) = true;
      end
    end
  end
end

%% print some info

if length(exper.subjects{1}) > 7
  tabchar = '\t';
else
  tabchar = '';
end

fprintf('\n');

% print out the trial counts for each subject
fprintf('Subject%s%s\n',sprintf(tabchar),sprintf(repmat('\t%s',1,length(events)),events{:}));
for ses = 1:length(exper.sessions)
  for sub = 1:length(exper.subjects)
    subStr = exper.subjects{sub};
    for evVal = 1:length(events)
      if length(events{evVal}) > 7
        tabchar_ev = '\t';
      else
        tabchar_ev = '';
      end
      subStr = cat(2,subStr,sprintf('\t%d%s',exper.nTrials.(events{evVal})(sub,ses),sprintf(tabchar_ev)));
    end
    fprintf('%s\n',subStr);
  end
end

% display who we're rejecting because of low trial counts
if sum(exper.nTrials.lowNum) > 0
  fprintf('\nSubjects with low trial counts:\n');
  for i = 1:length(exper.subjects)
    if exper.nTrials.lowNum(i)
      fprintf('%s\n',exper.subjects{i});
    end
  end
else
  fprintf('\nNo subjects with low trial counts.\n');
end

% print pre-defined subjects who we want to exclude
if ~isempty(exper.badBehSub)
  fprintf('\nSubjects with bad behavior:\n');
  for i = 1:length(exper.badBehSub)
    fprintf('%s\n',exper.badBehSub{i});
  end
else
  fprintf('\nNo subjects with bad behavior.\n');
end

% mark bad subjects and sessions
exper.badSub = zeros(length(exper.subjects),length(exper.sessions));
for ses = 1:length(exper.sessions)
  exper.badSub(:,ses) = sign(sum([ismember(exper.subjects,exper.badBehSub) exper.nTrials.lowNum(:,ses)],2));
end
exper.badSub = logical(exper.badSub);

% print a summary
fprintf('\nNumber of events included in EEG analyses (%d subjects; threshold: %d events):\n',sum(~exper.badSub),exper.nTrials.thresh);
for evVal = 1:length(events)
  meanTrials = mean(exper.nTrials.(events{evVal})(~exper.badSub));
  semTrials = std(exper.nTrials.(events{evVal})(~exper.badSub),0,1)/sqrt(sum(~exper.badSub));
  stdTrials = std(exper.nTrials.(events{evVal})(~exper.badSub),0,1);
  fprintf('%s:\tM=%.3f,\tSEM=%.3f,\tSD=%.3f\n',events{evVal},meanTrials,semTrials,stdTrials);
end

end

% %% figure out who to reject
% 
% for evVal = 1:length(events)
%   for sub = 1:length(exper.subjects)
%     for ses = 1:length(exper.sessions)
%       exper.nTrials.(events{evVal})(sub,ses) = size(ft_findcfg(data.(events{evVal}).sub(sub).ses(ses).data.cfg,'trl'),1);
%     end
%   end
% end
% 
% if exist('data_freq','var')
%   for evVal = 1:length(events)
%     for sub = 1:length(exper.subjects)
%       for ses = 1:length(exper.sessions)
%         exper.nTrials.(events{evVal})(sub,ses) = size(ft_findcfg(data_freq.(events{evVal}).sub(sub).ses(ses).data.cfg,'trl'),1);
%       end
%     end
%   end
% elseif exist('ga_freq','var')
%   % if we need to recreate the bad subjects list from a grand average file
%   goodSub = cell(size(ga_freq.(events{1}).cfg.previous,2),1);
%   for evVal = 1:length(events)
%     for sub = 1:size(ga_freq.(events{1}).cfg.previous,2)
%       datafile = ft_findcfg(ga_freq.(events{1}).cfg.previous{sub},'datafile');
%       [pathstr,datafile] = fileparts(datafile);
%       dividers = strfind(datafile,'_');
%       if isempty(dividers)
%         dividers = strfind(datafile,' ');
%       end
%       goodSub{sub} = datafile(strfind(datafile,exper.name):dividers(1)-1);
%       exper.nTrials.(events{evVal})(sub,ses) = size(ft_findcfg(ga_freq.(events{evVal}).cfg.previous{sub},'trl'),1);
%     end
%   end
%   exper.badSub = ~ismember(exper.subjects,goodSub);
%   [exper.subjects,subInd] = sort(cat(1,goodSub,exper.subjects(exper.badSub)));
%   for evVal = 1:length(events)
%     exper.nTrials.(events{evVal}) = exper.nTrials.(events{evVal})(subInd);
%   end
%   % add the extra event values onto the end of eventValues
%   if ~isempty(exper.eventValuesExtra)
%     for nVal = 1:length(exper.eventValuesExtra.newValue)
%       if ~ismember(exper.eventValuesExtra.newValue{nVal},events)
%         events = cat(2,events,exper.eventValuesExtra.newValue{nVal});
%       end
%     end
%     events = sort(events);
%   end
% end

