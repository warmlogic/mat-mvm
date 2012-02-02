function [conditions] = mm_ft_checkConditions(conditions,ana,condMethod)
%MM_FT_CHECKCONDITIONS Makes sure the conditions are set correctly, and
%make pairwise comparisons if desired
%
% [conditions] = mm_ft_checkConditions(conditions,ana,condMethod,data)
%
% inputs:
%   conditions = can be defined as one of three things:
%     1. your own comparisons: a cell containing one cell per comparison:
%        conditions = {{'A1','A2'},{'B1','B2'},{'C1','C2','C3'}};
%        in which case this function won't modify conditions
%
%     2. to make all pairwise comparisons within an event type, use:
%        conditions = {{''all_within_types''}};
%
%     3. to make all possible pairwise comparisons across event types, use:
%        conditions = {{''all_across_types''}};
%
%     NB: If you only have one event type, you can also use {'all'} in
%         conjunction with condMethod='pairwise', which will functionally
%         do the same thing as #2
%
%     NB: #3 could make a large amount of pairs, so watch out
%
%   ana = the analysis structure created in the experiment wrapper script;
%         ana.eventValues must contain the events; if applicable, events
%         should be separated by event type
%
%   condMethod = 'pairwise' for pairwise comparisons or 'single' for a flat
%                cell containing all condition names; 'check' is used to
%                check on a predefined array of event values and won't
%                create any comparisons
%
%     NB: condMethod='single' combined with conditions='all' or
%         'all_across_types' is the same as exper.eventValues (and
%         ana.eventValues when not using multiple event types);
%         condMethod='single' with conditions='all_within_types' is the
%         same as ana.eventValues
%
% output:
%  conditions = a cell array containing a cell array of event value
%               strings, with one cell array per type if desired

if nargin < 3 || isempty(condMethod)
  if ~iscell(conditions) && (~strcmp(conditions,'all') && ~strcmp(conditions,'all_across_types') && ~strcmp(conditions,'all_within_types'))
    condMethod = 'check';
%   elseif iscell(conditions) && (~strcmp(conditions,'all') && ~strcmp(conditions,'all_across_types') && ~strcmp(conditions,'all_within_types'))
%     condMethod = 'check';
  elseif iscell(conditions) && ~iscell(conditions{1}) && length(conditions) > 1
    condMethod = 'check';
  elseif iscell(conditions) && iscell(conditions{1}) && length(conditions{1}) > 1
    condMethod = 'check';
  else
    error('The variable ''condMethod'' was not defined; should be ''single'' or ''pairwise'' or ''check''.');
  end
  %fprintf('Guessing that condMethod=''%s''. Please see ''help mm_ft_checkConditions'' if this is not correct or for more info.',condMethod)
elseif nargin == 3 && ~ischar(condMethod)
  error('''condMethod'' should be a string');
end

if ~strcmp(condMethod,'single') && ~strcmp(condMethod,'pairwise') && ~strcmp(condMethod,'check')
  error('The variable ''condMethod'' was defined as ''%s''; should be ''single'' or ''pairwise'' or ''check''.',condMethod);
end

% make sure conditions is set inside a cell correctly
if ischar(conditions)
  conditions = {{conditions}};
elseif ~iscell(conditions{1})
  conditions = {conditions};
end

if strcmp(condMethod,'check') && length(conditions{1}) == 1 && (strcmp(conditions{1},'all') || strcmp(conditions{1},'all_across_types') || strcmp(conditions{1},'all_within_types'))
  error('condMethod=''check'' should not be used with conditions=''%s'' as no comparisons will be created.',conditions{1}{1});
%if strcmp(condMethod,'check') && length(conditions{1}) > 1
%  continue
end

if ~strcmp(condMethod,'check')
  % if we need to create the condition cells
  if strcmp(condMethod,'single')
    if length(ana.eventValues) == 1 && strcmp(conditions{1},'all')
      conditions = ana.eventValues;
    elseif length(ana.eventValues) > 1 && (strcmp(conditions{1},'all') || strcmp(conditions{1},'all_across_types'))
      conditions = cat(2,ana.eventValues{:});
    elseif length(ana.eventValues) > 1 && strcmp(conditions{1},'all_within_types')
      conditions = {};
      for typ = 1:length(ana.eventValues)
        conditions = cat(2,conditions,ana.eventValues(typ));
      end
    end
  elseif strcmp(condMethod,'pairwise')
    % see if we're comparing within or across types
    condWithinFlag = 0;
    condAcrossFlag = 0;
    for i = 1:length(conditions)
      % if only one entry in the cell, see what kind of combination to make
      if size(conditions{i},1) == 1 && size(conditions{i},2) == 1
        if strcmp(conditions{i},'all_within_types')
          condWithinFlag = 1;
        elseif strcmp(conditions{i},'all') && length(ana.eventValues) == 1
          condWithinFlag = 1;
        elseif strcmp(conditions{i},'all_across_types')
          condAcrossFlag = 1;
        else
          error('mm_ft_checkConditions:cfgConditionsFieldSizeOne','conditions was not defined properly.');
        end
      elseif size(conditions{i},1) == 1 && size(conditions{i},2) > 1
        condWithinFlag = 1;
      end
    end
    % if comparing within or across types, get the event combinations
    if condWithinFlag
      % PAIRWISE WITHIN EACH EVENT TYPE
      
      cond = {};
      for typ = 1:length(ana.eventValues)
        % find all the pairwise event combinations within this type
        cond = cat(1,cond,nchoosek(ana.eventValues{typ},2));
      end
      % preallocate
      conditions = {repmat({{}},1,size(cond,1))};
      % store each condition row as a separate cell
      for i = 1:size(cond,1)
        conditions{i} = cond(i,:);
      end
    elseif condAcrossFlag
      % PAIRWISE ACROSS ALL EVENT TYPES
      
      % get all of the events across all types together
      cond = cat(2,ana.eventValues{:});
      % find all the pairwise event combinations
      cond = nchoosek(cond,2);
      % preallocate
      conditions = {repmat({{}},1,size(cond,1))};
      % store each condition row as a separate cell
      for i = 1:size(cond,1)
        conditions{i} = cond(i,:);
      end
    end
  end
else
  % if we just need to check on the conditions that we've already made
  for i = 1:length(conditions)
    if size(conditions{i},1) > 1 && size(conditions{i},2) > 1
      error('mm_ft_checkConditions:cfgConditionsFieldLengthGT1','conditions was not defined properly; dim 2 length is greater than 1.');
    end
    % make sure the specified conditions are in ana.eventValues
    %if sum(~ismember(conditions{i},cat(2,ana.eventValues{:})))
    if sum(~ismember(cellflat(conditions{i}),cellflat(ana.eventValues))) ~= 0
      condNames = sprintf(repmat('%s ',1,length(conditions{i})),conditions{i}{:});
      anaEvValueNames = cellflat(ana.eventValues);
      anaEvValueNames = sprintf(repmat('%s ',1,length(cellflat(ana.eventValues))),anaEvValueNames{:});
      error('mm_ft_checkConditions:cfgConditionsNotIn_ana_eventValues','conditions %snot found in ana.eventValues: %s',condNames,anaEvValueNames);
    end
  end
end

end
