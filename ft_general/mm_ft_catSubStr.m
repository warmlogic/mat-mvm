function [ana_str] = mm_ft_catSubStr(cfg,exper,str_format)
%MM_FT_CATSUBSTR Concatenate strings of subject data for input to FieldTrip
%functions
%
% [ana_str] = mm_ft_catSubStr(cfg,exper,str_format)
%
% input:
%   cfg.conditions    = which condition fields to create in ana_str; should
%                       be a cell array of event values as strings
%   cfg.data_str      = the data struct name (e.g., 'data_tla',
%                       'data_freq', or just 'data')
%   cfg.is_ga         = grand average data or individual subject data?
%                       (1 or 0)
%   cfg.excludeBadSub = 1 or 0; excludes bad subjects by default
%   exper             = exper struct with subject, session, and badSub info
%   str_format        = 1 (new standard multisession format) or 2 (old format)
%
%
% output:
%   ana_str is a struct with fields corresponding to event values in
%   cfg.conditions. Each field contains a string to access the data of all
%   subjects for that event value. It is meant to serve as the input for
%   FieldTrip functions like ft_timelockstatistics, ft_freqstatistics, and
%   the grand average functions.
%   e.g., ana_str.CR == 'data_tla.CR.sub(1).ses(1).data,data_tla.CR.sub(2).ses(1).data'
%

if ~exist('str_format','var') || isempty(str_format)
  str_format = 2;
end

% make sure the conditions are stored in a cell
if ~iscell(cfg.conditions) && ischar(cfg.conditions)
  cfg.conditions = {cfg.conditions};
end

if nargin < 2
  if ~cfg.is_ga
    error('You must include the exper struct if using individual subject data.');
  end
end

if ~cfg.is_ga
  % exclude the bad subjects by default
  if ~isfield(cfg,'excludeBadSub')
    cfg.excludeBadSub = 1;
  end
end

if str_format == 1
  if iscell(cfg.conditions) && ~iscell(cfg.conditions{1})
    cfg.conditions = {cfg.conditions};
  end
  
  if ~cfg.is_ga
    % if we have individual subject data...
    
    % initialize to see if we've added the first subject
    firstOneDone = 0;
    
    % go through subjects, add strings for the ones we want
    for ses = 1:length(exper.sessions)
      for sub = 1:length(exper.subjects)
        if exper.badSub(sub,ses)
          if cfg.excludeBadSub
            % skip this subject if they're bad
            fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
            continue
          else
            % keep them in if they're bad
            fprintf('Including bad subject: %s\n',exper.subjects{sub});
          end
        else
          if ~firstOneDone
            % add the first subject; the string is formatetted differently
            for evVal = 1:length(cfg.conditions{ses})
              ana_str.(cfg.conditions{ses}{evVal}){ses} = sprintf('%s.%s.%s.sub(%d).data',cfg.data_str,exper.sesStr{ses},cfg.conditions{ses}{evVal},sub);
            end
            firstOneDone = 1;
          else
            for evVal = 1:length(cfg.conditions{ses})
              ana_str.(cfg.conditions{ses}{evVal}){ses} = sprintf('%s,%s.%s.%s.sub(%d).data',ana_str.(cfg.conditions{ses}{evVal}){ses},cfg.data_str,exper.sesStr{ses},cfg.conditions{ses}{evVal},sub);
            end
          end
        end
      end % sub
    end % ses
  else
    % if we have grand average data, then we just need the event value names
    for ses = 1:length(exper.sessions)
      ana_str = sprintf('%s.%s.%s',cfg.data_str,exper.sesStr{ses},cfg.conditions{ses}{1});
      if length(cfg.conditions{ses}) > 1
        for evVal = 2:length(cfg.conditions{ses})
          ana_str = sprintf('%s,%s.%s.%s',ana_str,cfg.data_str,exper.sesStr{ses},cfg.conditions{ses}{evVal});
        end
      end
    end
  end
elseif str_format == 2
  
  % concatenate the types together if they were separate
  if iscell(cfg.conditions) && iscell(cfg.conditions{1})
    cfg.conditions = cat(2,cfg.conditions{:});
  end
  
  if ~cfg.is_ga
    % if we have individual subject data...
    
    % initialize to see if we've added the first subject
    firstOneDone = 0;
    
    % go through subjects, add strings for the ones we want
    for sub = 1:length(exper.subjects)
      for ses = 1:length(exper.sessions)
        if exper.badSub(sub,ses)
          if cfg.excludeBadSub
            % skip this subject if they're bad
            fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
            continue
          else
            % keep them in if they're bad
            fprintf('Including bad subject: %s\n',exper.subjects{sub});
          end
        else
          if ~firstOneDone
            % add the first subject; the string is formatetted differently
            for evVal = 1:length(cfg.conditions)
              ana_str.(cfg.conditions{evVal}){ses} = sprintf('%s.%s.sub(%d).ses(%d).data',cfg.data_str,cfg.conditions{evVal},sub,ses);
            end
            firstOneDone = 1;
          else
            for evVal = 1:length(cfg.conditions)
              ana_str.(cfg.conditions{evVal}){ses} = sprintf('%s,%s.%s.sub(%d).ses(%d).data',ana_str.(cfg.conditions{evVal}){ses},cfg.data_str,cfg.conditions{evVal},sub,ses);
            end
          end
        end
      end
    end
  else
    % if we have grand average data, then we just need the event value names
    ana_str = sprintf('%s.%s',cfg.data_str,cfg.conditions{1});
    if length(cfg.conditions) > 1
      for evVal = 2:length(cfg.conditions)
        ana_str = sprintf('%s,%s.%s',ana_str,cfg.data_str,cfg.conditions{evVal});
      end
    end
  end
end % str_format

end % function
