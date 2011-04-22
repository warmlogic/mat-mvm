function [ana_str] = mm_ft_catSubStr(cfg,exper)
%MM_FT_CATSUBSTR Concatenate strings of subject data for input to FieldTrip
%functions
%
% [ana_str] = mm_ft_catSubStr(cfg,exper)
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
%
% output:
%   ana_str is a struct with fields corresponding to event values in
%   cfg.conditions. Each field contains a string to access the data of all
%   subjects for that event value. It is meant to serve as the input for
%   FieldTrip functions like ft_timelockstatistics, ft_freqstatistics, and
%   the grand average functions.
%   e.g., ana_str.CR == 'data_tla.CR.sub(1).ses(1).data,data_tla.CR.sub(2).ses(1).data'
%

% make sure the conditions are stored in a cell
if ~iscell(cfg.conditions) && ischar(cfg.conditions)
  cfg.conditions = {cfg.conditions};
end
% concatenate the types together if they were separate
if iscell(cfg.conditions) && iscell(cfg.conditions{1})
  cfg.conditions = cat(2,cfg.conditions{:});
end

if ~isfield(cfg,'excludeBadSub')
  cfg.excludeBadSub = 1;
end

if nargin < 2
  if ~cfg.is_ga
    error('Need to include exper struct if using individual subject data.');
  end
end

if ~cfg.is_ga
  firstOneDone = 0;
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      if exper.badSub(sub,ses)
        if cfg.excludeBadSub
          fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
          continue
        else
          fprintf('Including bad subject: %s\n',exper.subjects{sub});
        end
      else
        if ~firstOneDone
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
  ana_str = sprintf('%s.%s',cfg.data_str,cfg.conditions{1});
  if length(cfg.conditions) > 1
    for evVal = 2:length(cfg.conditions)
      ana_str = sprintf('%s,%s.%s',ana_str,cfg.data_str,cfg.conditions{evVal});
    end
  end
end

end
