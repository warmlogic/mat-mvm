function mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,data)
%MM_FT_PLOTTFR make (and save) Time-Frequency plots
%
%   mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,data)
%
% Inputs:
%   cfg_ft: parameters passed into the FT plotting function
%
%   cfg_plot.ftFxn      = FieldTrip plotting function to use. Supported
%                         functions: ft_singleplotTFR, ft_topoplotTFR, and
%                         ft_multiplotTFR
%   cfg_plot.conditions = Cell array containing cells of pairwise
%                         comparisons; Can be used for comparing a subset
%                         of events within a type.
%                         e.g., {{'T1a','T1c'}, {'T2a','T2c'}}, or it can
%                         be {{'all_within_types'}} or
%                         {{'all_across_types'}} to automatically create
%                         pairwise comparisons of event values. See
%                         MM_FT_CHECKCONDITIONS for more details.
%   cfg_plot.plotTitle  = 1 or 0. Whether to plot the title.
%   cfg_plot.subplot    = 1 or 0. Whether to make a subplot. cfg_ft.xlim
%                         can be a range of time values, otherwise 50ms
%                         steps between min and max. ft_topoplotER only.
%   cfg_plot.numCols    = If subplot == 1, the number of columns to plot
%   files.saveFigs      = 1 or 0. Whether to save the figures.
%
%   data                = output from ft_freqgrandaverage
%
% See also:
%   MM_FT_CHECKCONDITIONS

if ~isfield(cfg_ft,'parameter')
  error('Must specify cfg_ft.parameter, denoting the data to plot (e.g., ''powspctrm'' or ''cohspctrm'')');
end

if ~isfield(cfg_plot,'plotTitle')
  cfg_plot.plotTitle = 0;
end

cfg_plot.type = strrep(strrep(cfg_plot.ftFxn,'ft_',''),'plotTFR','');

if strcmp(cfg_plot.type,'multi') || strcmp(cfg_plot.type,'topo')
  % need a layout if doing a topo or multi plot
  if isfield(ana,'elec')
    cfg_ft.layout = ft_prepare_layout([],ana);
  else
    error('''ana'' struct must have ''elec'' field');
  end
  
  if ~isfield(cfg_plot,'roi')
    % use all channels in a topo or multi plot
    cfg_plot.roi = {'all'};
  end
  
  if strcmp(cfg_plot.type,'topo')
    if isfield(cfg_ft,'showlabels')
      % not allowed
      cfg_ft = rmfield(cfg_ft,'showlabels');
    end
    if isfield(cfg_ft,'markerfontsize')
      cfg_ft.markerfontsize = 9;
    end
  end
end

% make sure conditions are set correctly
if ~isfield(cfg_plot,'condMethod')
  if ~iscell(cfg_plot.conditions) && (strcmp(cfg_plot.conditions,'all') || strcmp(cfg_plot.conditions,'all_across_types') || strcmp(cfg_plot.conditions,'all_within_types'))
    cfg_plot.condMethod = 'single';
  elseif iscell(cfg_plot.conditions) && ~iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
    cfg_plot.condMethod = 'single';
  elseif iscell(cfg_plot.conditions) && iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions{1}) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
    cfg_plot.condMethod = 'single';
  else
    cfg_plot.condMethod = [];
  end
end

%cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);

% temporary hack for plotting a single subject and 1 evVal, because the
% data struct will have a field titled parameter and won't have the typical
% data.evVal.sub.ses.data structure
if ~isfield(data,cfg_ft.parameter)
  cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);
elseif isfield(data,cfg_ft.parameter) && length(cfg_plot.conditions{:}) > 1
  error('Cannot have more than one condition if data is only one evVal');
end

% make sure conditions are set up for the for loop
if ~isfield(cfg_plot,'types')
  cfg_plot.types = repmat({''},size(cfg_plot.conditions));
end

% set the channel information
if ~isfield(cfg_plot,'roi')
  error('Must specify either ROI names or channel names in cfg_plot.roi');
elseif isfield(cfg_plot,'roi')
  if ismember(cfg_plot.roi,ana.elecGroupsStr)
    % if it's in the predefined ROIs, get the channel numbers
    if strcmp(cfg_plot.type,'topo')
      cfg_ft.highlight = 'on';
      cfg_ft.highlightsize = 10;
      cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
    else
      cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
    end
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:});
  else
    % otherwise it should be the channel number(s) or 'all'
    if ~iscell(cfg_plot.roi)
      cfg_plot.roi = {cfg_plot.roi};
    end
    
    if strcmp(cfg_plot.type,'topo')
      if ~strcmp(cfg_plot.roi,'all')
        cfg_ft.highlight = 'on';
        cfg_ft.highlightsize = 10;
        cfg_ft.highlightchannel = cfg_plot.roi;
      end
    else
      cfg_ft.channel = cfg_plot.roi;
    end
    
    % set the string for the filename
    if isfield(cfg_ft,'cohrefchannel')
      cfg_plot.chan_str = [cfg_ft.cohrefchannel,'-',sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:})];
    else
      cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:});
    end
  end
end

% time - get this info for the figure name
if isfield(cfg_ft,'xlim')
  if strcmp(cfg_ft.xlim,'maxmin')
    cfg_ft.xlim = [min(data.(cfg_plot.conditions{1}{1}).time) max(data.(cfg_plot.conditions{1}{1}).time)];
  end
else
  cfg_ft.xlim = [min(data.(cfg_plot.conditions{1}{1}).time) max(data.(cfg_plot.conditions{1}{1}).time)];
end

% set parameters for the subplot
if isfield(cfg_plot,'subplot')
  if cfg_plot.subplot
    if ~strcmp(cfg_plot.type,'topo')
      fprintf('Subplot only works with topoplot! Changing to non-subplot.\n');
      cfg_plot.subplot = 0;
    elseif strcmp(cfg_plot.type,'topo')
      if length(cfg_ft.xlim) > 2
        % predefined time windows
        cfg_plot.timeS = cfg_ft.xlim;
      else
        % default: 50 ms time windows
        cfg_plot.timeS = (cfg_ft.xlim(1):0.05:cfg_ft.xlim(2));
      end
      
      if ~isfield(cfg_plot,'numCols')
        cfg_plot.numCols = 5;
      end
      if (length(cfg_plot.timeS)-1) < cfg_plot.numCols
        cfg_plot.numCols = (length(cfg_plot.timeS)-1);
      end
      cfg_plot.numRows = ceil((length(cfg_plot.timeS)-1)/cfg_plot.numCols);
      
      % a few settings to make the graphs viewable
      if ~isfield(cfg_ft,'comment')
        cfg_ft.comment = 'xlim';
      end
      cfg_ft.commentpos = 'title';
      cfg_ft.colorbar = 'no';
      cfg_ft.marker = 'on';
      cfg_ft.fontsize = 10;
      if isfield(cfg_ft,'markerfontsize')
        cfg_ft = rmfield(cfg_ft,'markerfontsize');
      end
      cfg_plot.plotTitle = 0;
    end
  end
else
  cfg_plot.subplot = 0;
end

% freq - get this info for the figure name
if isfield(cfg_ft,'ylim')
  if strcmp(cfg_ft.ylim,'maxmin')
    cfg_ft.ylim = [min(data.(cfg_plot.conditions{1}{1}).freq) max(data.(cfg_plot.conditions{1}{1}).freq)];
  end
else
  cfg_ft.ylim = [min(data.(cfg_plot.conditions{1}{1}).freq) max(data.(cfg_plot.conditions{1}{1}).freq)];
end

% find the indices for the times and frequencies that we want
%cfg_plot.timesel = find(data.(cfg_plot.conditions{1}{1}).time >= cfg_ft.xlim(1) & data.(cfg_plot.conditions{1}{1}).time <= cfg_ft.xlim(2));
%cfg_plot.freqsel = find(data.(cfg_plot.conditions{1}{1}).freq >= cfg_ft.ylim(1) & data.(cfg_plot.conditions{1}{1}).freq <= cfg_ft.ylim(2));

for typ = 1:length(cfg_plot.conditions)
  for evVal = 1:length(cfg_plot.conditions{typ})
    figure
    
    if cfg_plot.subplot
      for k = 1:length(cfg_plot.timeS)-1
        subplot(cfg_plot.numRows,cfg_plot.numCols,k);
        cfg_ft.xlim = [cfg_plot.timeS(k) cfg_plot.timeS(k+1)];
        if isfield(data,cfg_ft.parameter)
          % temporary hack for plotting a single subject and 1 evVal,
          % because the data struct will have a field titled parameter and
          % won't have the typical data.evVal.sub.ses.data structure
          feval(str2func(cfg_plot.ftFxn),cfg_ft,data);
        else
          feval(str2func(cfg_plot.ftFxn),cfg_ft,data.(cfg_plot.conditions{typ}{evVal}));
        end
      end
      % reset the xlim
      cfg_ft.xlim = [cfg_plot.timeS(1) cfg_plot.timeS(end)];
    else
      feval(str2func(cfg_plot.ftFxn),cfg_ft,data.(cfg_plot.conditions{typ}{evVal}));
    end
    
    if ~isempty(cfg_plot.types{typ})
      set(gcf,'Name',sprintf('%s, %s, %s, %.1f--%.1f Hz, %d--%d ms',cfg_plot.types{typ},strrep(cfg_plot.conditions{typ}{evVal},'_',''),strrep(cfg_plot.chan_str,'_',' '),cfg_ft.ylim(1),cfg_ft.ylim(2),round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000)))
    else
      set(gcf,'Name',sprintf('%s, %s, %.1f--%.1f Hz, %d--%d ms',strrep(cfg_plot.conditions{typ}{evVal},'_',''),strrep(cfg_plot.chan_str,'_',' '),cfg_ft.ylim(1),cfg_ft.ylim(2),round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000)))
    end
    
    if strcmp(cfg_ft.colorbar,'yes')
      cfg_plot.colorbar_str = '_cb';
    else
      cfg_plot.colorbar_str = '';
    end
    if cfg_plot.subplot
      cfg_plot.subplot_str = '_subplot';
    else
      cfg_plot.subplot_str = '';
    end
    if cfg_plot.plotTitle
      %title(sprintf('%s - %s, %.1f--%.1f Hz, %.1f--%.1f s',cfg_plot.condNames{1},cfg_plot.condNames{2},cfg_ft.ylim(1),cfg_ft.ylim(2),cfg_ft.xlim(1),cfg_ft.xlim(2)));
      title(sprintf('%s, %.1f--%.1f Hz, %.1f--%.1f s',strrep(cfg_plot.conditions{typ}{evVal},'_',''),cfg_ft.ylim(1),cfg_ft.ylim(2),cfg_ft.xlim(1),cfg_ft.xlim(2)));
      cfg_plot.title_str = '_title';
    else
      cfg_plot.title_str = '';
    end
    if ~isfield(files,'figFontName')
      files.figFontName = 'Helvetica';
    end
    publishfig(gcf,~cfg_plot.plotTitle,[],[],files.figFontName);
      
    if files.saveFigs
      if ~isempty(cfg_plot.types{typ})
        cfg_plot.figfilename = sprintf('tfr_%s_ga_%s_%s_%s%d_%d_%d_%d%s%s%s',cfg_plot.type,cfg_plot.types{typ},cfg_plot.conditions{typ}{evVal},cfg_plot.chan_str,round(cfg_ft.ylim(1)),round(cfg_ft.ylim(2)),round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.colorbar_str,cfg_plot.subplot_str,cfg_plot.title_str);
      else
        cfg_plot.figfilename = sprintf('tfr_%s_ga_%s_%s%d_%d_%d_%d%s%s%s',cfg_plot.type,cfg_plot.conditions{typ}{evVal},cfg_plot.chan_str,round(cfg_ft.ylim(1)),round(cfg_ft.ylim(2)),round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.colorbar_str,cfg_plot.subplot_str,cfg_plot.title_str);
      end
      dirs.saveDirFigsTFR = fullfile(dirs.saveDirFigs,['tfr_',cfg_plot.type]);
      if ~exist(dirs.saveDirFigsTFR,'dir')
        mkdir(dirs.saveDirFigsTFR)
      end
      
      if strcmp(files.figPrintFormat(1:2),'-d')
        files.figPrintFormat = files.figPrintFormat(3:end);
      end
      if ~isfield(files,'figPrintRes')
        files.figPrintRes = 150;
      end
      print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsTFR,cfg_plot.figfilename));
    end
    
  end % evVal
end % typ

end
