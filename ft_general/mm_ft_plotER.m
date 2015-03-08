function mm_ft_plotER(cfg_ft,cfg_plot,exper,ana,files,dirs,data,sesNum)
%MM_FT_PLOTER make (and save) event-related data plots
%
%   mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,data)
%
% Inputs:
%   cfg_ft: parameters passed into the FT plotting function
%
%   cfg_plot.ftFxn      = FieldTrip plotting function to use. Supported
%                         functions: ft_singleplotER, ft_topoplotER, and
%                         ft_multiplotER
%   cfg_plot.conditions = Cell array containing cells of pairwise
%                         comparisons; Can be used for comparing a subset
%                         of events within a type.
%                         e.g., {{'T1a','T1c'}, {'T2a','T2c'}}, or it can
%                         be {{'all_within_types'}} or
%                         {{'all_across_types'}} to automatically create
%                         pairwise comparisons of event values. See
%                         MM_FT_CHECKCONDITIONS for more details.
%   cfg_plot.plotTitle  = 1 or 0*. Whether to plot the title.
%   cfg_plot.plotLegend = 1 or 0*. Whether to plot the legend.
%   cfg_plot.subplot    = 1 or 0*. Whether to make a subplot. cfg_ft.xlim
%                         can be a range of time values, otherwise 50ms
%                         steps between min and max. Used with
%                         ft_topoplotER only.
%   cfg_plot.numCols    = If subplot == 1, the number of columns to plot
%   files.saveFigs      = 1 or 0. Whether to save the figures.
%
%   data                = output from ft_timelockgrandaverage
%
% See also:
%   MM_FT_CHECKCONDITIONS, MM_FT_PLOTERP

if ~isfield(cfg_ft,'parameter')
  error('Must specify cfg_ft.parameter, denoting the data to plot (e.g., ''avg'')');
end

if ~exist('sesNum','var')
  sesNum = 1;
end

cfg_plot.type = strrep(strrep(cfg_plot.ftFxn,'ft_',''),'plotER','');

% for automatically resizing figure windows
cfg_plot.screenXY = get(0,'ScreenSize');
cfg_plot.screenXY = cfg_plot.screenXY(3:4);

if ~isfield(cfg_plot,'plotTitle')
  cfg_plot.plotTitle = 0;
end
if ~isfield(cfg_plot,'plotLegend') || (strcmp(cfg_plot.type,'topo') || strcmp(cfg_plot.type,'multi'))
  cfg_plot.plotLegend = 0;
end

if ~isfield(files,'figFontName')
  files.figFontName = 'Helvetica';
end
if ~isfield(cfg_plot,'ticFontSize')
  cfg_plot.ticFontSize = [];
end
if ~isfield(cfg_plot,'labelFontSize')
  cfg_plot.labelFontSize = [];
end
if ~isfield(cfg_plot,'legendFontSize')
  cfg_plot.legendFontSize = [];
end

% check on the labels
if ~isfield(cfg_plot,'xlabel')
  cfg_plot.xlabel = 'Time (s)';
end
if ~isfield(cfg_plot,'ylabel')
  cfg_plot.ylabel = 'Voltage (\muV)';
end
cfg_plot.label_str = '';
if isfield(cfg_plot,'xlabel') && ~isempty(cfg_plot.xlabel)
  cfg_plot.label_str = cat(2,cfg_plot.label_str,'x');
end
if isfield(cfg_plot,'ylabel') && ~isempty(cfg_plot.ylabel)
  cfg_plot.label_str = cat(2,cfg_plot.label_str,'y');
end
if ~isempty(cfg_plot.label_str)
  cfg_plot.label_str = cat(2,'_',cfg_plot.label_str,'label');
end

if strcmp(cfg_plot.type,'single') || strcmp(cfg_plot.type,'multi')
  if ~isfield(cfg_ft,'fontsize')
    cfg_ft.fontsize = 9;
  end
  if ~isfield(cfg_ft,'linewidth')
    if strcmp(cfg_plot.type,'single')
      cfg_ft.linewidth = 2;
    elseif strcmp(cfg_plot.type,'multi')
      cfg_ft.linewidth = 1;
    end
  end
  if ~isfield(cfg_ft,'graphcolor')
    cfg_ft.graphcolor = 'rbkgcmyrbkgcmyrbkgcmy';
  end
  if ~isfield(cfg_ft,'linestyle')
    cfg_ft.linestyle = {'-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.'};
  end
  if ~isfield(cfg_plot,'x_bound')
    cfg_plot.x_bound = [];
  end
end
% not sure if this gets used
if ~isfield(cfg_plot,'excludeBadSub')
  cfg_plot.excludeBadSub = 1;
end

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
    % not sure fontsize does anything
    if ~isfield(cfg_ft,'fontsize')
      cfg_ft.fontsize = 10;
    end
    if ~isfield(cfg_ft,'markerfontsize')
      cfg_ft.markerfontsize = 9;
    end
  end
end

% % make sure conditions are set correctly
% if ~isfield(cfg_plot,'condMethod')
%   if ~iscell(cfg_plot.conditions) && (strcmp(cfg_plot.conditions,'all') || strcmp(cfg_plot.conditions,'all_across_types') || strcmp(cfg_plot.conditions,'all_within_types'))
%     cfg_plot.condMethod = 'single';
%   elseif iscell(cfg_plot.conditions) && ~iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
%     cfg_plot.condMethod = 'single';
%   elseif iscell(cfg_plot.conditions) && iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions{1}) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
%     cfg_plot.condMethod = 'single';
%   else
%     cfg_plot.condMethod = [];
%   end
% end
% cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);

% make sure conditions are set up for the for loop
if ~isfield(cfg_plot,'types')
  cfg_plot.types = repmat({''},size(cfg_plot.conditions));
end

if cfg_plot.plotLegend == 1
  % for the legend labels in line plots
  if ~isfield(cfg_plot,'rename_conditions')
    cfg_plot.rename_conditions = cfg_plot.conditions;
  end
end

% % temporary hack for plotting a single subject and 1 evVal, because the
% % data struct will have a field titled parameter and won't have the typical
% % data.ses.evVal.sub.data structure
% if ~isfield(data,cfg_ft.parameter)
%   cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);
% elseif isfield(data,cfg_ft.parameter) && length(cfg_plot.conditions{:}) > 1
%   error('Cannot have more than one condition if data is only one evVal');
% end

% % get the label info for this data struct
% if cfg_plot.is_ga && isfield(data.(cfg_plot.conditions{1}{1}),'label')
%   lab = data.(cfg_plot.conditions{1}{1}).label;
% elseif cfg_plot.is_ga && isfield(data.(cfg_plot.conditions{1}{1}),'labelcmb')
%   labcmb = data.(cfg_plot.conditions{1}{1}).labelcmb;
% elseif ~cfg_plot.is_ga && isfield(data.(cfg_plot.conditions{1}{1}),'label')
%   lab = data.(cfg_plot.conditions{1}{1}).sub(1).ses(1).data.label;
% elseif ~cfg_plot.is_ga && isfield(data.(cfg_plot.conditions{1}{1}),'labelcmb')
%   labcmb = data.(cfg_plot.conditions{1}{1}).sub(1).ses(1).data.labelcmb;
% end

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
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:});
  end
end

% time - get this info for the figure name
if isfield(cfg_ft,'xlim')
  if strcmp(cfg_ft.xlim,'maxmin')
    cfg_ft.xlim = [min(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{1}{1}).time) max(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{1}{1}).time)];
  end
else
  cfg_ft.xlim = [min(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{1}{1}).time) max(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{1}{1}).time)];
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
      if ~isfield(cfg_ft,'fontsize')
        cfg_ft.fontsize = 10;
      end
      if isfield(cfg_ft,'markerfontsize')
        cfg_ft = rmfield(cfg_ft,'markerfontsize');
      end
      cfg_plot.plotTitle = 0;
    end
  end
else
  cfg_plot.subplot = 0;
end

% each entry in cfg_plot.conditions is a cell containing this type of event
for typ = 1:length(cfg_plot.conditions)
  if strcmp(cfg_plot.type,'topo')
    % ft_topoplotER requires zlim, not ylim
    if ~isfield(cfg_ft,'zlim') && isfield(cfg_ft,'ylim')
      cfg_ft.zlim = cfg_ft.ylim;
      %cfg_ft = rmfield(cfg_ft,'ylim');
    end
    
    for evVal = 1:length(cfg_plot.conditions{typ})
      figure
      
      if cfg_plot.subplot
        for k = 1:length(cfg_plot.timeS)-1
          subplot(cfg_plot.numRows,cfg_plot.numCols,k);
          cfg_ft.xlim = [cfg_plot.timeS(k) cfg_plot.timeS(k+1)];
          if isfield(data,cfg_ft.parameter)
            % temporary hack for plotting a single subject and 1 evVal,
            % because the data struct will have a field titled parameter and
            % won't have the typical data.(exper.sesStr{sesNum}).evVal.sub.data structure
            feval(str2func(cfg_plot.ftFxn),cfg_ft,data);
          else
            feval(str2func(cfg_plot.ftFxn),cfg_ft,data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{evVal}));
          end
        end
        % reset the xlim
        cfg_ft.xlim = [cfg_plot.timeS(1) cfg_plot.timeS(end)];
      else
        feval(str2func(cfg_plot.ftFxn),cfg_ft,data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{evVal}));
        
        if cfg_plot.plotTitle
          title(sprintf('%s, %s, %.1f--%.1f s',strrep(cfg_plot.conditions{typ}{evVal},'_','-'),strrep(cfg_plot.chan_str,'_',' '),cfg_ft.xlim(1),cfg_ft.xlim(2)));
        end
        
        publishfig(gcf,~cfg_plot.plotTitle,cfg_plot.ticFontSize,cfg_plot.labelFontSize,files.figFontName,cfg_plot.legendFontSize);
        if exist('tightfig','file')
          tightfig(gcf);
        end
        
      end % subplot
      set(gcf,'Name',sprintf('%s, %s, %.1f--%.1f s',strrep(cfg_plot.conditions{typ}{evVal},'_',''),strrep(cfg_plot.chan_str,'_',' '),cfg_ft.xlim(1),cfg_ft.xlim(2)));
      if cfg_plot.plotLegend
        cfg_plot.legend_str = '_legend';
      else
        cfg_plot.legend_str = '';
      end
      if cfg_plot.subplot
        cfg_plot.subplot_str = '_subplot';
      else
        cfg_plot.subplot_str = '';
      end
      if cfg_plot.plotTitle
        cfg_plot.title_str = '_title';
      else
        cfg_plot.title_str = '';
      end
      
      if files.saveFigs
        if ~isempty(cfg_plot.types{typ})
          cfg_plot.figfilename = sprintf('tla_%s_ga_%s_%s_%s%d_%d%s%s%s%s',cfg_plot.type,cfg_plot.types{typ},cfg_plot.conditions{typ}{evVal},cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.legend_str,cfg_plot.subplot_str,cfg_plot.title_str,cfg_plot.label_str);
        else
          cfg_plot.figfilename = sprintf('tla_%s_ga_%s_%s%d_%d%s%s%s%s',cfg_plot.type,cfg_plot.conditions{typ}{evVal},cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.legend_str,cfg_plot.subplot_str,cfg_plot.title_str,cfg_plot.label_str);
        end
        dirs.saveDirFigsER = fullfile(dirs.saveDirFigs,['tla_',cfg_plot.type]);
        if ~exist(dirs.saveDirFigsER,'dir')
          mkdir(dirs.saveDirFigsER)
        end
        
        if strcmp(files.figPrintFormat(1:2),'-d')
          files.figPrintFormat = files.figPrintFormat(3:end);
        end
        if ~isfield(files,'figPrintRes')
          files.figPrintRes = 150;
        end
        print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsER,cfg_plot.figfilename));
      end
      
      % get the figure's current position and size
      cfg_plot.pos = get(gcf, 'Position');
      % get the height x width ratio
      hwRatio = cfg_plot.pos(3) / cfg_plot.pos(4);
      % % square figure
      % cfg_plot.figSize = [ceil(min(cfg_plot.screenXY) * 0.85) ceil(min(cfg_plot.screenXY) * 0.85)];
      % maintain figure height x width ratio
      cfg_plot.figSize = [ceil(min(cfg_plot.screenXY) * 0.85) ceil(min(cfg_plot.screenXY) * 0.85 * hwRatio)];
      % resize the figure window
      set(gcf, 'Units', 'pixels', 'Position', [ceil(cfg_plot.pos(1) * 0.6), cfg_plot.pos(2), cfg_plot.figSize(2), cfg_plot.figSize(1)]);
    end % for evVal
    
  elseif strcmp(cfg_plot.type,'single') || strcmp(cfg_plot.type,'multi')
    
    % get all the event values together for this type
    cfg = [];
    cfg.conditions = cfg_plot.conditions{typ};
    cfg.data_str = 'data';
    cfg.is_ga = cfg_plot.is_ga;
    %ana_str = mm_ft_catSubStr(cfg);
    ana_str = mm_catSubStr_multiSes(cfg,exper,sesNum);
    
    figure
    
    if strcmp(cfg_plot.type,'single')
      eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana_str));
      hold on
      
      if ischar(cfg_ft.ylim) && strcmp(cfg_ft.ylim,'maxmin')
        timesel = data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{1}).time >= cfg_ft.xlim(1) & data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{1}).time <= cfg_ft.xlim(2);
        voltmin = min(min(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{1}).(cfg_ft.parameter)(ismember(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{1}).label,cfg_ft.channel),timesel)));
        voltmax = max(max(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{1}).(cfg_ft.parameter)(ismember(data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{1}).label,cfg_ft.channel),timesel)));
      else
        voltmin = cfg_ft.ylim(1);
        voltmax = cfg_ft.ylim(2);
      end
      
      plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
      plot([0 0],[voltmin voltmax],'k--'); % vertical
      
      if cfg_plot.plotLegend
        legend(strrep(cfg_plot.rename_conditions{typ},'_','-'),'Location',cfg_plot.legendloc);
      end
      
      if ~isempty(cfg_plot.x_bound)
        % can't save EPS files with alpha so plot the bounds lines instead
        %h = fill([cfg_plot.x_bound(1),cfg_plot.x_bound(2),cfg_plot.x_bound(2),cfg_plot.x_bound(1)],[cfg_ft.ylim(1),cfg_ft.ylim(1),cfg_ft.ylim(2),cfg_ft.ylim(2)],cfg_plot.fillcolor);
        %set(h,'FaceAlpha',cfg_plot.fillalpha);
        %set(h,'EdgeColor',cfg_plot.filledge)
        
        plot([cfg_plot.x_bound(1) cfg_plot.x_bound(1)],[voltmin voltmax],'k'); % vertical
        plot([cfg_plot.x_bound(2) cfg_plot.x_bound(2)],[voltmin voltmax],'k'); % vertical
      end
      hold off
      
      xlabel(cfg_plot.xlabel);
      ylabel(cfg_plot.ylabel);
      if ~isempty(cfg_plot.types{typ})
        set(gcf,'Name',[strrep(cfg_plot.chan_str,'_',' '),' ',cfg_plot.types{typ}])
      else
        set(gcf,'Name',strrep(cfg_plot.chan_str,'_',' '))
      end
      
    elseif strcmp(cfg_plot.type,'multi')
      eval(sprintf('ft_multiplotER(cfg_ft,%s);',ana_str));
      %feval(str2func(cfg_plot.ftFxn),cfg_ft,data.(exper.sesStr{sesNum}).(cfg_plot.conditions{typ}{evVal}));
    end
    
    if cfg_plot.plotTitle
      if strcmp(cfg_plot.roi,'all')
        title(sprintf('%s %.1f--%.1f s',sprintf(repmat('%s,',1,length(cfg_plot.conditions{typ})),strrep(cfg_plot.conditions{typ}{:}),'_','-'),cfg_ft.xlim(1),cfg_ft.xlim(2)));
      else
        title(sprintf('%s %s, %.1f--%.1f s',sprintf(repmat('%s,',1,length(cfg_plot.conditions{typ})),strrep(cfg_plot.conditions{typ}{:}),'_','-'),strrep(cfg_plot.chan_str,'_',' '),cfg_ft.xlim(1),cfg_ft.xlim(2)));
      end
      %title(strrep(cfg_plot.chan_str,'_',' '));
    end
    if ~cfg_plot.subplot
      publishfig(gcf,~cfg_plot.plotTitle,cfg_plot.ticFontSize,cfg_plot.labelFontSize,files.figFontName,cfg_plot.legendFontSize);
    end
    if exist('tightfig','file')
      tightfig(gcf);
    end
    
    if cfg_plot.plotLegend
      cfg_plot.legend_str = '_legend';
    else
      cfg_plot.legend_str = '';
    end
    if cfg_plot.subplot
      cfg_plot.subplot_str = '_subplot';
    else
      cfg_plot.subplot_str = '';
    end
    if cfg_plot.plotTitle
      cfg_plot.title_str = '_title';
    else
      cfg_plot.title_str = '';
    end
    
    if files.saveFigs
      if ~isempty(cfg_plot.types{typ})
        cfg_plot.figfilename = sprintf('tla_%s_ga_%s_%s%s%d_%d%s%s%s%s',cfg_plot.type,cfg_plot.types{typ},sprintf(repmat('%s_',1,length(cfg_plot.conditions{typ})),cfg_plot.conditions{typ}{:}),cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.legend_str,cfg_plot.subplot_str,cfg_plot.title_str,cfg_plot.label_str);
      else
        cfg_plot.figfilename = sprintf('tla_%s_ga_%s%s%d_%d%s%s%s%s',cfg_plot.type,sprintf(repmat('%s_',1,length(cfg_plot.conditions{typ})),cfg_plot.conditions{typ}{:}),cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.legend_str,cfg_plot.subplot_str,cfg_plot.title_str,cfg_plot.label_str);
      end
      dirs.saveDirFigsER = fullfile(dirs.saveDirFigs,['tla_',cfg_plot.type]);
      if ~exist(dirs.saveDirFigsER,'dir')
        mkdir(dirs.saveDirFigsER)
      end
      
      if strcmp(files.figPrintFormat(1:2),'-d')
        files.figPrintFormat = files.figPrintFormat(3:end);
      end
      if ~isfield(files,'figPrintRes')
        files.figPrintRes = 150;
      end
      print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsER,cfg_plot.figfilename));
    end
    
    % get the figure's current position and size
    cfg_plot.pos = get(gcf, 'Position');
    % get the height x width ratio
    hwRatio = cfg_plot.pos(3) / cfg_plot.pos(4);
    % % square figure
    % cfg_plot.figSize = [ceil(min(cfg_plot.screenXY) * 0.85) ceil(min(cfg_plot.screenXY) * 0.85)];
    % maintain figure height x width ratio
    cfg_plot.figSize = [ceil(min(cfg_plot.screenXY) * 0.85) ceil(min(cfg_plot.screenXY) * 0.85 * hwRatio)];
    % resize the figure window
    set(gcf, 'Units', 'pixels', 'Position', [ceil(cfg_plot.pos(1) * 0.6), cfg_plot.pos(2), cfg_plot.figSize(2), cfg_plot.figSize(1)]);
  end % topo/single/multi
end

end
