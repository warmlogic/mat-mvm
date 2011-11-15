function mm_ft_contrasttopoER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla)
%MM_FT_CONTRASTTOPOER plot (and save) contast topoplots of ERP data
%   

if ~isfield(cfg_ft,'parameter')
  error('Must specify cfg_ft.parameter, denoting the data to plot (e.g., ''avg'' or ''individual'')');
end

if cfg_plot.plotColorbar
  cfg_plot.colorbar_str = '_cb';
else
  cfg_plot.colorbar_str = '';
end

% initialize for storing the contrast topoplots
cont_topo = [];

% get the label info
%lab = ga_tla.(ana.events.values{1}{1}).label;

% set the channel information
if ~isfield(cfg_plot,'roi')
  error('Must specify either ROI names or channel names in cfg_plot.roi');
elseif isfield(cfg_plot,'roi')
  if ismember(cfg_plot.roi,ana.elecGroupsStr)
    % if it's in the predefined ROIs, get the channel numbers
    cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
    % find the channel indices for averaging
    %cfg_plot.chansel = ismember(lab,cfg_ft.channel);
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:});
  else
    % otherwise it should be the channel number(s) or 'all'
    if ~iscell(cfg_plot.roi)
      cfg_plot.roi = {cfg_plot.roi};
    end
    
    % find the channels to highlight
    if ~strcmp(cfg_plot.roi,'all')
      cfg_ft.highlightchannel = cfg_plot.roi;
    end
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:});
  end
end

cfg_ft.layout = ft_prepare_layout([],ana);

for typ = 1:length(ana.events.values)
  if ~isfield(cfg_plot,'conditions') || strcmp(cfg_plot.conditions{typ},'all')
    % find all the pairwise event combinations
    cfg_plot.cond = nchoosek(ana.events.values{typ},2);
  else
    cfg_plot.cond = cfg_plot.conditions{typ};
  end
  % % get the names for plotting
  % cfg_plot.conditionNames = nchoosek(ana.events.names{typ},2);
  
  for c = 1:size(cfg_plot.cond,1)
    vs_str = sprintf('%s%s',cfg_plot.cond{c,1},sprintf(repmat('vs%s',1,size(cfg_plot.cond,2)-1),cfg_plot.cond{c,2:end}));
    
    % create contrast
    cont_topo.(vs_str) = ga_tla.(cfg_plot.cond{c,1});
    cont_topo.(vs_str).(cfg_ft.parameter) = ga_tla.(cfg_plot.cond{c,1}).(cfg_ft.parameter) - ga_tla.(cfg_plot.cond{c,2}).(cfg_ft.parameter);
    
    % make a plot
    figure
    ft_topoplotER(cfg_ft,cont_topo.(vs_str));
    set(gcf,'Name',sprintf('%s - %s, %.1f--%.1f s, %s',cfg_plot.cond{c,1},cfg_plot.cond{c,2},cfg_ft.xlim(1),cfg_ft.xlim(2),strrep(cfg_plot.chan_str,'_',' ')))
    if cfg_plot.plotColorbar
      h = colorbar;
      set(get(h,'YLabel'),'string','Voltage (\muV)');
    end
    if cfg_plot.plotTitle
      %title(sprintf('%s - %s, %.1f--%.1f s',cfg_plot.conditionNames{c,1},cfg_plot.conditionNames{c,2},cfg_ft.xlim(1),cfg_ft.xlim(2)));
      title(sprintf('%s - %s, %.1f--%.1f s',cfg_plot.cond{c,1},cfg_plot.cond{c,2},cfg_ft.xlim(1),cfg_ft.xlim(2)));
    end
    if ~isfield(files,'figFontName')
      files.figFontName = 'Helvetica';
    end
    publishfig(gcf,~cfg_plot.plotTitle,[],[],files.figFontName);
    
    if files.saveFigs
      if ~isempty(ana.events.types{typ})
        cfg_plot.figfilename = sprintf('tla_topo_ga_%s_%s_%s%d_%d%s',ana.events.types{typ},vs_str,cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.colorbar_str);
      else
        cfg_plot.figfilename = sprintf('tla_topo_ga_%s_%s%d_%d%s',vs_str,cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.colorbar_str);
      end
      dirs.saveDirFigsTopo = fullfile(dirs.saveDirFigs,'tla_topo');
      if ~exist(dirs.saveDirFigsTopo,'dir')
        mkdir(dirs.saveDirFigsTopo)
      end
      
      if strcmp(files.figPrintFormat(1:2),'-d')
        files.figPrintFormat = files.figPrintFormat(3:end);
      end
      if ~isfield(files,'figPrintRes')
        files.figPrintRes = 150;
      end
      print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsTopo,cfg_plot.figfilename));
    end
  end
end

end

