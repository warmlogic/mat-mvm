function mm_ft_plotERP(cfg_ft,cfg_plot,ana,exper,files,dirs,data)
%MM_FT_PLOTERP make (and save) ERP plots
%
% WARNING: This is an old function! Use mm_ft_plotER instead.
%
% See also: MM_FT_PLOTER

warning([mfilename,':oldFxn'],'%s is an old function! Use mm_ft_plotER instead.',mfilename);

if ~isfield(cfg_ft,'parameter')
  error('Must define cfg_ft.parameter (e.g., ''avg'')');
end
if ~isfield(cfg_ft,'fontsize')
  cfg_ft.fontsize = 9;
end
if ~isfield(cfg_ft,'linewidth')
  cfg_ft.linewidth = 2;
end
if ~isfield(cfg_ft,'graphcolor')
  cfg_ft.graphcolor = 'rbkgcmy';
end
if ~isfield(cfg_ft,'linestyle')
  cfg_ft.linestyle = {'-','--','-.','-','--','-.','-'};
end
if ~isfield(cfg_plot,'plotTitle')
  cfg_plot.plotTitle = 0;
end
if ~isfield(cfg_plot,'plotLegend')
  cfg_plot.plotLegend = 0;
end
if ~isfield(cfg_plot,'excludeBadSub')
  cfg_plot.excludeBadSub = 1;
end

%cfg_plot.fillcolor = [.8,.8,.8];
%cfg_plot.fillalpha = 0.3;
%cfg_plot.filledge = 'none';

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
cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);
% make sure conditions are set up for the for loop
if ~isfield(cfg_plot,'types')
  cfg_plot.types = repmat({''},size(cfg_plot.conditions));
end

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
    cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
    % find the channel indices for averaging
%     cfg_plot.chansel = ismember(lab,cfg_ft.channel);
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:});
  else
    % otherwise it should be the channel number(s) or 'all'
    if ~iscell(cfg_plot.roi)
      cfg_plot.roi = {cfg_plot.roi};
    end
    cfg_ft.channel = cfg_plot.roi;
    
%     % find the channel indices for averaging
%     if strcmp(cfg_ft.channel,'all')
%       cfg_plot.chansel = ismember(lab,ft_channelselection(cfg_ft.channel,lab));
%     else
%       cfg_plot.chansel = ismember(lab,cfg_ft.channel);
%     end
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_ft.channel)),cfg_ft.channel{:});
  end
end

% each entry in cfg_plot.conditions is a cell containing this type of event
for typ = 1:length(cfg_plot.conditions)
  % get all the event values together for this type
  cfg = [];
  cfg.conditions = cfg_plot.conditions{typ};
  cfg.data_str = 'data';
  cfg.is_ga = cfg_plot.is_ga;
  cfg.excludeBadSub = cfg_plot.excludeBadSub;
  ana_str = mm_ft_catSubStr(cfg,exper);
  
  figure
  eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana_str));
  hold on
  
  plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
  plot([0 0],[cfg_ft.ylim(1) cfg_ft.ylim(2)],'k--'); % vertical
  
  if cfg_plot.plotLegend
    legend(strrep(cfg_plot.conditions{typ},'_','-'),'Location',cfg_plot.legendloc);
  end
  
  % can't save EPS files with alpha so plot the bounds lines instead
  %h = fill([cfg_plot.x_bound(1),cfg_plot.x_bound(2),cfg_plot.x_bound(2),cfg_plot.x_bound(1)],[cfg_ft.ylim(1),cfg_ft.ylim(1),cfg_ft.ylim(2),cfg_ft.ylim(2)],cfg_plot.fillcolor);
  %set(h,'FaceAlpha',cfg_plot.fillalpha);
  %set(h,'EdgeColor',cfg_plot.filledge)
  
  plot([cfg_plot.x_bound(1) cfg_plot.x_bound(1)],[cfg_ft.ylim(1) cfg_ft.ylim(2)],'k'); % vertical
  plot([cfg_plot.x_bound(2) cfg_plot.x_bound(2)],[cfg_ft.ylim(1) cfg_ft.ylim(2)],'k'); % vertical
  hold off
  
  xlabel('Time (s)');
  ylabel('Voltage (\muV)');
  if ~isempty(cfg_plot.types{typ})
    set(gcf,'Name',[strrep(cfg_plot.chan_str,'_',' '),' ',cfg_plot.types{typ}])
  else
    set(gcf,'Name',strrep(cfg_plot.chan_str,'_',' '))
  end
  if cfg_plot.plotTitle
    title(strrep(cfg_plot.chan_str,'_',' '));
  end
  %title(strrep(cfg_plot.chan_str,'_',' '))
  %axis ij % negative up
  if ~isfield(files,'figFontName')
    files.figFontName = 'Helvetica';
  end
  publishfig(gcf,~cfg_plot.plotTitle,[],[],files.figFontName);
  
  if cfg_plot.plotLegend
    cfg_plot.legend_str = '_legend';
  else
    cfg_plot.legend_str = '';
  end
  if files.saveFigs
    if ~isempty(cfg_plot.types{typ})
      cfg_plot.figfilename = sprintf('tla_erp_ga_%s_%s%s%d_%d%s',cfg_plot.types{typ},sprintf(repmat('%s_',1,length(cfg.conditions)),cfg.conditions{:}),cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.legend_str);
    else
      cfg_plot.figfilename = sprintf('tla_erp_ga_%s%s%d_%d%s',sprintf(repmat('%s_',1,length(cfg.conditions)),cfg.conditions{:}),cfg_plot.chan_str,round(cfg_ft.xlim(1)*1000),round(cfg_ft.xlim(2)*1000),cfg_plot.legend_str);
    end
    dirs.saveDirFigsERP = fullfile(dirs.saveDirFigs,'tla_erp');
    if ~exist(dirs.saveDirFigsERP,'dir')
      mkdir(dirs.saveDirFigsERP)
    end
    
    if strcmp(files.figPrintFormat(1:2),'-d')
      files.figPrintFormat = files.figPrintFormat(3:end);
    end
    if ~isfield(files,'figPrintRes')
      files.figPrintRes = 150;
    end
    print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsERP,cfg_plot.figfilename));
  end
end

end
