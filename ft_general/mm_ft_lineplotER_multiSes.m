function mm_ft_lineplotER_multiSes(cfg_ft,cfg_plot,ana,exper,sesNum,data)
%MM_FT_LINEPLOTER_MULTISES Line plot of ERP data
%   

warning('Do not use this function for now.');

% if ~isfield(cfg_ft,'parameter')
%   error('Must define cfg_ft.parameter for ft_singleplotER');
% end
if ~isfield(cfg_ft,'fontsize')
  cfg_ft.fontsize = 9;
end
if ~isfield(cfg_ft,'linewidth')
  cfg_ft.linewidth = 2;
end
if ~isfield(cfg_ft,'graphcolor')
  cfg_ft.graphcolor = 'rbkgcmyrbkgcmyrbkgcmy';
end
if ~isfield(cfg_ft,'linestyle')
  cfg_ft.linestyle = {'-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.'};
end
if ~isfield(cfg_ft,'markerstyle')
  cfg_ft.linestyle = {'o','o','o','o','o','o','o','s','s','s','s','s','s','s'};
end
if ~isfield(cfg_plot,'excludeBadSub')
  cfg_plot.excludeBadSub = 1;
end

if ~isfield(cfg_plot,'axisxy')
  cfg_plot.axisxy = false;
end

if ~isfield(cfg_plot,'roi')
  error('Must specify either ROI names or channel names in cfg_plot.roi');
elseif isfield(cfg_plot,'roi')
  if ismember(cfg_plot.roi,ana.elecGroupsStr)
    % if it's in the predefined ROIs, get the channel numbers
    cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
  else
    % otherwise it should be the channel number(s)
    cfg_ft.channel = cfg_plot.roi;
  end
end

% % make sure it's set up for the for loop
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

% for typ = 1:length(cfg_plot.conditions)
%   cfg = [];
%   cfg.conditions = cfg_plot.conditions{typ};
%   cfg.data_str = 'data';
%   cfg.is_ga = cfg_plot.is_ga;
%   cfg.excludeBadSub = cfg_plot.excludeBadSub;
%   ana_str = mm_ft_catSubStr(cfg,exper);
%   
%   figure
%   eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana_str));
%   legend(strrep(cfg_plot.conditions{typ},'_','-'),'Location',cfg_plot.legendloc);
%   set(gcf,'Name',sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}));
%   
%   plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
%   plot([0 0],[cfg_ft.ylim(1) cfg_ft.ylim(2)],'k--'); % vertical
% end

cfg_ga = [];
cfg_ga.parameter = 'avg';
cfg_ga.latency = cfg_ft.latency;
cfg_ga.keepindividual = 'no';
cfg_ga.channel = cfg_ft.channel;

for ev = 1:length(cfg_plot.conditions)
  cfg = [];
  cfg.conditions = cfg_plot.conditions{ev};
  cfg.data_str = 'data';
  cfg.is_ga = cfg_plot.is_ga;
  cfg.excludeBadSub = cfg_plot.excludeBadSub;
  ana_str = mm_catSubStr_multiSes2(cfg,exper,sesNum);
  
  ga = eval(sprintf('ft_timelockgrandaverage(cfg_ga,%s);',ana_str.(cfg_plot.conditions{ev})));
  
  %[data_sel1,data_sel2] = eval(sprintf('ft_selectdata_new(cfg_ft,%s);',ana_str.(cfg_plot.conditions{ev})));
  
  figure
  eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana_str));
  legend_conds = repmat(cfg_plot.conditions,1,length(sesNum));
  if length(sesNum) > 1
    condInd = 0;
    for ses = 1:length(sesNum)
      %legend(strrep(cfg_plot.conditions,'_','-'),'Location',cfg_plot.legendloc);
      for c = 1:length(cfg_plot.conditions)
        condInd = condInd + 1;
        legend_conds{condInd} = sprintf('%s-%s',legend_conds{condInd},exper.sesStr{sesNum(ses)});
      end
    end
  end
  legend(strrep(legend_conds,'_','-'),'Location',cfg_plot.legendloc);
  set(gcf,'Name',sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}));
  
  plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
  plot([0 0],[cfg_ft.ylim(1) cfg_ft.ylim(2)],'k--'); % vertical
  
  if cfg_plot.axisxy
    %axis xy;
    set(gca,'YDir','reverse');
  end
end

end

