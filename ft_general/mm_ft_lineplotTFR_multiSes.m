function mm_ft_lineplotTFR_multiSes(cfg_ft,cfg_plot,ana,exper,files,dirs,sesNum,data)
%MM_FT_LINEPLOTTFR_MULTISES Line plot of Time-Freq data
%   

% warning('Do not use this function for now.');

if ~isfield(cfg_ft,'parameter')
  cfg_ft.parameter = 'powspctrm';
end


% check on the labels
if ~isfield(cfg_plot,'xlabel')
  cfg_plot.xlabel = 'Conditions';
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

% some settings for plotting
if ~isfield(cfg_plot,'plot_order')
  cfg_plot.plot_order = cfg_plot.conditions;
end
% for the x-tick labels in the line plots
if ~isfield(cfg_plot,'rename_conditions')
  cfg_plot.rename_conditions = cfg_plot.plot_order;
end


% if   s_rc_mark = 'ro';
%   s_fo_mark = 'mo';
%   m_rc_mark = 'bx';
%   m_fo_mark = 'cx';

if ~isfield(cfg_plot,'linespec')
  % cfg_plot.linespec = 'k--o';
  cfg_plot.linespec = {'bo','cx','ro','mx'};
end
if ~isfield(cfg_plot,'markcolor')
  %cfg_plot.markcolor = {'w','k','w','k','w','k','w','k'};
  %cfg_plot.markcolor = {'w','w','w','w','w','w','w','w'};
  cfg_plot.markcolor = {'none','none','none','none','none','none','none','none','none'};
end

if ~isfield(cfg_plot,'plotLegend')
  cfg_plot.plotLegend = 0;
elseif isfield(cfg_plot,'plotLegend') && cfg_plot.plotLegend
  if ~isfield(cfg_plot,'legendtext')
    %cfg_plot.legendtext = {'Data'};
    cfg_plot.legendtext = cfg_plot.rename_conditions;
  end
  if ~isfield(cfg_plot,'legendloc')
    cfg_plot.legendloc = 'NorthEast';
  end
end

if ~isfield(cfg_plot,'axisxy')
  cfg_plot.axisxy = false;
end


% set up how the lines will look
if ~isfield(cfg_plot,'linewidth')
  cfg_plot.linewidth = 3;
end
if ~isfield(cfg_plot,'marksize')
  cfg_plot.marksize = 10;
end
if ~isfield(cfg_plot,'errwidth')
  cfg_plot.errwidth = 2;
end
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
if ~isfield(cfg_plot,'removeErrBarEnds')
  cfg_plot.removeErrBarEnds = true;
end


% if ~isfield(cfg_ft,'fontsize')
%   cfg_ft.fontsize = 9;
% end
% if ~isfield(cfg_ft,'graphcolor')
%   cfg_ft.graphcolor = 'rbkgcmyrbkgcmyrbkgcmy';
% end
% if ~isfield(cfg_ft,'linestyle')
%   cfg_ft.linestyle = {'-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.'};
% end
% if ~isfield(cfg_ft,'markerstyle')
%   cfg_ft.markerstyle = {'o','o','o','o','o','o','o','s','s','s','s','s','s','s'};
% end
% if ~isfield(cfg_plot,'excludeBadSub')
%   cfg_plot.excludeBadSub = 1;
% end


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

% cfg_ga = [];
% cfg_ga.parameter = 'avg';
% cfg_ga.latency = cfg_ft.latency;
% cfg_ga.keepindividual = 'no';
% cfg_ga.channel = cfg_ft.channel;
% 
% for ev = 1:length(cfg_plot.conditions)
%   cfg = [];
%   cfg.conditions = cfg_plot.conditions{ev};
%   cfg.data_str = 'data';
%   cfg.is_ga = cfg_plot.is_ga;
%   cfg.excludeBadSub = cfg_plot.excludeBadSub;
%   ana_str = mm_catSubStr_multiSes2(cfg,exper,sesNum);
%   
%   ga = eval(sprintf('ft_timelockgrandaverage(cfg_ga,%s);',ana_str.(cfg_plot.conditions{ev})));
%   
%   %[data_sel1,data_sel2] = eval(sprintf('ft_selectdata_new(cfg_ft,%s);',ana_str.(cfg_plot.conditions{ev})));
%   
%   figure
%   eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana_str));
%   legend_conds = repmat(cfg_plot.conditions,1,length(sesNum));
%   if length(sesNum) > 1
%     condInd = 0;
%     for ses = 1:length(sesNum)
%       %legend(strrep(cfg_plot.conditions,'_','-'),'Location',cfg_plot.legendloc);
%       for c = 1:length(cfg_plot.conditions)
%         condInd = condInd + 1;
%         legend_conds{condInd} = sprintf('%s-%s',legend_conds{condInd},exper.sesStr{sesNum(ses)});
%       end
%     end
%   end
%   legend(strrep(legend_conds,'_','-'),'Location',cfg_plot.legendloc);
%   set(gcf,'Name',sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}));
%   
%   plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
%   plot([0 0],[cfg_ft.ylim(1) cfg_ft.ylim(2)],'k--'); % vertical
%   
%   if cfg_plot.axisxy
%     %axis xy;
%     set(gca,'YDir','reverse');
%   end
% end

numSub = length(exper.subjects) - sum(exper.badSub(:,sesNum));
cfg_ana = struct;

% get times, data, SEM
for evVal = 1:length(cfg_plot.conditions)
  ev = cfg_plot.conditions{evVal};
  %cfg_ana.values.(ev) = nan(numSub,length(exper.sessions));
  cfg_ana.values.(ev) = nan(numSub,1);
  goodSubInd = 0;
  for sub = 1:length(exper.subjects)
    %for ses = 1:length(exper.sessions)
    if exper.badSub(sub,sesNum)
      fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
      continue
    else
      goodSubInd = goodSubInd + 1;
      
      % get the right channels (on an individual subject basis)
      if ismember(cfg_plot.roi,ana.elecGroupsStr)
        cfg_ana.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
        cfg_ana.chansel = ismember(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.label,cfg_ana.channel);
        cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:});
      else
        % find the channel indices for averaging
        cfg_ana.chansel = ismember(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.label,cfg_plot.roi);
        cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_ft.channel)),cfg_ft.channel{:});
      end
      
      cfg_ana.timesel.(ev) = false(size(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time));
      tbeg = nearest(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time,cfg_plot.latency(1));
      tend = nearest(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time,cfg_plot.latency(2));
      cfg_ana.timesel.(ev)(tbeg:tend) = true;
      
      
      %cfg_ana.timesel.(ev) = find(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time >= cfg_plot.latency(1) & data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time <= cfg_plot.latency(2));
      %cfg_ana.values.(ev)(goodSubInd,sesNum) = mean(mean(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.(cfg_ft.parameter)(cfg_ana.chansel,cfg_ana.timesel.(ev)),1),2);
      cfg_ana.values.(ev)(goodSubInd) = mean(mean(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.(cfg_ft.parameter)(cfg_ana.chansel,cfg_ana.timesel.(ev)),1),2);
    end
    %end % ses
  end % sub
  cfg_ana.sem.(ev) = std(cfg_ana.values.(ev))/sqrt(length(cfg_ana.values.(ev)));
end % evVal



% do the mean amplitude line plots
if ~isfield(cfg_plot,'ylim')
  cfg_plot.ylim = eval(sprintf('[floor(min([%s])) ceil(max([%s]))]',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:}),sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:})));
elseif isfield(cfg_plot,'ylim') && strcmp(cfg_plot.ylim,'minmax')
  cfg_plot.ylim = eval(sprintf('[floor(min([%s])) ceil(max([%s]))]',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:}),sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:})));
end

if ~isfield(cfg_plot,'xlim')
  cfg_plot.xlim = [0 length(cfg_plot.conditions)+0.5];
end



figure
% plot the lines
% eval(sprintf('plot([%s],cfg_plot.linespec,''LineWidth'',cfg_plot.linewidth);',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:})));
hold on

plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal
plot([0 0],[cfg_plot.ylim(1) cfg_plot.ylim(2)],'k--'); % vertical

h_d = nan(1,length(cfg_plot.plot_order));
for c = 1:length(cfg_plot.plot_order)
  % errorbars
  h = errorbar(c,mean(cfg_ana.values.(cfg_plot.plot_order{c}),1),cfg_ana.sem.(cfg_plot.plot_order{c}),cfg_plot.linespec{c},'LineWidth',cfg_plot.errwidth);
  % remove errorbar ends
  if cfg_plot.removeErrBarEnds
    chil = get(h,'Children');
    xdata = get(chil(2),'XData');
    ydata = get(chil(2),'YData');
    xdata(cfg_plot.errBarEndMarkerInd) = NaN;
    ydata(cfg_plot.errBarEndMarkerInd) = NaN;
    set(chil(2),'XData',xdata);
    set(chil(2),'YData',ydata);
    set(h,'Children',chil);
  end
  % plot the markers
  h_d(c) = plot(c,mean(cfg_ana.values.(cfg_plot.plot_order{c}),1),cfg_plot.linespec{c},'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor{c});
end
if cfg_plot.plotLegend
  legend(h_d,cfg_plot.legendtext,'Location',cfg_plot.legendloc);
end

hold off

set(gcf,'Name',sprintf('%s, %.1fs--%.1f s',strrep(cfg_plot.chan_str,'_',' '),cfg_plot.latency(1),cfg_plot.latency(2)))

% make it look good
axis([.5 (length(cfg_plot.rename_conditions) + .5) cfg_plot.ylim(1) cfg_plot.ylim(2)])
xlabel(cfg_plot.xlabel);
ylabel(cfg_plot.ylabel);
set(gca,'XTick',(1:length(cfg_plot.rename_conditions)))
if ~isempty(cfg_plot.xlabel)
  set(gca,'XTickLabel',strrep(cfg_plot.rename_conditions,'_','-'))
else
  set(gca,'XTickLabel',repmat({''},1,length(cfg_plot.rename_conditions)))
end
set(gca,'YTick',(cfg_plot.ylim(1):.5:cfg_plot.ylim(2)))
axis square

if cfg_plot.axisxy
  %axis xy;
  set(gca,'YDir','reverse');
end

if ~isfield(files,'figFontName')
  files.figFontName = 'Helvetica';
end
publishfig(gcf,0,[],[],files.figFontName);
if files.saveFigs
  cfg_plot.figfilename = sprintf('tla_line_ga_%s%s%d_%d%s',sprintf(repmat('%s_',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:}),cfg_plot.chan_str,round(cfg_plot.latency(1)*1000),round(cfg_plot.latency(2)*1000),cfg_plot.label_str);
  dirs.saveDirFigsLine = fullfile(dirs.saveDirFigs,'tla_line');
  if ~exist(dirs.saveDirFigsLine,'dir')
    mkdir(dirs.saveDirFigsLine)
  end
  
  if strcmp(files.figPrintFormat(1:2),'-d')
    files.figPrintFormat = files.figPrintFormat(3:end);
  end
  if ~isfield(files,'figPrintRes')
    files.figPrintRes = 150;
  end
  print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsLine,cfg_plot.figfilename));
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



end

