function mm_ft_avgplotTFR_multiSes(cfg_ft,cfg_plot,ana,exper,files,dirs,sesNum,data)
%MM_FT_AVGPLOTTFR_MULTISES - Plots of average Time-Frequency data, support for multiple time windows
%

if ~isfield(cfg_ft,'parameter')
  cfg_ft.parameter = 'powspctrm';
end

if ~isfield(cfg_plot,'title')
  cfg_plot.title = '';
end

% check on the labels
if ~isfield(cfg_plot,'xlabel')
  cfg_plot.xlabel = 'Conditions';
end
if ~isfield(cfg_plot,'ylabel')
  cfg_plot.ylabel = 'Power';
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

if ~isfield(cfg_plot,'marker')
  cfg_plot.marker = {'o','x','o','x','s','^','s','^'};
  %cfg_plot.marker = {'o','x','s','^','o','x','d','v'};
end

if ~isfield(cfg_plot,'markeredgecolor')
  %cfg_plot.markcolor = {'w','k','w','k','w','k','w','k'};
  %cfg_plot.markcolor = {'w','w','w','w','w','w','w','w'};
  cfg_plot.markeredgecolor = {'b','c','r','m','b','c','r','m'};
end

if ~isfield(cfg_plot,'markerfacecolor')
  %cfg_plot.markcolor = {'w','k','w','k','w','k','w','k'};
  %cfg_plot.markcolor = {'w','w','w','w','w','w','w','w'};
  cfg_plot.markerfacecolor = repmat({'none'},1,length(cfg_plot.conditions));
end

if ~isfield(cfg_plot,'linecolor')
  %cfg_plot.markcolor = {'w','k','w','k','w','k','w','k'};
  %cfg_plot.markcolor = {'w','w','w','w','w','w','w','w'};
  cfg_plot.linecolor = repmat({'none'},1,length(cfg_plot.conditions));
end

if length(cfg_plot.marker) < length(cfg_plot.conditions)
  error('Not enough cfg_plot.marker (%d) specified for the number of conditions (%d)',length(cfg_plot.marker),length(cfg_plot.conditions));
end
if length(cfg_plot.markeredgecolor) < length(cfg_plot.conditions)
  error('Not enough cfg_plot.markeredgecolor (%d) specified for the number of conditions (%d)',length(cfg_plot.markeredgecolor),length(cfg_plot.conditions));
end

% if ~isfield(cfg_plot,'linespec')
%   % cfg_plot.linespec = 'k--o';
%   cfg_plot.linespec = {'bo','cx','ro','mx','bs','c^','rs','m^'};
% end
% if ~isfield(cfg_plot,'markcolor')
%   %cfg_plot.markcolor = {'w','k','w','k','w','k','w','k'};
%   %cfg_plot.markcolor = {'w','w','w','w','w','w','w','w'};
%   cfg_plot.markcolor = {'none','none','none','none','none','none','none','none','none'};
% end
% 
% if length(cfg_plot.linespec) < length(cfg_plot.conditions)
%   error('Not enough cfg_plot.linespec (%d) specified for the number of conditions (%d)',length(cfg_plot.linespec),length(cfg_plot.conditions));
% end
% if length(cfg_plot.markcolor) < length(cfg_plot.conditions)
%   error('Not enough cfg_plot.markcolor (%d) specified for the number of conditions (%d)',length(cfg_plot.markcolor),length(cfg_plot.conditions));
% end

if ~isfield(cfg_plot,'condNamesAtBottom')
  cfg_plot.condNamesAtBottom = true;
end

if ~isfield(cfg_plot,'plotLegend')
  cfg_plot.plotLegend = false;
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
if ~verLessThan('matlab', '8.4')
  cfg_plot.removeErrBarEnds = false;
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

numSub = length(exper.subjects) - sum(exper.badSub(:,sesNum));

if ~isfield(cfg_plot,'xlim')
  cfg_plot.xlim = [0 (length(cfg_plot.conditions)*size(cfg_plot.latency,1))+0.5];
  if size(cfg_plot.latency,1) > 1
    cfg_plot.xlim(2) = cfg_plot.xlim(2) + size(cfg_plot.latency,1) - 1;
  end
end

figure
% plot the lines
% eval(sprintf('plot([%s],cfg_plot.linespec,''LineWidth'',cfg_plot.linewidth);',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:})));
hold on

plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal
if cfg_plot.xlim(1) ~= 0
  plot([0 0],[cfg_plot.ylim(1) cfg_plot.ylim(2)],'k--'); % vertical
end

lat_str = '';
xIndCounter = 0;
for lat = 1:size(cfg_plot.latency,1)
  lat_str = sprintf('%s_%d_%d',lat_str,round(cfg_plot.latency(lat,1)*1000),round(cfg_plot.latency(lat,2)*1000));
  
  if lat > 1
    plot([xIndCounter+1 xIndCounter+1],[cfg_plot.ylim(1) cfg_plot.ylim(2)],'k-'); % vertical
    xIndCounter = xIndCounter + 1;
  end
  
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
          if length(cfg_plot.roi) <= 10
            cfg_plot.chan_str = sprintf(repmat('_%s',1,length(cfg_plot.roi)),cfg_plot.roi{:});
          else
            cfg_plot.chan_str = sprintf('_%dROI',length(cfg_plot.roi));
          end
        else
          % find the channel indices for averaging
          cfg_ana.chansel = ismember(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.label,cfg_plot.roi);
          if length(cfg_plot.roi) <= 10
            cfg_plot.chan_str = sprintf(repmat('_%s',1,length(cfg_ft.channel)),cfg_ft.channel{:});
          else
            cfg_plot.chan_str = sprintf('_%dROI',length(cfg_plot.roi));
          end
        end
        
        cfg_ana.timesel.(ev) = false(size(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time));
        tbeg = nearest(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time,cfg_plot.latency(lat,1));
        tend = nearest(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.time,cfg_plot.latency(lat,2));
        cfg_ana.timesel.(ev)(tbeg:tend) = true;
        
        cfg_ana.freqsel.(ev) = false(size(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.freq));
        fbeg = nearest(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.freq,cfg_plot.freqs(1));
        fend = nearest(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.freq,cfg_plot.freqs(2));
        cfg_ana.freqsel.(ev)(fbeg:fend) = true;
        
        cfg_ana.values.(ev)(goodSubInd) = mean(mean(mean(data.(exper.sesStr{sesNum}).(ev).sub(sub).data.(cfg_ft.parameter)(cfg_ana.chansel,cfg_ana.freqsel.(ev),cfg_ana.timesel.(ev)),3),2),1);
      end
      %end % ses
    end % sub
    cfg_ana.sem.(ev) = std(cfg_ana.values.(ev)) ./ sqrt(length(cfg_ana.values.(ev)));
  end % evVal
  
  % do the mean amplitude line plots
  if ~isfield(cfg_plot,'ylim')
    cfg_plot.ylim = eval(sprintf('[floor(min([%s])) ceil(max([%s]))]',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:}),sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:})));
  elseif isfield(cfg_plot,'ylim') && strcmp(cfg_plot.ylim,'minmax')
    cfg_plot.ylim = eval(sprintf('[floor(min([%s])) ceil(max([%s]))]',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:}),sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.conditions)),cfg_plot.conditions{:})));
  end
  
  h_d = nan(1,length(cfg_plot.plot_order));
  for c = 1:length(cfg_plot.plot_order)
    xIndCounter = xIndCounter + 1;
    
    % errorbars
    %h = errorbar(xIndCounter,mean(cfg_ana.values.(cfg_plot.plot_order{c}),1),cfg_ana.sem.(cfg_plot.plot_order{c}),cfg_plot.linespec{c},'LineWidth',cfg_plot.errwidth);
    h = errorbar(xIndCounter,mean(cfg_ana.values.(cfg_plot.plot_order{c}),1),cfg_ana.sem.(cfg_plot.plot_order{c}),'Color',cfg_plot.markeredgecolor{c},'LineWidth',cfg_plot.errwidth);
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
    %h_d(c) = plot(xIndCounter,mean(cfg_ana.values.(cfg_plot.plot_order{c}),1),cfg_plot.linespec{c},'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor{c});
    h_d(c) = plot(xIndCounter,mean(cfg_ana.values.(cfg_plot.plot_order{c}),1),'Marker',cfg_plot.marker{c},'MarkerSize',cfg_plot.markersize,'MarkerEdgeColor',cfg_plot.markeredgecolor{c},'MarkerFaceColor',cfg_plot.markerfacecolor{c},'LineWidth',cfg_plot.linewidth,'Color',cfg_plot.linecolor{c});
    
    if cfg_plot.condNamesAtBottom
      text(xIndCounter,cfg_plot.ylim(1)-.2, cfg_plot.rename_conditions{c}, 'Rotation', 15, 'FontSize', 18, 'HorizontalAlignment', 'center');
    end
  end
end

if ~isempty(cfg_plot.title)
  title(cfg_plot.title);
end

if cfg_plot.plotLegend
  legend(h_d,cfg_plot.legendtext,'Location',cfg_plot.legendloc);
end

hold off

set(gcf,'Name',sprintf('%s, %s ms',strrep(cfg_plot.chan_str,'_',' '),lat_str(2:end)))

% make it look good
axis([.5 (xIndCounter + .5) cfg_plot.ylim(1) cfg_plot.ylim(2)])
xlabel(cfg_plot.xlabel);
ylabel(cfg_plot.ylabel);
set(gca,'XTick',(1:xIndCounter))
if ~isempty(cfg_plot.xlabel)
  set(gca,'XTickLabel',strrep(cfg_plot.rename_conditions,'_','-'))
else
  set(gca,'XTickLabel',repmat({''},1,xIndCounter))
end
set(gca,'YTick',(cfg_plot.ylim(1):.5:cfg_plot.ylim(2)))
% axis square

if cfg_plot.axisxy
  %axis xy;
  set(gca,'YDir','reverse');
end

if ~isfield(files,'figFontName')
  files.figFontName = 'Helvetica';
end
publishfig(gcf,0,[],[],files.figFontName);
if exist('tightfig','file')
  tightfig(gcf);
end

screenXY = get(0,'ScreenSize');
screenXY = screenXY(3:4);
% get the figure's current position and size
pos = get(gcf, 'Position');
% get the height x width ratio
hwRatio = pos(3) / pos(4);
% multiplier = 0.85;
multiplier = 0.75;
% % square figure
% figSize = [ceil(min(screenXY) * multiplier) ceil(min(screenXY) * multiplier)];
% maintain figure height x width ratio
figSize = [ceil(min(screenXY) * multiplier) ceil(min(screenXY) * multiplier * hwRatio)];
% resize the figure window
set(gcf, 'Units', 'pixels', 'Position', [ceil(pos(1) * 0.6), pos(2), figSize(2), figSize(1)]);

if files.saveFigs
  cfg_plot.figfilename = sprintf('tfr_avg_ga_%s%d_%d_%s%s%s',sprintf(repmat('%s_',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:}),round(cfg_plot.freqs(1)),round(cfg_plot.freqs(2)),lat_str(2:end),cfg_plot.chan_str,cfg_plot.label_str);
  dirs.saveDirFigsLine = fullfile(dirs.saveDirFigs,'tfr_avg');
  if ~exist(dirs.saveDirFigsLine,'dir')
    mkdir(dirs.saveDirFigsLine)
  end
  
  if strcmp(files.figPrintFormat(1:2),'-d')
    files.figPrintFormat = files.figPrintFormat(3:end);
  end
  if ~isfield(files,'figPrintRes')
    files.figPrintRes = 150;
  end
  %print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsLine,cfg_plot.figfilename));
  screen2file(fullfile(dirs.saveDirFigsLine,cfg_plot.figfilename),files);
end

end
