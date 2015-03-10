function space2_pow_avgPlot(conditions,rename_conditions,freqs,latencies,sigElecs,ana,exper,files,dirs,data_pow)

cfg_ft = [];
cfg_ft.parameter = 'powspctrm';

cfg_plot = [];

cfg_plot.conditions = conditions;
cfg_plot.plot_order = cfg_plot.conditions;
cfg_plot.rename_conditions = rename_conditions;

cfg_plot.condNamesAtBottom = false;

cfg_plot.ylabel = 'Z-Transformed Power';

cfg_plot.marker = {'o','x','s','^','o','x','d','v'};
cfg_plot.markersize = 20;

% cfg_plot.markeredgecolor = {'b','c','r','m','r','m','r','m'};
if exist('linspecer','file')
  if length(cfg_plot.conditions) == 8
    thisColor = linspecer(length(cfg_plot.conditions));
  else
    % add in OnePres
    thisColor = linspecer(8);
    thisColor = [[0 0 0]; [0 1 0]; thisColor];
    cfg_plot.marker = cat(2,'o','x',cfg_plot.marker);
  end
  cfg_plot.markeredgecolor = cell(1,size(thisColor,1));
  for i = 1:size(thisColor,1)
    cfg_plot.markeredgecolor{i} = thisColor(i,:);
  end
else
  cfg_plot.markeredgecolor = {'b','c','r','m','r','m','r','m'};
end

% cfg_plot.linecolor = repmat({'none'},1,length(cfg_plot.conditions));
% cfg_plot.markerfacecolor = repmat({'none'},1,length(cfg_plot.conditions));

% cfg_plot.freqs = ana.freq.theta;
% cfg_plot.freqs = ana.freq.alpha_lower;
% cfg_plot.freqs = ana.freq.alpha_upper;
% cfg_plot.freqs = ana.freq.beta_lower;

cfg_plot.freqs = freqs;
cfg_plot.roi = sigElecs;

cfg_plot.legendloc = 'NorthEast';

cfg_plot.ylim = [-1 1];

cfg_plot.latency = latencies;
% cfg_plot.latency = [0 0.5; 0.5 1.0];
% cfg_plot.latency = [0 0.333; 0.333 0.666; 0.666 1.0];
% cfg_plot.latency = [0 0.5];
% cfg_plot.latency = [0.5 1.0];

cfg_plot.plotLegend = true;
cfg_plot.xlabel = '';

sesNum = 1;
mm_ft_avgplotTFR_multiSes(cfg_ft,cfg_plot,ana,exper,files,dirs,sesNum,data_pow);

end
