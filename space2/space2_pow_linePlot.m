function space2_pow_linePlot(conditions,rename_conditions,freqs,latencies,sigElecs,theseColors,linestyle,ana,exper,files,dirs,ga_pow)

cfg = [];
cfg.parameter = 'powspctrm';

cfg.times = latencies;

cfg.conditions = conditions;
cfg.rename_conditions = rename_conditions;

% cfg.graphcolor = 'bcrm';
% cfg.linestyle = {'-','--','-','--'};

cfg.graphcolor = theseColors;
cfg.linestyle = linestyle;

cfg.freqs = freqs;
% cfg.freqs = ana.freq.theta;
% cfg.freqs = ana.freq.alpha_lower;
% cfg.freqs = ana.freq.alpha_upper;
% cfg.freqs = ana.freq.beta_lower;

cfg.rois = {sigElecs};

cfg.plotTitle = false;
cfg.plotLegend = true;
cfg.legendloc = 'NorthEast';

cfg.plotErrorBars = false;
cfg.eb_transp = true;

cfg.plotClusSig = false;
% cfg.clusAlpha = 0.1;
% %cfg.clusTimes = cfg.times;
% % cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
% %cfg.clusTimes = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
% cfg.clusTimes = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
% cfg.clusLimits = false;

cfg.linewidth = 2;
% cfg.limitlinewidth = 0.5;
% cfg.textFontSize = 10;

cfg.yminmax = [-1 1];
%cfg.yminmax = [-0.6 0.6];
%cfg.yminmax = [-0.5 0.2];
cfg.nCol = 3;

% % whole power
% cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
mm_ft_lineTFR(cfg,ana,exper,files,dirs,ga_pow);

end
