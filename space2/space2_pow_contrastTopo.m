function space2_pow_contrastTopo(conditions,freqs,latencies,sigElecs,ana,exper,files,dirs,ga_pow)

% files.saveFigs = true;
% files.figPrintFormat = 'png';

cfg_plot = [];
cfg_plot.plotTitle = false;

cfg_plot.conditions = conditions;

cfg_ft = [];

cfg_ft.ylim = freqs;
% cfg_ft.ylim = ana.freq.alpha_lower;
% cfg_ft.ylim = ana.freq.alpha_upper;
% cfg_ft.ylim = ana.freq.beta_lower;

cfg_ft.parameter = 'powspctrm';

% will mask electrodes not in cfg_plot.roi
% cfg_ft.maskparameter = '';
cfg_ft.maskparameter = 'datamask';
cfg_plot.maskvalue = 0;

cfg_ft.interactive = 'no';
% cfg_ft.interactive = 'yes';

%cfg_ft.colormap = hot(64);
cfg_ft.colormap = jet(64);
cfg_ft.colorbar = 'yes';

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

%cfg_plot.ftFxn = 'ft_singleplotTFR';
cfg_plot.ftFxn = 'ft_topoplotTFR';
cfg_plot.zlabel = 'Z-Transformed Power';
%cfg_ft.marker = 'on';
cfg_ft.marker = 'off';
% cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.comment = 'no';
%cfg_ft.xlim = [0.5 0.8]; % time
cfg_plot.subplot = 0;
% cfg_ft.xlim = [0 0.5]; % time
% cfg_ft.xlim = [0.52 1.0]; % time
cfg_ft.xlim = latencies; % time

% cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time
%cfg_ft.xlim = (0:0.05:1.0); % time
%cfg_plot.roi = {'PS'};

cfg_ft.shading = 'interp';
% cfg_ft.style = 'fill';

cfg_plot.roi = sigElecs;

% cfg_ft.zlim = [-0.5 0.5]; % color
cfg_ft.zlim = [-1.0 1.0]; % color

% cfg_plot.ftFxn = 'ft_multiplotTFR';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';

ses = 1;
mm_ft_contrastTFR(cfg_ft,cfg_plot,ana,exper,files,dirs,ga_pow,ses);

end
