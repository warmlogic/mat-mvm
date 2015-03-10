function sigElecs = space2_pow_sigElecs(conditions,freqs,latencies,clusTimes,nSigComparisons,ana,exper,files,dirs,ga_pow)

% files.saveFigs = false;
% files.figPrintFormat = 'png';

cfg = [];
cfg.parameter = 'powspctrm';

% what to plot
cfg.times = latencies;

cfg.freqs = freqs;

cfg.plotTitle = true;
cfg.plotLegend = true;

cfg.plotErrorBars = false;
cfg.eb_transp = true;

cfg.plotClusSig = true;
% cfg.clusAlpha = 0.1;
cfg.clusAlpha = 0.05;
cfg.clusTimes = clusTimes; % from actual cluster stats
% cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';

cfg.clusLimits = true;

cfg.linewidth = 2;
cfg.limitlinewidth = 0.5;

cfg.textFontSize = 10;

%cfg.ylim = [-0.6 0.6];
%cfg.ylim = [-0.5 0.2];
cfg.nCol = 3;

% =====================================================================
% finding significantly different electrodes across the entire scalp
% cfg.rois = {{'all'}};
cfg.rois = {{'C'}};
cfg.conditions = conditions;

% cfg.conditions = {{'word_rc_spac_p2','word_rc_mass_p2'}};
% cfg.conditions = {{'img_rc_spac_p2','img_rc_mass_p2'}};
% cfg.clusTimes = cfg_ana.latencies; % actual stats
cfg.clusTimes = [0.02:0.1:0.92; 0.1:0.1:1.0]'; % both halves
% cfg.clusTimes = [0.02:0.1:0.42; 0.1:0.1:0.5]'; % first half
% cfg.clusTimes = [0.52:0.1:0.92; 0.6:0.1:1.0]'; % second half

% cfg.freqs = ana.freq.theta;
% cfg.freqs = ana.freq.alpha_lower;
% cfg.freqs = ana.freq.alpha_upper;
% cfg.freqs = ana.freq.beta_lower;
cfg.sigElecAnyTime = true;
% =====================================================================

% whole power
cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-300_-100';
cfg.clusDirStr = '_avgT_avgF';
cfg.ylabel = 'Z-Trans Pow';
sigElecsAcrossComparisons = mm_ft_lineTFR(cfg,ana,exper,files,dirs,ga_pow);

% =====================================================================
% if size(sigElecsAcrossComparisons.pos{1},2) > 1
%   % nSigComparisons = 5;
%   
%   % significantly different in half of comparisons
%   nSigComparisons = floor(size(sigElecsAcrossComparisons.pos{1},2) / 2);
%   
%   % % significantly different in one fourth of comparisons
%   % nSigComparisons = floor(size(sigElecsAcrossComparisons.pos{1},2) / 4);
%   
%   % % significantly different in 5 of 28 of comparisons
%   % nSigComparisons = floor(size(sigElecsAcrossComparisons.pos{1},2) / 5);
% else
%   nSigComparisons = 1;
% end

sesNum = 1;
sigElecs = ga_pow.(exper.sesStr{sesNum}).(ana.eventValues{1}{1}{1}).label(sum(sigElecsAcrossComparisons.pos{1} + sigElecsAcrossComparisons.neg{1},2) >= nSigComparisons);
% sigElecs = sigElecs(~ismember(sigElecs,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'eyes_below','eyes_above','eyes_horiz','periph20'})}))));
sigElecs = sigElecs(~ismember(sigElecs,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'eyes_below','eyes_above','eyes_horiz','E49','E48','E43','E38','E32','E21','E17','E14','E1','E121','E120','E119','E113'})}))));

end
