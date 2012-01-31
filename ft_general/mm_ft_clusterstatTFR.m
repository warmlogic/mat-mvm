function [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data)
%MM_FT_CLUSTERSTATTFR Do FT cluster statistics for time-frequency data
%
%   [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data)
%
% Will compare the event values in cfg_ana.conditions using FieldTrip's
% non-parametric cluster permutation analysis method.
%
% cfg_ana.conditions should be a cell containing one cell for each
% comparison, which contains a row of (comma separated) strings of the
% events to compare. Instead it can be {{'all_within_types'}} or
% {{'all_across_types'}} (or {{'all'}} if there is only one type) to
% automatically create pairwise comparisons of event values. See
% MM_FT_CHECKCONDITIONS for more details.
%
% See also:
%   MM_FT_CHECKCONDITIONS

if ~isfield(cfg_ft,'parameter')
  error('Must specify the cfg_ft.parameter, denoting the parameter to test (e.g., ''powspctrm'' or ''cohspctrm'')');
end

% % set the channel information
% if ~isfield(cfg_ana,'roi')
%   error('Must specify either ROI names or channel names in cfg_ana.roi');
% elseif isfield(cfg_ana,'roi')
%   if ismember(cfg_ana.roi,ana.elecGroupsStr)
%     % if it's in the predefined ROIs, get the channel numbers
%     cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi)});
%   else
%     % otherwise it should be the channel number(s) or 'all'
%     if ~iscell(cfg_ana.roi)
%       cfg_ana.roi = {cfg_ana.roi};
%     end
%     cfg_ft.channel = cfg_ana.roi;
%   end
% end

if ~isfield(cfg_ft,'avgoverchan')
  cfg_ft.avgoverchan = 'no';
end
if ~isfield(cfg_ft,'avgovertime')
  cfg_ft.avgovertime = 'no';
end
if ~isfield(cfg_ft,'avgoverfreq')
  cfg_ft.avgoverfreq = 'no';
end

cfg_ft.elec = ana.elec;
if strcmp(cfg_ft.avgoverchan,'no')
  cfg_ft.method = 'distance';
  cfg_ft.neighbours = ft_prepare_neighbours(cfg_ft);
end

% use the Monte Carlo Method to calculate the significance probability
cfg_ft.method = 'montecarlo';
cfg_ft.correctm = 'cluster';
% test statistic that will be evaluated under the permutation distribution
cfg_ft.clusterstatistic = 'maxsum';
% alpha level of the sample-specific test statistic that will be used for
% thresholding
if ~isfield(cfg_ft,'clusteralpha')
  cfg_ft.clusteralpha = 0.05;
end
% alpha level of the permutation test
if ~isfield(cfg_ft,'alpha')
  cfg_ft.alpha = 0.025;
end
% number of draws from the permutation distribution
if ~isfield(cfg_ft,'numrandomization')
  cfg_ft.numrandomization = 1000;
end
% minimum number of neighborhood channels that is required for a selected
% sample to be included in the clustering algorithm (default = 0)
if ~isfield(cfg_ft,'minnbchan')
  cfg_ft.minnbchan = 2;
end

% exclude the bad subjects from the subject count
cfg_ana.numSub = length(exper.subjects) - sum(exper.badSub);

% make sure cfg_ana.conditions is set correctly
if ~isfield(cfg_ana,'condMethod') || isempty(cfg_ana.condMethod)
  if ~iscell(cfg_ana.conditions) && (strcmp(cfg_ana.conditions,'all') || strcmp(cfg_ana.conditions,'all_across_types') || strcmp(cfg_ana.conditions,'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  elseif iscell(cfg_ana.conditions) && ~iscell(cfg_ana.conditions{1}) && length(cfg_ana.conditions) == 1 && (strcmp(cfg_ana.conditions{1},'all') || strcmp(cfg_ana.conditions{1},'all_across_types') || strcmp(cfg_ana.conditions{1},'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  elseif iscell(cfg_ana.conditions) && iscell(cfg_ana.conditions{1}) && length(cfg_ana.conditions{1}) == 1 && (strcmp(cfg_ana.conditions{1},'all') || strcmp(cfg_ana.conditions{1},'all_across_types') || strcmp(cfg_ana.conditions{1},'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  else
    cfg_ana.condMethod = [];
  end
end
cfg_ana.conditions = mm_ft_checkConditions(cfg_ana.conditions,ana,cfg_ana.condMethod);

% set the saving directory
dirs.saveDirClusStat = fullfile(dirs.saveDirProc,sprintf('tfr_stat_clus_%d_%d',cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000));
if ~exist(dirs.saveDirClusStat,'dir')
  mkdir(dirs.saveDirClusStat)
end

for cnd = 1:length(cfg_ana.conditions)
  % initialize for storing the cluster statistics
  stat_clus = [];
  
  % initialize to note whether we found a cluster in this comparison; will
  % save stat_clus regardless of whether foundclus == 1
  foundclus = 0;
  
  % set the number of conditions that we're testing
  cfg_ana.numConds = size(cfg_ana.conditions{cnd},2);
  
  % get the strings of all the subjects in the conditions we're testing
  cfg = [];
  cfg.conditions = cfg_ana.conditions{cnd};
  cfg.data_str = 'data';
  cfg.is_ga = 0;
  ana_str = mm_ft_catSubStr(cfg,exper);
  
  % set some info to use later
  ses = 1;
  vs_str = sprintf('%s%s',cfg_ana.conditions{cnd}{1},sprintf(repmat('vs%s',1,cfg_ana.numConds-1),cfg_ana.conditions{cnd}{2:end}));
  
  % put in all the subjects for the first conditions
  subj_str = sprintf('%s',ana_str.(cfg_ana.conditions{cnd}{1}){ses});
  % do the subsequent conditions
  for i = 2:cfg_ana.numConds
    subj_str = cat(2,subj_str,sprintf(',%s',ana_str.(cfg_ana.conditions{cnd}{i}){ses}));
  end
  
  % make the design matrix
  cfg_ft.design = zeros(2,cfg_ana.numSub*cfg_ana.numConds);
  % set the unit and independent variables
  cfg_ft.uvar = 1; % the 1st row in cfg_ft.design contains the units of observation (subject number)
  cfg_ft.ivar = 2; % the 2nd row in cfg_ft.design contains the independent variable (data/condition)
  for i = 1:cfg_ana.numSub
    for j = 1:cfg_ana.numConds
      subIndex = i + ((j - 1)*cfg_ana.numSub);
      cfg_ft.design(1,subIndex) = i; % UO/DV (subject #s)
      cfg_ft.design(2,subIndex) = j; % IV (condition #s)
    end
  end
  
  % use the dependent samples T-statistic as a measure to evaluate the
  % effect at the sample level
  if cfg_ana.numConds == 2
    cfg_ft.statistic = 'depsamplesT';
    % test tails: -1 = left tail, 0 = two tail, 1 = right tail
    cfg_ft.tail = 0;
    cfg_ft.clustertail = 0;
  elseif cfg_ana.numConds > 2
    cfg_ft.statistic = 'depsamplesF';
    % test tails: -1 = left tail, 0 = two tail, 1 = right tail
    cfg_ft.tail = 1;
    cfg_ft.clustertail = 1;
  end
  
  % run the nonparametric cluster statistics
  stat_clus.(vs_str) = eval(sprintf('ft_freqstatistics(cfg_ft,%s);',subj_str));
  
  % save it if there were positive or negative clusters
  if isfield(stat_clus.(vs_str),'posclusters')
    if ~isempty(stat_clus.(vs_str).posclusters)
      foundclus = 1;
    end
  end
  if isfield(stat_clus.(vs_str),'negclusters')
    if ~isempty(stat_clus.(vs_str).negclusters)
      foundclus = 1;
    end
  end
  
  % note whether there positive or negative clusters
  stat_clus.(vs_str).foundclus = foundclus;
  
  saveFile = fullfile(dirs.saveDirClusStat,sprintf('tfr_stat_clus_%s_%.1f_%.1f_%d_%d.mat',vs_str,cfg_ft.frequency(1),cfg_ft.frequency(2),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000));
  fprintf('Saving %s\n',saveFile);
  save(saveFile,'stat_clus');
end

end
