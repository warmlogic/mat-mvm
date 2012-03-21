function [cfg_ana] = mm_ft_printVolt(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data)
%mm_ft_printVolt - print voltages to a file
% 
% See also:
%   MM_FT_CHECKCONDITIONS

if ~isfield(cfg_ft,'parameter')
  error('Must specify cfg_ft.parameter, denoting the data to test (e.g., ''avg'' or ''individual'')');
end

cfg_ft.method = 'analytic';
cfg_ft.statistic = 'depsamplesT';
cfg_ft.computestat = 'yes';
cfg_ft.computecritval = 'yes';
cfg_ft.computeprob = 'yes';
if ~isfield(cfg_ft,'tail')
  cfg_ft.tail = 0; % -1=left, 0=both, 1=right
end
if ~isfield(cfg_ft,'alpha')
  cfg_ft.alpha = 0.05;
end
if ~isfield(cfg_ft,'correctm')
  cfg_ft.correctm = 'no';
end

if ~isfield(cfg_plot,'individ_plots')
  cfg_plot.individ_plots = 0;
end
if ~isfield(cfg_plot,'line_plots')
  cfg_plot.line_plots = 0;
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

if ~isfield(cfg_ana,'excludeBadSub')
  cfg_ana.excludeBadSub = 1;
elseif isfield(cfg_ana,'excludeBadSub') && cfg_ana.excludeBadSub ~= 1
  fprintf('Must exclude bad subjects. Setting cfg_ana.excludeBadSub = 1;');
  cfg_ana.excludeBadSub = 1;
end

% get the label info for this data struct
if isfield(data.(exper.eventValues{1}).sub(1).ses(1).data,'label');
  lab = data.(exper.eventValues{1}).sub(1).ses(1).data.label;
else
  error('label information not found in data struct');
end

% set the channel information
if ~isfield(cfg_ana,'roi')
  error('Must specify either ROI names or channel names in cfg_ana.roi');
elseif isfield(cfg_ana,'roi')
  if ismember(cfg_ana.roi,ana.elecGroupsStr)
    % if it's in the predefined ROIs, get the channel numbers
    cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi)});
    % find the channel indices for averaging
    cfg_ana.chansel = ismember(lab,cfg_ft.channel);
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:});
  else
    % otherwise it should be the channel number(s) or 'all'
    if ~iscell(cfg_ana.roi)
      cfg_ana.roi = {cfg_ana.roi};
    end
    cfg_ft.channel = cfg_ana.roi;
    
    % find the channel indices for averaging
    if strcmp(cfg_ft.channel,'all')
      cfg_ana.chansel = ismember(lab,ft_channelselection(cfg_ft.channel,lab));
    else
      cfg_ana.chansel = ismember(lab,cfg_ft.channel);
    end
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_ft.channel)),cfg_ft.channel{:});
  end
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

% collect all conditions into one cell
allConds = unique(cat(2,cfg_ana.conditions{:}));

% some settings for plotting
if cfg_plot.line_plots == 1
  if ~isfield(cfg_plot,'plot_order')
    cfg_plot.plot_order = allConds;
  end
  % for the x-tick labels in the line plots
  if ~isfield(cfg_plot,'rename_conditions')
    cfg_plot.rename_conditions = cfg_plot.plot_order;
  end
end

% get times, data, SEM
for evVal = 1:length(allConds)
  ev = allConds{evVal};
  cfg_ana.values.(ev) = nan(cfg_ana.numSub,length(exper.sessions));
  goodSubInd = 0;
  for sub = 1:length(exper.subjects)
    ses = 1;
    %for ses = 1:length(exper.sessions)
    if exper.badSub(sub,ses)
      fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
      continue
    else
      goodSubInd = goodSubInd + 1;
      cfg_ana.timesel.(ev) = find(data.(ev).sub(sub).ses(ses).data.time >= cfg_ft.latency(1) & data.(ev).sub(sub).ses(ses).data.time <= cfg_ft.latency(2));
      cfg_ana.values.(ev)(goodSubInd,ses) = mean(mean(data.(ev).sub(sub).ses(ses).data.(cfg_ft.parameter)(cfg_ana.chansel,cfg_ana.timesel.(ev)),1),2);
    end
    %end % ses
  end % sub
  cfg_ana.sem.(ev) = std(cfg_ana.values.(ev))/sqrt(length(cfg_ana.values.(ev)));
  
  fprintf(cfg_ana.fid,'%s\t%s%s\n',ev,cfg_ana.roi{:},sprintf(repmat('\t%.4f',1,length(cfg_ana.values.(ev))),cfg_ana.values.(ev)));
  
end % evVal



end
