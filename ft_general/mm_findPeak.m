function [subjectPeaks] = mm_findPeak(cfg,ana,exper,data,cfg_plot)

% [subjectPeaks] = mm_findPeak(cfg,ana,exper,data,cfg_plot)
%
% TODO: output peak timing of subjects (if cfg.is_ga=false)
%
% cfg_plot is optional
%
% example usage:
%
% cfg = [];
% cfg.conditions = cellflat(ana.eventValues{1});
% 
% cfg.datadim = 'elec';
% cfg.roi = {'center101'};
% % cfg.roi = {'LPS','RPS'};
% cfg.latency = [0.3 0.8];
% 
% % cfg.datadim = 'time';
% % cfg.roi = {'E62'};
% % % cfg.roi = {'LPS'};
% % cfg.latency = [0 1.0];
% 
% cfg.is_ga = true;
% cfg.sesNum = 1;
% 
% cfg.plotit = true;
% cfg.voltlim = [-3 3];
% % cfg.voltlim = [-1 5];
% 
% % only for datadim='elec' and datadim='peak2peak'
% cfg.plottype = 'topo';
% % cfg.plottype = 'multi';
% 
% % % only for datadim='peak2peak'
% % cfg.datadim = 'peak2peak';
% % cfg.roi = {'posterior'};
% % cfg.pospeak = [0.08 0.14];
% % cfg.negpeak = [0.14 0.2];

%% Set defaults

if nargin < 5
  cfg_plot = [];
end

if all(~ismember(cfg.datadim,{'elec','time','peak2peak'}))
  error('cfg.datadim was set to ''%s''. options are ''elec'' and ''time'' and ''peak2peak''');
end

if ~isfield(cfg,'datadim') || isempty(cfg.datadim)
  cfg.datadim = 'elec';
  fprintf('cfg.datadim unset (options are ''elec'' and ''time''). Finding peak on dimension: ''%s''\n',cfg.datadim);
end
if strcmp(cfg.datadim,'peak2peak')
  if ~isfield(cfg,'pospeak') || ~isfield(cfg,'negpeak') 
    error('if doing ''%s'' need to set fields cfg.pospeak and cfg.negpeak',cfg.datadim);
  end
end

if ~isfield(cfg,'conditions')
  cfg.conditions = cellflat(ana.eventValues);
end

if ~isfield(cfg,'is_ga')
  cfg.is_ga = true;
  %fprintf('cfg.is_ga unset. defaulting to %d\n',cfg.is_ga);
end

if ~isfield(cfg,'excludeBadSub')
  cfg.excludeBadSub = true;
  %fprintf('cfg.excludeBadSub unset. defaulting to %d\n',cfg.excludeBadSub);
end

if ~isfield(cfg,'sesNum')
  cfg.sesNum = 1;
  fprintf('cfg.sesNum unset. defaulting to %d\n',cfg.sesNum);
end

if ~isfield(cfg,'voltlim')
  cfg.voltlim = [-5 5];
end

if ~isfield(cfg,'order')
  cfg.order = 'descend';
end

if ~isfield(cfg,'output')
  cfg.output = false;
  subjectPeaks = [];
end

if ~isfield(cfg,'plotit')
  cfg.plotit = true;
end
if cfg.plotit
  if strcmp(cfg.datadim,'elec') && cfg.plotit && ~isfield(cfg,'plottype')
    cfg.plottype = 'topo';
  end
  
  if ~isfield(cfg_plot,'interactive')
    cfg_plot.interactive = 'no';
  end
end

%% Set the electrodes and times

cfg_ft = [];
if ~isfield(cfg,'roi')
  cfg.roi = 'center101';
  fprintf('cfg.roi unset. defaulting to ''%s''\n',cfg.roi);
end
cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
if ~strcmp(cfg.datadim,'peak2peak')
  if ~isfield(cfg,'latency')
    cfg.latency = 'all';
    fprintf('cfg.latency unset. defaulting to ''%s''\n',cfg.latency);
  end
  
  cfg_ft.latency = cfg.latency;
end

%% Get the data to be processed

% set up the data string
cfg.data_str = 'data';
if cfg.is_ga
  ana_str = mm_catSubStr_multiSes2(cfg,exper,cfg.sesNum);
else
  cfg.separateSubjects = true;
  ana_str = mm_catSubStr_multiSes2(cfg,exper,cfg.sesNum);
end

cond_str = sprintf(repmat(' %s',1,length(cfg.conditions)),cfg.conditions{:});
fprintf('Finding peak across conditions:%s\n',cond_str);

if length(cfg.conditions) > 1
  if cfg.is_ga
    ga_allCond = eval(sprintf('ft_timelockgrandaverage([],%s);',ana_str));
  else
    fn = fieldnames(ana_str);
    % average each subject
    ga_allSub = cell(1,length(fn));
    for i = 1:length(fn)
      ga_allSub{i} = eval(sprintf('ft_timelockgrandaverage([],%s);',ana_str.(fn{i})));
    end
    % make the grand average
    ana_str = sprintf(repmat(',ga_allSub{%d}',1,length(ga_allSub)),1:length(ga_allSub));
    ana_str = ana_str(2:end);
    ga_allCond = eval(sprintf('ft_timelockgrandaverage([],%s);',ana_str));
  end
elseif length(cfg.conditions) == 1
  if cfg.is_ga
    ga_allCond = eval(ana_str);
  else
    fn = fieldnames(ana_str);
    % get each subject
    ga_allSub = cell(1,length(fn));
    for i = 1:length(fn)
      ga_allSub{i} = eval(ana_str.(fn{i}));
    end
    % make the grand average
    ana_str = sprintf(repmat(',ga_allSub{%d}',1,length(ga_allSub)),1:length(ga_allSub));
    ana_str = ana_str(2:end);
    ga_allCond = eval(sprintf('ft_timelockgrandaverage([],%s);',ana_str));
  end
end

%% Find the peaks

if strcmp(cfg.datadim,'elec')
  % do this first to find the peak electrode/ROI
  
  fprintf('Selecting data of interest...');
  ga_data = ft_selectdata_new(cfg_ft,ga_allCond);
  fprintf('Done.\n');
  
  [y,i] = sort(mean(ga_data.avg,2),1,cfg.order);
  
  fprintf('labels in ''%s'' order:\n',cfg.order);
  disp(ga_data.label(i)');
  fprintf('voltages in ''%s'' order:\n',cfg.order);
  disp(y');
  
  if cfg.plotit
    fprintf('Plotting...');
    
    if ~isfield(cfg_plot,'layout') && isfield(ana,'elec')
      cfg_plot.layout = ft_prepare_layout([],ana);
    end
    if ~isfield(cfg_plot,'xlim')
      if strcmp(cfg.latency,'all')
        cfg_plot.xlim = 'maxmin';
      else
        cfg_plot.xlim = cfg_ft.latency;
      end
    end
    
    if strcmp(cfg.plottype,'topo')
      % topoplot
      
      cfg_plot.highlightchannel = cfg_ft.channel;
      cfg_plot.channel = 'all';
      cfg_plot.highlight = 'labels';
      cfg_plot.colorbar = 'yes';
      
      %cfg_ft.highlight = 'on';
      %cfg_ft.highlightsize = 10;
      %cfg_ft.highlightsymbol = '*';
      
      if ~isfield(cfg_plot,'zlim')
        cfg_plot.zlim = cfg.voltlim;
      end
      if ~isfield(cfg_plot,'marker')
        if isfield(cfg,'marker')
          cfg_plot.marker = cfg.marker;
        else
          cfg_plot.marker = 'on';
          %cfg_ft.marker = 'label';
        end
      end
      if ~isfield(cfg_plot,'markerfontsize')
        cfg_plot.markerfontsize = 10;
      end
      
      figure
      ft_topoplotER(cfg_plot,ga_allCond);
      
    elseif strcmp(cfg.plottype,'multi')
      % multiplot
      if ~isfield(cfg_plot,'ylim')
        cfg_plot.ylim = cfg.voltlim;
      end
      if ~isfield(cfg_plot,'showlabels')
        cfg_plot.showlabels = 'yes';
      end
      if ~isfield(cfg_plot,'showoutline')
        cfg_plot.showoutline = 'yes';
      end
      if ~isfield(cfg_plot,'fontsize')
        cfg_plot.fontsize = 9;
      end
      
      figure
      ft_multiplotER(cfg_plot,ga_allCond);
    end
    fprintf('Done.\n');
  end
  
elseif strcmp(cfg.datadim,'time')
  % do this once you have nailed down which is the peak electrode
  
  cfg_t = [];
  cfg_t.latency = cfg.latency;
  cfg_t.channel = cfg_ft.channel;
  ga_data = ft_selectdata_new(cfg_t,ga_allCond);
  % average across all channels
  thisData = mean(ga_data.avg,1);
  
  [y,i] = sort(thisData,2,cfg.order);
  
  nMS = 20;
  if length(i) > nMS
    i = i(1:nMS);
    y = y(1:nMS);
  end
  fprintf('For channels (averaged):\n');
  disp(cfg_ft.channel);
  fprintf('In time range: %.3f sec to %.3f sec\n\n',cfg_t.latency(1),cfg_t.latency(2));
  fprintf('First %d peak samples (sec) in ''%s'' order:\n',length(i),cfg.order);
  disp(ga_data.time(i));
  fprintf('Voltages (uV) at peak milliseconds in ''%s'' order:\n',cfg.order);
  disp(y);
  
  if cfg.plotit
    fprintf('Plotting...');
    
    cfg_plot.channel = cfg_ft.channel;
      
    if ~isfield(cfg_plot,'xlim')
      if strcmp(cfg.latency,'all')
        %cfg_plot.xlim = 'maxmin';
        cfg_plot.xlim = [ga_allCond.time(1) ga_allCond.time(end)];
      else
        cfg_plot.xlim = cfg.latency;
      end
    end
    
    if ~isfield(cfg_plot,'ylim')
      cfg_plot.ylim = cfg.voltlim;
    end
    
    figure
    ft_singleplotER(cfg_plot,ga_allCond);
    hold on
    plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal
    plot([0 0],[cfg_plot.ylim(1) cfg_plot.ylim(2)],'k--'); % vertical
    hold off

    fprintf('Done.\n');
  end
  
elseif strcmp(cfg.datadim,'peak2peak')
  
  cfg_p2p = [];
  cfg_p2p.latency = cfg.pospeak;
  ga_pos = ft_selectdata_new(cfg_p2p,ga_allCond);
  
  cfg_p2p.latency = cfg.negpeak;
  ga_neg = ft_selectdata_new(cfg_p2p,ga_allCond);
  
  theseChans = ismember(ga_pos.label,cfg_ft.channel);
  pos_max = max(ga_pos.avg,[],2);
  neg_min = min(ga_neg.avg,[],2);
  
  %posnegdiff = pos_max - abs(neg_min);
  posnegdiff = pos_max - neg_min;
  
  [y,i] = sort(posnegdiff(theseChans,:),1,cfg.order);
  
  fprintf('labels in ''%s'' order:\n',cfg.order);
  disp(cfg_ft.channel(i));
  fprintf('voltages in ''%s'' order:\n',cfg.order);
  disp(y');
  
  if cfg.plotit
    fprintf('Plotting...');
    
    ga_diff = ga_pos;
    ga_diff.avg = posnegdiff;
    ga_diff.var = zeros(size(posnegdiff));
    ga_diff.time = 0;
    
    if ~isfield(cfg_plot,'layout') && isfield(ana,'elec')
      cfg_plot.layout = ft_prepare_layout([],ana);
    end
    
    if strcmp(cfg.plottype,'topo')
      % topoplot
      
      cfg_plot.xlim = 'maxmin';
      
      cfg_plot.highlightchannel = cfg_ft.channel;
      cfg_plot.channel = 'all';
      cfg_plot.highlight = 'labels';
      cfg_plot.colorbar = 'yes';
      
      %cfg_ft.highlight = 'on';
      %cfg_ft.highlightsize = 10;
      %cfg_ft.highlightsymbol = '*';
      
      if ~isfield(cfg_plot,'zlim')
        cfg_plot.zlim = cfg.voltlim;
      end
      if ~isfield(cfg_plot,'marker')
        if isfield(cfg,'marker')
          cfg_plot.marker = cfg.marker;
        else
          cfg_plot.marker = 'on';
          %cfg_ft.marker = 'label';
        end
      end
      if ~isfield(cfg_plot,'markerfontsize')
        cfg_plot.markerfontsize = 10;
      end
      
      figure
      ft_topoplotER(cfg_plot,ga_diff);
      
    elseif strcmp(cfg.plottype,'multi')
      % multiplot
      
      cfg_plot.xlim = [min([cfg.pospeak cfg.negpeak]) max([cfg.pospeak cfg.negpeak])];
      
      if ~isfield(cfg_plot,'ylim')
        cfg_plot.ylim = cfg.voltlim;
      end
      if ~isfield(cfg_plot,'showlabels')
        cfg_plot.showlabels = 'yes';
      end
      if ~isfield(cfg_plot,'showoutline')
        cfg_plot.showoutline = 'yes';
      end
      if ~isfield(cfg_plot,'fontsize')
        cfg_plot.fontsize = 9;
      end
      
      figure
      ft_multiplotER(cfg_plot,ga_allCond);
    end
    fprintf('Done.\n');
  end
  
end

if cfg.plotit
  set(gcf,'Name',cond_str);
end

