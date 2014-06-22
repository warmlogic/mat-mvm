function mm_findPeak(cfg,ana,exper,data,cfg_plot)

% function mm_findPeak(cfg,ana,exper,data,cfg_plot)
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
% cfg.latency = [0.5 0.8];
% 
% % cfg.datadim = 'time';
% % cfg.roi = {'E91'};
% % cfg.latency = [0 1.0];
% 
% cfg.is_ga = true;
% cfg.sesNum = 1;
% 
% cfg.plotit = true;
% cfg.voltlim = [-3 3];
% % cfg.voltlim = [-1 5];
% 
% % only for datadim = 'elec'
% cfg.plottype = 'topo';
% % cfg.plottype = 'multi';
% 
% % only for datadim='peak2peak'
% % cfg.pospeak = [0.1 0.15];
% % cfg.negpeak = [0.15 0.2];

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
  cfg.conditions = cellflat(ana.eventValues{1});
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

% set the electrodes and times
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

% set up the data string
cfg.data_str = 'data';
ana_str = mm_catSubStr_multiSes2(cfg,exper,cfg.sesNum);

fprintf('Finding peak across conditions:%s\n',sprintf(repmat(' %s',1,length(cfg.conditions)),cfg.conditions{:}));

fprintf('Concatenating...');
allCond = eval(sprintf('ft_appenddata([],%s);',ana_str));
fprintf('Done.\n');

fprintf('Converting...');
ga_allCond = ft_timelockanalysis([],allCond);
fprintf('Done.\n');

if strcmp(cfg.datadim,'elec')
  % do this first to find the peak electrode/ROI
  
  fprintf('Selecting data of interest...');
  ga_data = ft_selectdata_new(cfg_ft,ga_allCond);
  fprintf('Done.\n');
  
  if ~isfield(cfg,'order')
    cfg.order = 'descend';
  end
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
      ft_topoplotER(cfg_plot,ga_data);
      
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
  
  if cfg.plotit
    fprintf('Plotting...');
    
    cfg_plot.channel = cfg_ft.channel;
      
    if ~isfield(cfg_plot,'xlim')
      if strcmp(cfg.latency,'all')
        cfg_plot.xlim = 'maxmin';
      else
        cfg_plot.xlim = cfg.latency;
      end
    end
    
    if ~isfield(cfg_plot,'ylim')
      cfg_plot.ylim = cfg.voltlim;
    end
    
    figure
    ft_singleplotER(cfg_plot,ga_allCond);
    fprintf('Done.\n');
  end
  
elseif strcmp(cfg.datadim,'peak2peak')
  
  cfg_ft.latency = cfg.pospeak;
  ga_pos = ft_selectdata_new(cfg_ft,ga_allCond);
  
  cfg_ft.latency = cfg.negpeak;
  ga_neg = ft_selectdata_new(cfg_ft,ga_allCond);
  
  pos_max = max(ga_pos.avg,[],2);
  neg_min = min(ga_neg.avg,[],2);
  
  posnegdiff = pos_max - abs(neg_min);
  
  if ~isfield(cfg,'order')
    cfg.order = 'descend';
  end
  
  [y,i] = sort(posnegdiff,1,cfg.order);
  
  fprintf('labels in ''%s'' order:\n',cfg.order);
  disp(ga_allCond.label(i)');
  fprintf('voltages in ''%s'' order:\n',cfg.order);
  disp(y');
  
end

