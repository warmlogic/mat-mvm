function mm_findPeak(cfg,ana,exper,data,cfg_ft)

% time
%
% electrode
%
% the only real field option for cfg_ft is cfg_ft.channel

if nargin < 6
  cfg_ft = [];
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
  fprintf('Finding peak across all conditions:%s\n',sprintf(repmat(' %s',1,length(cfg_plot.conditions)),cfg_plot.conditions{:}));
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

if ~isfield(cfg_ft,'channel')
  if ~isfield(cfg,'roi')
    cfg.roi = 'center109';
    fprintf('cfg.roi unset. defaulting to ''%s''\n',cfg.roi);
  end
  cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
end

if ~isfield(cfg,'latency') || (isfield(cfg,'latency') && isempty(cfg.latency))
  if ~strcmp(cfg.datadim,'peak2peak')
    cfg.latency = 'all';
    fprintf('cfg.latency unset. defaulting to ''%s''\n',cfg.latency);
    
    %if ~isfield(cfg_ft,'latency')
    cfg_ft.latency = cfg.latency;
    %end
  elseif strcmp(cfg.datadim,'peak2peak')
    cfg_ft.latency = 'all';
  end
end

cfg.data_str = 'data';

ana_str = mm_catSubStr_multiSes2(cfg,exper,cfg.sesNum);

% % %cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'posterior'})});
% % %cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS', 'RPS', 'LPI', 'RPI', 'LPS2', 'RPS2', 'LPI2', 'RPI2', 'PS', 'PI'})});
% cfg_ft.channel{end+1} = 'E81';

% cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'RPI2'})});
% cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'RPI3'})});

ga_allCond = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',ana_str));

if ~isfield(cfg,'plotit')
  cfg.plotit = true;
end

if strcmp(cfg.datadim,'elec') && cfg.plotit && ~isfield(cfg,'plottype')
  cfg.plottype = 'multi';
end

if strcmp(cfg.datadim,'elec')
  % do this first to find the peak electrode/ROI
  
  ga_data = ft_selectdata_new(cfg,ga_allCond);
  
  if ~isfield(cfg,'order')
    cfg.order = 'descend';
  end
  [y,i] = sort(mean(ga_data.avg,2),1,cfg.order);
  
  fprintf('labels in ''%s'' order:\n',cfg.order);
  disp(ga_data.label(i)');
  fprintf('voltages in ''%s'' order:\n',cfg.order);
  disp(y');
  
  if cfg.plotit
    if ~isfield(cfg_ft,'layout') && isfield(ana,'elec')
      cfg_ft.layout = ft_prepare_layout([],ana);
    end
    if ~isfield(cfg_ft,'xlim')
      if strcmp(cfg.latency,'all')
        cfg_ft.xlim = 'maxmin';
      else
        cfg_ft.xlim = cfg.latency;
      end
    end
    
    if strcmp(cfg.plottype,'topo')
      % topoplot
      
      cfg_ft.highlightchannel = cfg_ft.channel;
      cfg_ft.channel = 'all';
      cfg.highlightsymbol = '*';
      
      if ~isfield(cfg_ft,'zlim') && isfield(cfg,'zlim')
        cfg_ft.zlim = cfg.zlim;
      end
      if ~isfield(cfg_ft,'marker')
        if isfield(cfg,'marker')
          cfg.marker = cfg_ft.marker;
        else
          cfg_ft.marker = 'on';
          %cfg_ft.marker = 'label';
        end
      end
      if ~isfield(cfg_ft,'markerfontsize')
        cfg_ft.markerfontsize = 9;
      end
      
      figure
      ft_topoplotER(cfg_ft,ga_allCond);
      
    elseif strcmp(cfg.plottype,'multi')
      % multiplot
      if ~isfield(cfg_ft,'ylim')
        if isfield(cfg,'ylim')
          cfg_ft.ylim = cfg.ylim;
        else
          cfg_ft.ylim = [-5 5];
        end
      end
      if ~isfield(cfg_ft,'showlabels')
        cfg_ft.showlabels = 'yes';
      end
      if ~isfield(cfg_ft,'interactive')
        % cfg_ft.interactive = 'yes';
        cfg_ft.interactive = 'no';
      end
      if ~isfield(cfg_ft,'showoutline')
        cfg_ft.showoutline = 'yes';
      end
      if ~isfield(cfg_ft,'fontsize')
        cfg_ft.fontsize = 9;
      end
      
      figure
      ft_multiplotER(cfg_ft,ga_allCond);
    end
  end
  
elseif strcmp(cfg.datadim,'time')
  % do this once you have nailed down which is the peak electrode
  
  if cfg.plotit
    figure
    ft_singleplotER(cfg_ft,ga_allCond);
  end
  
elseif strcmp(cfg.datadim,'peak2peak')
  
  %   pospeak = [0.1 0.15];
  %
  %   negpeak = [0.15 0.2];
  
  
  cfg_sd = [];
  cfg_sd.latency = cfg.pospeak;
  ga_pos = ft_selectdata_new(cfg_sd,ga_allCond);
  
  
  cfg_sd = [];
  cfg_sd.latency = cfg.negpeak;
  ga_neg = ft_selectdata_new(cfg_sd,ga_allCond);
  
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


