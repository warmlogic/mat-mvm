function [peakInfo] = mm_findPeak(cfg,ana,exper,data,cfg_plot)

% [peakInfo] = mm_findPeak(cfg,ana,exper,data,cfg_plot)
%
% cfg_plot is optional
%
% example usage:
%
% cfg = [];
% cfg.conditions = cellflat(ana.eventValues{1});
% 
% % Find the peak electrode (averages across the latency window)
% cfg.datadim = 'elec';
% cfg.roi = {'center101'};
% % cfg.roi = {'LPS','RPS'};
% cfg.latency = [0.4 0.8]; % LPC
% 
% % % Once you've found the peak electrode(s), find the peak time point
% % % (averages across the electrodes). After that, home in on a window for
% % % a subsequent analysis (demonstrated here).
% % cfg.datadim = 'time';
% % cfg.order = 'descend'; % descend = positive peaks first
% % cfg.roi = {'E62'};
% % % cfg.roi = {'LPS'};
% % lpcPeak = 0.592;
% % % cfg.latency = [lpcPeak-0.05 lpcPeak+0.05]; % around GA peak +/- 50
% % cfg.latency = [lpcPeak-0.1 lpcPeak+0.1]; % around GA peak +/- 100
% % % To get the average within your specified latency, set
% % % cfg.avgovertime=true (only for datadim='time'); otherwise voltage at
% % % peak time sample is output
% 
% % % If finding negative peaks (e.g., N1, N400), use order='ascend'
% % cfg.order = 'ascend'; % ascend = negative peaks first
%
% % cfg.is_ga = true;
% cfg.is_ga = false;
% cfg.outputSubjects = true;
% cfg.sesNum = 1;
% 
% cfg.plotit = true;
% cfg.voltlim = [-3 3];
% % cfg.voltlim = [-1 5];
% 
% % only for datadim='elec' and datadim='peak2peak'
% cfg.plottype = 'topo'; % (default)
% % cfg.plottype = 'multi';
% 
% % % You can also measure a peak to peak difference with pospeak/negpeak
% % cfg.datadim = 'peak2peak';
% % cfg.order = 'descend'; % descend = positive differences first
% % cfg.roi = {'posterior'};
% % cfg.pospeak = [0.08 0.14];
% % cfg.negpeak = [0.14 0.2];
%
% [peakInfo] = mm_findPeak(cfg,ana,exper,data_tla);

%% Set defaults

if nargin < 5
  cfg_plot = [];
end

if ~isfield(cfg,'datadim') || isempty(cfg.datadim)
  cfg.datadim = 'elec';
  error('cfg.datadim unset (options are ''elec'', ''time'', and ''peak2peak'').');
end

if ~ismember(cfg.datadim,{'elec','time','peak2peak'})
  error('cfg.datadim was set to ''%s''. options are ''elec'' and ''time'' and ''peak2peak''');
end

if strcmp(cfg.datadim,'time')
  if ~isfield(cfg,'avgovertime')
    cfg.avgovertime = false;
  end
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
  % cfg.is_ga = true;
  % % fprintf('cfg.is_ga unset. defaulting to %d\n',cfg.is_ga);
  error('cfg.is_ga unset. please specify whether you are using grand average data or not');
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

peakInfo = [];
peakInfo.datadim = cfg.datadim;
if ~isfield(cfg,'outputSubjects')
  cfg.outputSubjects = false;
end
if cfg.is_ga && cfg.outputSubjects
  cfg.outputSubejcts = false;
end
if cfg.outputSubjects
  if ~isfield(cfg,'nPoints_sub')
    cfg.nPoints_sub = 10;
  end
end

if ~isfield(cfg,'plotit')
  cfg.plotit = true;
end
if cfg.plotit
  if ismember(cfg.datadim,{'elec','peak2peak'}) && cfg.plotit && ~isfield(cfg,'plottype')
    cfg.plottype = 'topo';
  end
  
  if ~isfield(cfg_plot,'interactive')
    cfg_plot.interactive = 'yes';
  end
end

if ~isfield(cfg,'rmPreviousCfg')
  cfg.rmPreviousCfg = true;
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
      if cfg.rmPreviousCfg
        ga_allSub{i}.cfg = rmfield(ga_allSub{i}.cfg,'previous');
      end
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

if cfg.rmPreviousCfg
  ga_allCond.cfg = rmfield(ga_allCond.cfg,'previous');
end

%% Find the peaks

if strcmp(cfg.datadim,'elec')
  % do this first to find the peak electrode/ROI
  
  fprintf('Selecting data of interest...');
  ga_data = ft_selectdata(cfg_ft,ga_allCond);
  fprintf('Done.\n');
  
  % sort the channels by voltage after averaging across this time window
  [y,i] = sort(mean(ga_data.avg,2),1,cfg.order);
  
  peakInfo.channel = ga_data.label(i)';
  peakInfo.latency = cfg_ft.latency;
  peakInfo.voltage = y';
  
  fprintf('In averaged time range: %.3f sec to %.3f sec\n\n',cfg_ft.latency(1),cfg_ft.latency(2));
  fprintf('channels in ''%s'' order:\n',cfg.order);
  disp(peakInfo.channel);
  fprintf('voltages in ''%s'' order:\n',cfg.order);
  disp(peakInfo.voltage);
  
  if cfg.plotit
    %fprintf('Plotting...');
    if ~isfield(cfg_plot,'layout') && isfield(ana,'elec')
      cfg_plot.layout = ft_prepare_layout([],ana);
    end
    if ~isfield(cfg_plot,'xlim')
      if ischar(cfg_ft.latency) && strcmp(cfg_ft.latency,'all')
        %cfg_plot.xlim = 'maxmin';
        cfg_plot.xlim = [ga_allCond.time(1) ga_allCond.time(end)];
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
    drawnow;
    %fprintf('Done.\n');
  end
  
  % Save individual subject data
  if ~cfg.is_ga && cfg.outputSubjects
    fprintf('Getting data for individual subjects...\n');
    
    if length(cfg_ft.channel) < cfg.nPoints_sub
      cfg.nPoints_sub = length(cfg_ft.channel);
    end
    
    peakInfo.subjects.channel = cell(length(ga_allSub),cfg.nPoints_sub);
    peakInfo.subjects.voltage = nan(length(ga_allSub),cfg.nPoints_sub);
    
    for gs = 1:length(ga_allSub)
      sub_data = ft_selectdata(cfg_ft,ga_allSub{gs});
      [y,i] = sort(mean(sub_data.avg,2),1,cfg.order);
      
      channels = sub_data.label(i);
      peakInfo.subjects.channel(gs,:) = channels(1:cfg.nPoints_sub)';
      peakInfo.subjects.voltage(gs,:) = y(1:cfg.nPoints_sub)';
    end
    fprintf('Done.\n');
  end
  
elseif strcmp(cfg.datadim,'time')
  % do this once you have nailed down which is the peak electrode
  
  cfg_sd = [];
  cfg_sd.latency = cfg_ft.latency;
  cfg_sd.channel = cfg_ft.channel;
  ga_data = ft_selectdata(cfg_sd,ga_allCond);
  % sort time samples by voltage after averaging across all channels
  [y,i] = sort(mean(ga_data.avg,1),2,cfg.order);
  
  peakInfo.channel = cfg_ft.channel;
  peakInfo.latency = ga_data.time(i);
  peakInfo.voltage = y;
  
  if cfg.avgovertime
    % average over the entire latency (so reordering with sort doesn't
    % matter)
    peakInfo.voltage = mean(peakInfo.voltage);
  end
  
  % just choose a subset of timepoints and voltages
  nPoints_ga = 20;
  if length(peakInfo.latency) > nPoints_ga
    % give the sorted latency even though voltage may be averaged
    peakInfo.latency = peakInfo.latency(1:nPoints_ga);
    if ~cfg.avgovertime
      % only output voltage at the top n timepoints
      peakInfo.voltage = peakInfo.voltage(1:nPoints_ga);
    end
  end
  
  fprintf('For channels (averaged):\n');
  disp(peakInfo.channel);
  fprintf('In time range: %.3f sec to %.3f sec\n\n',cfg_sd.latency(1),cfg_sd.latency(2));
  fprintf('First %d peak samples (sec) in ''%s'' order:\n',length(peakInfo.latency),cfg.order);
  disp(peakInfo.latency);
  fprintf('Voltages (uV) at peak milliseconds in ''%s'' order:\n',cfg.order);
  disp(peakInfo.voltage);
  
  if cfg.plotit
    %fprintf('Plotting...');
    cfg_plot.channel = cfg_ft.channel;
      
    if ~isfield(cfg_plot,'xlim')
      if ischar(cfg_ft.latency) && strcmp(cfg_ft.latency,'all')
        %cfg_plot.xlim = 'maxmin';
        cfg_plot.xlim = [ga_allCond.time(1) ga_allCond.time(end)];
      else
        cfg_plot.xlim = cfg_ft.latency;
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
    
    drawnow
    %fprintf('Done.\n');
  end
  
  % Save individual subejct data
  if ~cfg.is_ga && cfg.outputSubjects
    fprintf('Getting data for individual subjects...\n');
    
    if length(peakInfo.latency) < cfg.nPoints_sub
      cfg.nPoints_sub = length(peakInfo.latency);
    end
    
    peakInfo.subjects.latency = nan(length(ga_allSub),cfg.nPoints_sub);
    if ~cfg.avgovertime
      peakInfo.subjects.voltage = nan(length(ga_allSub),cfg.nPoints_sub);
    else
      peakInfo.subjects.voltage = nan(length(ga_allSub),1);
    end
    
    for gs = 1:length(ga_allSub)
      sub_data = ft_selectdata(cfg_ft,ga_allSub{gs});
      [y,i] = sort(mean(sub_data.avg,1),2,cfg.order);
      
      % give the sorted latency even though voltage may be averaged
      sub_latency = sub_data.time(i);
      peakInfo.subjects.latency(gs,:) = sub_latency(1:cfg.nPoints_sub)';
      if ~cfg.avgovertime
        % only output voltage at the top n timepoints
        peakInfo.subjects.voltage(gs,:) = y(1:cfg.nPoints_sub)';
      else
        % average over the entire latency (so reordering with sort doesn't
        % matter)
        peakInfo.subjects.voltage(gs) = mean(y);
      end
    end
    fprintf('Done.\n');
  end
  
elseif strcmp(cfg.datadim,'peak2peak')
  
  cfg_p2p = [];
  cfg_p2p.latency = cfg.pospeak;
  ga_pos = ft_selectdata(cfg_p2p,ga_allCond);
  
  cfg_p2p.latency = cfg.negpeak;
  ga_neg = ft_selectdata(cfg_p2p,ga_allCond);
  
  
  pos_max = max(ga_pos.avg,[],2);
  neg_min = min(ga_neg.avg,[],2);
  
  %posnegdiff = pos_max - abs(neg_min);
  posnegdiff = pos_max - neg_min;
  
  [y,i] = sort(posnegdiff,1,cfg.order);
  
  sortedChannels = ga_pos.label(i);
  theseChans = ismember(sortedChannels,cfg_ft.channel);
  
  peakInfo.channel = sortedChannels(theseChans);
  peakInfo.latency = [cfg.pospeak; cfg.negpeak];
  peakInfo.voltage = y';
  
  fprintf('channels in ''%s'' order:\n',cfg.order);
  disp(peakInfo.channel);
  fprintf('voltages in ''%s'' order:\n',cfg.order);
  disp(peakInfo.voltage);
  
  if cfg.plotit
    %fprintf('Plotting...');
    
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
      
      cfg_plot.xlim = [peakInfo.latency(1,1) peakInfo.latency(2,2)];
      
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
    drawnow;
    %fprintf('Done.\n');
  end
  
  % Save individual subject data
  if ~cfg.is_ga && cfg.outputSubjects
    fprintf('Getting data for individual subjects...\n');
    
    if length(cfg_ft.channel) < cfg.nPoints_sub
      cfg.nPoints_sub = length(cfg_ft.channel);
    end
    
    peakInfo.subjects.voltage = nan(length(ga_allSub),cfg.nPoints_sub);
    
    for gs = 1:length(ga_allSub)
      cfg_p2p = [];
      cfg_p2p.latency = cfg.pospeak;
      sub_pos = ft_selectdata(cfg_p2p,ga_allSub{gs});
      
      cfg_p2p.latency = cfg.negpeak;
      sub_neg = ft_selectdata(cfg_p2p,ga_allSub{gs});
      
      theseChans = ismember(sub_pos.label,cfg_ft.channel);
      pos_max = max(sub_pos.avg,[],2);
      neg_min = min(sub_neg.avg,[],2);
      
      %posnegdiff = pos_max - abs(neg_min);
      posnegdiff = pos_max - neg_min;
      
      [y,i] = sort(posnegdiff(theseChans,:),1,cfg.order);
      
      channels = cfg_ft.channel(i);
      peakInfo.subjects.channel(gs,:) = channels(1:cfg.nPoints_sub)';
      peakInfo.subjects.voltage(gs,:) = y(1:cfg.nPoints_sub)';
    end
    fprintf('Done.\n');
  end
  
end

if cfg.plotit
  set(gcf,'Name',cond_str);
end

