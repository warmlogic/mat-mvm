function mm_ft_lineTFR(cfg,ana,files,dirs,data)
%MM_FT_LINETFR make (and save) Time-Frequency plots
%
%   mm_ft_lineTFR(cfg,ana,files,dirs,data)
%
%
% See also:
%   MM_FT_CHECKCONDITIONS

% OLD
%
% Inputs:
%   cfg: parameters passed into the FT plotting function
%
%   cfg.ftFxn      = FieldTrip plotting function to use. Supported
%                         functions: ft_singleplotTFR, ft_topoplotTFR, and
%                         ft_multiplotTFR
%   cfg.conditions = Cell array containing cells of pairwise
%                         comparisons; Can be used for comparing a subset
%                         of events within a type.
%                         e.g., {{'T1a','T1c'}, {'T2a','T2c'}}, or it can
%                         be {{'all_within_types'}} or
%                         {{'all_across_types'}} to automatically create
%                         pairwise comparisons of event values. See
%                         MM_FT_CHECKCONDITIONS for more details.
%   cfg.plotTitle  = 1 or 0. Whether to plot the title.
%   cfg.subplot    = 1 or 0. Whether to make a subplot. cfg.xlim
%                         can be a range of time values, otherwise 50ms
%                         steps between min and max. ft_topoplotER only.
%   cfg.numCols    = If subplot == 1, the number of columns to plot
%   files.saveFigs      = 1 or 0. Whether to save the figures.
%
%   data                = output from ft_freqgrandaverage

cfg.type = 'line';

if ~isfield(cfg,'parameter')
  error('Must specify cfg.parameter, denoting the data to plot (e.g., ''powspctrm'' or ''cohspctrm'')');
end

if isfield(ana,'elec')
  cfg.layout = ft_prepare_layout([],ana);
else
  error('''ana'' struct must have ''elec'' field');
end

if ~isfield(cfg,'rois')
  error('Must set cfg.rois');
  % % use all channels in a topo or multi plot
  % cfg.roi = {{'all'}};
end

if ~iscell(cfg.rois)
  cfg.rois = {cfg.rois};
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plotting setup %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg,'plotTitle')
  cfg.plotTitle = false;
end

if ~isfield(cfg,'xlabel')
  cfg.xlabel = 'Time (s)';
end
if ~isfield(cfg,'ylabel')
  cfg.ylabel = 'Z-Transformed Power';
end

if ~isfield(cfg,'plotLegend')
  cfg.plotLegend = false;
end
if cfg.plotLegend
  % for the legend labels in line plots
  if ~isfield(cfg,'rename_conditions')
    cfg.rename_conditions = cfg.conditions;
  end
  if ~isfield(cfg,'legendloc')
    cfg.legendloc = 'SouthWest';
  end
end

if ~isfield(cfg,'yminmax')
  cfg.yminmax = [-0.5 0.5];
end
if ~isfield(cfg,'xminmax')
  cfg.xminmax = [min(cfg.times(:)) max(cfg.times(:))];
end

if ~isfield(cfg,'linewidth')
  cfg.linewidth = 2;
end
if ~isfield(cfg,'graphcolor')
  cfg.graphcolor = 'rbkgcmy';
end
if ~isfield(cfg,'linestyle')
  cfg.linestyle = {'-','--','-.','-','--','-.','-'};
end

if ~isfield(cfg,'plotClusSig')
  cfg.plotClusSig = false;
end
if cfg.plotClusSig
  if ~isfield(cfg,'clusDirStr')
    cfg.clusDirStr = '';
  end
  if ~isfield(cfg,'clusAlpha')
    cfg.clusAlpha = 0.05;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% make sure conditions are set correctly
if ~isfield(cfg,'condMethod')
  if ~iscell(cfg.conditions) && (strcmp(cfg.conditions,'all') || strcmp(cfg.conditions,'all_across_types') || strcmp(cfg.conditions,'all_within_types'))
    cfg.condMethod = 'single';
  elseif iscell(cfg.conditions) && ~iscell(cfg.conditions{1}) && length(cfg.conditions) == 1 && (strcmp(cfg.conditions{1},'all') || strcmp(cfg.conditions{1},'all_across_types') || strcmp(cfg.conditions{1},'all_within_types'))
    cfg.condMethod = 'single';
  elseif iscell(cfg.conditions) && iscell(cfg.conditions{1}) && length(cfg.conditions{1}) == 1 && (strcmp(cfg.conditions{1},'all') || strcmp(cfg.conditions{1},'all_across_types') || strcmp(cfg.conditions{1},'all_within_types'))
    cfg.condMethod = 'single';
  else
    cfg.condMethod = [];
  end
end

%cfg.conditions = mm_ft_checkConditions(cfg.conditions,ana,cfg.condMethod);

% temporary hack for plotting a single subject and 1 evVal, because the
% data struct will have a field titled parameter and won't have the typical
% data.evVal.sub.ses.data structure
if ~isfield(data,cfg.parameter)
  cfg.conditions = mm_ft_checkConditions(cfg.conditions,ana,cfg.condMethod);
elseif isfield(data,cfg.parameter) && length(cfg.conditions{:}) > 1
  error('Cannot have more than one condition if data is only one evVal');
end

% make sure conditions are set up for the for loop
if ~isfield(cfg,'types')
  cfg.types = repmat({''},size(cfg.conditions));
end

% for automatically resizing figure windows
cfg.screenXY = get(0,'ScreenSize');
cfg.screenXY = cfg.screenXY(3:4);

% % time (x-axis) - get this info for the figure name
% if isfield(cfg,'xlim')
%   if strcmp(cfg.xlim,'maxmin')
%     cfg.xlim = [min(data.(cfg.conditions{1}{1}).time) max(data.(cfg.conditions{1}{1}).time)];
%   end
% else
%   cfg.xlim = [min(data.(cfg.conditions{1}{1}).time) max(data.(cfg.conditions{1}{1}).time)];
% end
%
% % amplitude (y-axis)
%
%
% % freq (separate column) - get this info for the figure name
% if isfield(cfg,'ylim')
%   if strcmp(cfg.ylim,'maxmin')
%     cfg.ylim = [min(data.(cfg.conditions{1}{1}).freq) max(data.(cfg.conditions{1}{1}).freq)];
%   end
% else
%   cfg.ylim = [min(data.(cfg.conditions{1}{1}).freq) max(data.(cfg.conditions{1}{1}).freq)];
% end

nROIs = length(cfg.rois);
nFreq = size(cfg.freqs,1);
nTime = size(cfg.times,1);

for typ = 1:length(cfg.conditions)
  % get all the condition names together
  nCond = length(cfg.conditions{typ});
  cond_str = sprintf(repmat('%s_',1,nCond),cfg.conditions{typ}{:});
  if cfg.plotClusSig
    condCombos = nchoosek(1:nCond,2);
  end
  
  for f = 1:nFreq
    figure
    
    % subplot setup
    if ~isfield(cfg,'nCol');
      cfg.nCol = 2;
    end
    cfg.nRow = ceil(nROIs / cfg.nCol);
    
    % get all the ROI names together
    chan_str_all = cellflat(cfg.rois);
    chan_str_all = sprintf(repmat('%s_',1,length(chan_str_all)),chan_str_all{:});
    
    % initialize
    h = nan(nROIs,nCond);
    dataVec = nan(nCond,nROIs,nTime);
    
    for r = 1:nROIs
      cfg.roi = cfg.rois{r};
      
      subplot(cfg.nRow,cfg.nCol,r);
      
      plot(cfg.xminmax,[0 0],'k--'); % horizontal
      hold on;
      plot([0 0],cfg.yminmax,'k--'); % vertical
      
      % set the channel information
      if ismember(cfg.roi,ana.elecGroupsStr)
        % if it's in the predefined ROIs, get the channel numbers
        cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
        % set the string for the filename
        chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
      else
        % otherwise it should be the channel number(s) or 'all'
        cfg.channel = cfg.roi;
        
        % set the string for the filename
        if isfield(cfg,'cohrefchannel')
          chan_str = [cfg.cohrefchannel,'-',sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:})];
        else
          chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
        end
      end
      
      for evVal = 1:nCond
        for t = 1:nTime
          % find the indices for the chans, freqs, and times that we want
          chansel = ismember(data.(cfg.conditions{typ}{evVal}).label,cfg.channel);
          freqsel = data.(cfg.conditions{typ}{evVal}).freq >= cfg.freqs(f,1) & data.(cfg.conditions{typ}{evVal}).freq <= cfg.freqs(f,2);
          timesel = data.(cfg.conditions{typ}{evVal}).time >= cfg.times(t,1) & data.(cfg.conditions{typ}{evVal}).time <= cfg.times(t,2);
          
          dataVec(evVal,r,t) = nanmean(nanmean(nanmean(data.(cfg.conditions{typ}{evVal}).(cfg.parameter)(chansel,freqsel,timesel),3),2),1);
        end % t
        
        % TODO: don't plot at the average time, plot at the first time and
        % add on the final time point manually?
        h(r,evVal) = plot(mean(cfg.times,2),squeeze(dataVec(evVal,r,:)),[cfg.graphcolor(evVal),cfg.linestyle{evVal}],'LineWidth',cfg.linewidth);
      end % evVal
      
      xlim(cfg.xminmax);
      ylim(cfg.yminmax);
      xlabel(cfg.xlabel);
      ylabel(cfg.ylabel);
      hold off
      
      if ~isempty(cfg.types{typ})
        set(gcf,'Name',sprintf('%s, %s, %s, %.1f--%.1f Hz',cfg.types{typ},strrep(cond_str,'_',' '),strrep(chan_str_all,'_',' '),cfg.freqs(f,1),cfg.freqs(f,2)))
      else
        set(gcf,'Name',sprintf('%s, %s, %.1f--%.1f Hz',strrep(cond_str,'_',' '),strrep(chan_str_all,'_',' '),cfg.freqs(f,1),cfg.freqs(f,2)))
      end
      
      if cfg.plotTitle
        title(sprintf('%s, %.1f--%.1f Hz',strrep(chan_str,'_',' '),cfg.freqs(f,1),cfg.freqs(f,2)));
        cfg.title_str = '_title';
      else
        cfg.title_str = '';
      end
      
      if cfg.plotLegend && r == nROIs
        legend(h(r,:),strrep(cfg.rename_conditions{typ},'_',''),'Location',cfg.legendloc);
        cfg.legend_str = '_legend';
      else
        cfg.legend_str = '';
      end
    end % r
    
    % Find the significant clusters
    if cfg.plotClusSig
      foundpos = zeros(nROIs,nTime);
      foundneg = zeros(nROIs,nTime);
      for t = 1:size(cfg.clusTimes,1)
        dirs.saveDirClusStat = fullfile(dirs.saveDirProc,sprintf('tfr_stat_clus_%d_%d%s',cfg.clusTimes(t,1)*1000,cfg.clusTimes(t,2)*1000,cfg.clusDirStr));
        for evVal = 1:nCond
          vs_str = sprintf('%svs%s',cfg.conditions{typ}{condCombos(evVal,1)},cfg.conditions{typ}{condCombos(evVal,2)});
          savedFile = fullfile(dirs.saveDirClusStat,sprintf('tfr_stat_clus_%s_%.1f_%.1f_%d_%d.mat',vs_str,cfg.freqs(f,1),cfg.freqs(f,2),cfg.clusTimes(t,1)*1000,cfg.clusTimes(t,2)*1000));
          if exist(savedFile,'file')
            %fprintf('Loading %s\n',savedFile);
            load(savedFile);
          else
            warning([mfilename,':FileNotFound'],'No stat_clus file found for %s: %s. Going to next comparison.\n',vs_str,savedFile);
            continue
          end
          
          if ~isfield(stat_clus.(vs_str),'posclusters') && ~isfield(stat_clus.(vs_str),'negclusters')
            %fprintf('%s:\tNo positive or negative clusters\n',vs_str);
            continue
          end
          
          for r = 1:nROIs
            cfg.roi = cfg.rois{r};
            subplot(cfg.nRow,cfg.nCol,r);
            
            % set the channel information
            if ismember(cfg.roi,ana.elecGroupsStr)
              % if it's in the predefined ROIs, get the channel numbers
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
              % set the string for the filename
              chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
            else
              % otherwise it should be the channel number(s) or 'all'
              cfg.channel = cfg.roi;
              
              % set the string for the filename
              if isfield(cfg,'cohrefchannel')
                chan_str = [cfg.cohrefchannel,'-',sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:})];
              else
                chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
              end
            end
            
            % find the positive clusters
            if ~isempty(stat_clus.(vs_str).posclusters)
              sigpos = [];
              for iPos = 1:length(stat_clus.(vs_str).posclusters)
                sigpos(iPos) = stat_clus.(vs_str).posclusters(iPos).prob < cfg.clusAlpha;
              end
              sigpos = find(sigpos == 1);
              
              posCLM = squeeze(stat_clus.(vs_str).posclusterslabelmat);
              sigposCLM = zeros(size(posCLM,1),size(sigpos,2));
              for iPos = 1:length(sigpos)
                sigposCLM(:,iPos) = (posCLM == sigpos(iPos));
                
                if sum(ismember(cfg.channel,data.(cfg.conditions{typ}{evVal}).label(sigposCLM(:,iPos) > 0))) > 0
                  %fprintf('%s:\t***Found positive clusters***\n',vs_str);
                  foundpos(r,t) = foundpos(r,t) + 1;
                  text(mean(cfg.clusTimes(t,:),2) - 0.05,max(max([dataVec(:,r,t) dataVec(:,r,t)])) + (0.05 * foundpos(r,t)),...
                    sprintf('%s>%s',cfg.conditions{typ}{condCombos(evVal,1)},cfg.conditions{typ}{condCombos(evVal,2)}));
                end
              end % iPos
            end
            
            % find the negative clusters
            if ~isempty(stat_clus.(vs_str).negclusters)
              signeg = [];
              for iNeg = 1:length(stat_clus.(vs_str).negclusters)
                signeg(iNeg) = stat_clus.(vs_str).negclusters(iNeg).prob < cfg.clusAlpha;
              end
              signeg = find(signeg == 1);
              
              negCLM = squeeze(stat_clus.(vs_str).negclusterslabelmat);
              signegCLM = zeros(size(negCLM,1),size(signeg,2));
              for iNeg = 1:length(signeg)
                signegCLM(:,iNeg) = (negCLM == signeg(iNeg));
                
                if sum(ismember(cfg.channel,data.(cfg.conditions{typ}{evVal}).label(signegCLM(:,iNeg) > 0))) > 0
                  %fprintf('%s:\t***Found negative clusters***\n',vs_str);
                  foundneg(r,t) = foundneg(r,t) + 1;
                  text(mean(cfg.clusTimes(t,:),2) - 0.05,min(min([dataVec(:,r,t) dataVec(:,r,t)])) - (0.05 * foundneg(r,t)),...
                    sprintf('%s<%s',cfg.conditions{typ}{condCombos(evVal,1)},cfg.conditions{typ}{condCombos(evVal,2)}));
                end
              end % iNeg
            end
            hold off
          end % r
        end % evVal
      end % t
    end % if plotClusSig
    
    if ~isfield(files,'figFontName')
      files.figFontName = 'Helvetica';
    end
    publishfig(gcf,~cfg.plotTitle,[],[],files.figFontName);
    
    if files.saveFigs
      if ~isempty(cfg.types{typ})
        cfg.figfilename = sprintf('tfr_%s_ga_%s_%s%s%d_%d%s%s',cfg.type,cfg.types{typ},cond_str,chan_str_all,round(cfg.xminmax(1)*1000),round(cfg.xminmax(2)*1000),cfg.legend_str,cfg.title_str);
      else
        cfg.figfilename = sprintf('tfr_%s_ga_%s%s%d_%d%s%s',cfg.type,cond_str_all,chan_str,round(cfg.xminmax(1)*1000),round(cfg.xminmax(2)*1000),cfg.legend_str,cfg.title_str);
      end
      dirs.saveDirFigsTFR = fullfile(dirs.saveDirFigs,['tfr_',cfg.type]);
      if ~exist(dirs.saveDirFigsTFR,'dir')
        mkdir(dirs.saveDirFigsTFR)
      end
      
      if strcmp(files.figPrintFormat(1:2),'-d')
        files.figPrintFormat = files.figPrintFormat(3:end);
      end
      if ~isfield(files,'figPrintRes')
        files.figPrintRes = 150;
      end
      print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsTFR,cfg.figfilename));
    end
    
    % get the figure's current position and size
    cfg.pos = get(gcf, 'Position');
    % get the height x width ratio
    hwRatio = cfg.pos(3) / cfg.pos(4);
    % % square figure
    % cfg.figSize = [ceil(min(cfg.screenXY) * 0.85) ceil(min(cfg.screenXY) * 0.85)];
    % maintain figure height x width ratio
    cfg.figSize = [ceil(min(cfg.screenXY) * 0.85) ceil(min(cfg.screenXY) * 0.85 * hwRatio)];
    % resize the figure window
    set(gcf, 'Units', 'pixels', 'Position', [ceil(cfg.pos(1) * 0.6), cfg.pos(2), cfg.figSize(2), cfg.figSize(1)]);
    
  end % frq
end % typ

end
