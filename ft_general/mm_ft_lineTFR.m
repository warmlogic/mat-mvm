function mm_ft_lineTFR(cfg,ana,files,dirs,data)
%MM_FT_LINETFR make (and save) Time-Frequency plots
%
%   mm_ft_lineTFR(cfg,ana,files,dirs,data)
%
% NB: figure doesn't save very well.
%
% NB: When plotting cluster permutation test significant differences, the
%     directionality labels will get printed for that ROI if even just one
%     channel in the ROI was below cfg.clusAlpha (i.e., it doesn't mean the
%     whole ROI shows a significant effect).
%
% Also, expect to have similar significance patterns across adjacent ROIs
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

if ~isfield(cfg,'type')
  cfg.type = 'line';
end
cfg.type = strrep(cfg.type,' ','');

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
    cfg.legendloc = 'NorthWest';
  end
end

if ~isfield(cfg,'yminmax')
  cfg.setylim = true;
else
  cfg.setylim = false;
end
if ~isfield(cfg,'xminmax')
  cfg.xlim = [min(reshape(cfg.times,[],1)) max(reshape(cfg.times,[],1))];
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

% cluster significance stuff
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
  if cfg.setylim
    textSpaceVert = 0.05;
    textSpaceHorz = 0.05;
  else
    textSpaceVert = 0.075;
    textSpaceHorz = 0.05;
  end
  if ~isfield(cfg,'clusSize')
    cfg.clusSize = [0.01 0.05 0.1 0.2 0.3];
  end
  if ~isfield(cfg,'clusSymb');
    cfg.clusSymb = ['*','x','+','o','.']; % [0.01 0.05 0.1 0.2 0.3]
  end
  if ~isfield(cfg,'clusLimits')
    cfg.clusLimits = false;
  end
  if ~isfield(cfg,'clusLimColor')
    cfg.clusLimColor = [0 0 1];
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

% subplot setup
if ~isfield(cfg,'nCol') &&  ~isfield(cfg,'nRow')
  cfg.nCol = 3;
  cfg.nRow = ceil(nROIs / cfg.nCol);
  fprintf('Using default cfg.nCol=%d; with %d ROIs, cfg.nRow=%d.\n',cfg.nCol,nROIs,cfg.nRow);
elseif isfield(cfg,'nCol') &&  ~isfield(cfg,'nRow')
  cfg.nRow = ceil(nROIs / cfg.nCol);
elseif ~isfield(cfg,'nCol') &&  isfield(cfg,'nRow')
  cfg.nCol = ceil(nROIs / cfg.nRow);
elseif isfield(cfg,'nCol') &&  isfield(cfg,'nRow')
  if (cfg.nCol * cfg.nRow) < nROIs
    error('Plotting %d ROIs, but you did not specify enough rows (cfg.nRow=%d) and/or columns (cfg.nCol=%d).',nROIs,cfg.nRow,cfg.nROI);
  end
end

for typ = 1:length(cfg.conditions)
  % get all the condition names together
  nCond = length(cfg.conditions{typ});
  cond_str = sprintf(repmat('%s_',1,nCond),cfg.conditions{typ}{:});
  if cfg.plotClusSig
    condCombos = nchoosek(1:nCond,2);
  end
  
  for f = 1:nFreq
    figure
    
    % get all the ROI names together
    chan_str_all = cellflat(cfg.rois);
    chan_str_all = sprintf(repmat('%s_',1,length(chan_str_all)),chan_str_all{:});
    
    % initialize
    h = nan(nROIs,nCond);
    dataVec = nan(nCond,nROIs,nTime);
    
    for r = 1:nROIs
      cfg.roi = cfg.rois{r};
      
      subplot(cfg.nRow,cfg.nCol,r);
      
      plot(cfg.xlim,[0 0],'k--'); % horizontal
      hold on
      
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
      
      if cfg.setylim
        cfg.ylim = [(min(reshape(dataVec,[],1)) - 0.1) (max(reshape(dataVec,[],1)) + 0.1)];
      end
      plot([0 0],cfg.ylim,'k--'); % vertical
      
      xlim(cfg.xlim);
      ylim(cfg.ylim);
      
      xlabel(cfg.xlabel);
      ylabel(cfg.ylabel);
      
      hold off
      
      if ~isempty(cfg.types{typ})
        type_str = sprintf('%s, ',cfg.types{typ});
      else
        type_str = '';
      end
      
      if cfg.plotClusSig
        alpha_str = sprintf(', a=%.2f',cfg.clusAlpha);
      else
        alpha_str = '';
      end
      
      set(gcf,'Name',sprintf('%s%s, %s, %.1f--%.1f Hz%s',type_str,strrep(cond_str,'_',' '),strrep(chan_str_all,'_',' '),cfg.freqs(f,1),cfg.freqs(f,2),alpha_str))
      
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
      % initialize
      foundpos = zeros(nROIs,nTime);
      foundneg = zeros(nROIs,nTime);
      
      for t = 1:size(cfg.clusTimes,1)
        dirs.saveDirClusStat = fullfile(dirs.saveDirProc,sprintf('tfr_stat_clus_%d_%d%s',round(cfg.clusTimes(t,1)*1000),round(cfg.clusTimes(t,2)*1000),cfg.clusDirStr));
        if ~exist(dirs.saveDirClusStat,'dir')
          error('%s not found!\nMake sure this directory is accessible for loading cluster permutation test results.',dirs.saveDirClusStat);
        end
        for evVal = 1:nCond
          vs_str = sprintf('%svs%s',cfg.conditions{typ}{condCombos(evVal,1)},cfg.conditions{typ}{condCombos(evVal,2)});
          savedFile = fullfile(dirs.saveDirClusStat,sprintf('tfr_stat_clus_%s_%.1f_%.1f_%d_%d.mat',vs_str,cfg.freqs(f,1),cfg.freqs(f,2),round(cfg.clusTimes(t,1)*1000),round(cfg.clusTimes(t,2)*1000)));
          if exist(savedFile,'file')
            %fprintf('Loading %s\n',savedFile);
            load(savedFile);
            %fprintf('No stat_clus file found for %s. Checking reverse condition order.\n',vs_str);
            cond1 = cfg.conditions{typ}{condCombos(evVal,1)};
            cond2 = cfg.conditions{typ}{condCombos(evVal,2)};
          else
            % try the reverse order
            vs_str = sprintf('%svs%s',cfg.conditions{typ}{condCombos(evVal,2)},cfg.conditions{typ}{condCombos(evVal,1)});
            savedFile = fullfile(dirs.saveDirClusStat,sprintf('tfr_stat_clus_%s_%.1f_%.1f_%d_%d.mat',vs_str,cfg.freqs(f,1),cfg.freqs(f,2),round(cfg.clusTimes(t,1)*1000),round(cfg.clusTimes(t,2)*1000)));
            if exist(savedFile,'file')
              %fprintf('Loading %s\n',savedFile);
              load(savedFile);
              cond1 = cfg.conditions{typ}{condCombos(evVal,2)};
              cond2 = cfg.conditions{typ}{condCombos(evVal,1)};
            else
              warning([mfilename,':FileNotFound'],'No stat_clus file found for %s (or reverse order): %s. Going to next comparison.\n',vs_str,savedFile);
              continue
            end
          end
          
          if ~isfield(stat_clus.(vs_str),'posclusters') && ~isfield(stat_clus.(vs_str),'negclusters')
            %fprintf('%s:\tNo positive or negative clusters\n',vs_str);
            continue
          end
          
          for r = 1:nROIs
            % choose the right roi/subplot
            cfg.roi = cfg.rois{r};
            subplot(cfg.nRow,cfg.nCol,r);
            
            % set the channel information
            if ismember(cfg.roi,ana.elecGroupsStr)
              % if it's in the predefined ROIs, get the channel numbers
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
            else
              % otherwise it should be the channel number(s)
              cfg.channel = cfg.roi;
            end
            
            if cfg.clusLimits
              if cfg.setylim
                cfg.ylim = [(min(reshape(dataVec,[],1)) - 0.1) (max(reshape(dataVec,[],1)) + 0.1)];
              end
              clusLimColorL = cfg.clusLimColor;
              clusLimColorR = cfg.clusLimColor;
              if t == 1
                clusLimColorL = [0 0 0];
              elseif t == size(cfg.clusTimes,1)
                clusLimColorR = [0 0 0];
              end
              hold on
              plot([cfg.clusTimes(t,1) cfg.clusTimes(t,1)],cfg.ylim,'Color',clusLimColorL); % vertical
              plot([cfg.clusTimes(t,2) cfg.clusTimes(t,2)],cfg.ylim,'Color',clusLimColorR); % vertical
              hold off
            end
            
            % find the positive clusters
            sigpos = [];
            if ~isempty(stat_clus.(vs_str).posclusters)
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
                  
                  % debug
                  fprintf('%.1fto%.1fHz: %.2fto%.2fms, %5s, %s: Pos clus #%d (%d chan), p=%.3f (found=%d)\n',...
                    cfg.freqs(f,1),cfg.freqs(f,2),cfg.clusTimes(t,1),cfg.clusTimes(t,2),cell2mat(cfg.roi),vs_str,sigpos(iPos),...
                    sum(ismember(cfg.channel,data.(cfg.conditions{typ}{evVal}).label(sigposCLM(:,iPos) > 0))),...
                    stat_clus.(vs_str).posclusters(iPos).prob,foundpos(r,t));
                  
                  clus_symb = cfg.clusSymb(find(stat_clus.(vs_str).posclusters(iPos).prob < cfg.clusSize,1,'first'));
                  clus_str = sprintf('%s>%s:%s',cond1,cond2,clus_symb);
                  text(mean(cfg.clusTimes(t,:),2) - (textSpaceHorz * 2),max(reshape(dataVec(:,r,t),[],1)) + (textSpaceVert * foundpos(r,t)),...
                    clus_str);
                end
              end % iPos
            end
            
            % find the negative clusters
            signeg = [];
            if ~isempty(stat_clus.(vs_str).negclusters)
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
                  
                  % debug
                  fprintf('%.1fto%.1fHz: %.2fto%.2fms, %5s, %s: Neg clus #%d (%d chan), p=%.3f (found=%d)\n',...
                    cfg.freqs(f,1),cfg.freqs(f,2),cfg.clusTimes(t,1),cfg.clusTimes(t,2),cell2mat(cfg.roi),vs_str,signeg(iNeg),...
                    sum(ismember(cfg.channel,data.(cfg.conditions{typ}{evVal}).label(signegCLM(:,iNeg) > 0))),...
                    stat_clus.(vs_str).negclusters(iNeg).prob,foundneg(r,t));
                  
                  clus_symb = cfg.clusSymb(find(stat_clus.(vs_str).negclusters(iNeg).prob < cfg.clusSize,1,'first'));
                  clus_str = sprintf('%s<%s:%s',cond1,cond2,clus_symb);
                  text(mean(cfg.clusTimes(t,:),2) - (textSpaceHorz * 2),min(reshape(dataVec(:,r,t),[],1)) - (textSpaceVert * foundneg(r,t)),...
                    clus_str);
                end
              end % iNeg
            end
            %hold off
          end % r
        end % evVal
      end % t
    end % if plotClusSig
    
    if ~isfield(files,'figFontName')
      files.figFontName = 'Helvetica';
    end
    
    if files.saveFigs
      if ~isempty(cfg.types{typ})
        type_str = sprintf('%s_',cfg.types{typ});
      else
        type_str = '';
      end
      
      if cfg.plotClusSig
        alpha_str = sprintf('%.2f',cfg.clusAlpha);
        alpha_str = sprintf('_a%s',alpha_str(end-1:end));
      else
        alpha_str = '';
      end
      
      cfg.figfilename = sprintf('tfr_%s_ga_%s%s%s%d_%d_%d_%d%s%s%s',cfg.type,type_str,cond_str,chan_str_all,round(cfg.xlim(1)*1000),round(cfg.xlim(2)*1000),round(cfg.freqs(f,1)),round(cfg.freqs(f,2)),alpha_str,cfg.legend_str,cfg.title_str);
      
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
    
    publishfig(gcf,~cfg.plotTitle,[],[],files.figFontName);
    
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
