function sigElecsAcrossComparisons = mm_ft_lineTFR(cfg,ana,exper,files,dirs,data,sesNum)
%MM_FT_LINETFR make (and save) Time-Frequency plots
%
%   mm_ft_lineTFR(cfg,ana,exper,files,dirs,data)
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
% NB: for the output sigElecsAcrossComparisons, this will be "true" for
%     electrodes that are significantly different across a given comparison
%     using the entire scalp, not just the electrodes defined in cfg.rois.
%     You can use this to select subsets of electrodes across a bunch of
%     comparisons under the assumption that they are relatively important
%     if they are significantly different in a given number of comparisons
%     (see space2_pow.m for an example).
%
% See also:
%   MM_FT_CHECKCONDITIONS

if ~exist('sesNum','var')
  sesNum = 1;
end

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
  % cfg.rois = {{'all'}};
end

if ~iscell(cfg.rois)
  cfg.rois = {cfg.rois};
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plotting setup %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg,'plotErrorBars')
  % requires running ft_freqgrandaverage with cfg.keepindividual='yes'
  cfg.plotErrorBars = false;
end
if cfg.plotErrorBars
  if ~isfield(cfg,'eb_transp')
    cfg.eb_transp = true;
  end
  if ndims(data.(exper.sesStr{1}).(cfg.conditions{1}{1}).(cfg.parameter)) ~= 4
    error('To plot error bars, need to run ft_freqgrandaverage with cfg.keepindividual=''yes'' so data.%s is subj X elec X freq X time.',cfg.parameter)
  end
  cfg.eb_str = '_eb';
else
  cfg.eb_str = '';
end

if ~isfield(cfg,'plotTitle')
  cfg.plotTitle = false;
end

if ~isfield(cfg,'textFontSize')
  cfg.textFontSize = 10;
end

if ~isfield(cfg,'ticFontSize')
  cfg.ticFontSize = 20;
end

if ~isfield(cfg,'labelFontSize')
  cfg.labelFontSize = 24;
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
if ~isfield(cfg,'limitlinewidth')
  cfg.limitlinewidth = 0.5;
end
if ~isfield(cfg,'graphcolor')
  %cfg.graphcolor = 'rbkgcmy';
  cfg.graphcolor = 'rbkgcmyrbkgcmyrbkgcmy';
end
if ~isfield(cfg,'linestyle')
  %cfg.linestyle = {'-','--','-.','-','--','-.','-'};
  cfg.linestyle = {'-','-','-','-','-','-','-','--','--','--','--','--','--','--','-.','-.','-.','-.','-.','-.','-.'};
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

% % make sure conditions are set correctly
% if ~isfield(cfg,'condMethod')
%   if ~iscell(cfg.conditions) && (strcmp(cfg.conditions,'all') || strcmp(cfg.conditions,'all_across_types') || strcmp(cfg.conditions,'all_within_types'))
%     cfg.condMethod = 'single';
%   elseif iscell(cfg.conditions) && ~iscell(cfg.conditions{1}) && length(cfg.conditions) == 1 && (strcmp(cfg.conditions{1},'all') || strcmp(cfg.conditions{1},'all_across_types') || strcmp(cfg.conditions{1},'all_within_types'))
%     cfg.condMethod = 'single';
%   elseif iscell(cfg.conditions) && iscell(cfg.conditions{1}) && length(cfg.conditions{1}) == 1 && (strcmp(cfg.conditions{1},'all') || strcmp(cfg.conditions{1},'all_across_types') || strcmp(cfg.conditions{1},'all_within_types'))
%     cfg.condMethod = 'single';
%   else
%     cfg.condMethod = [];
%   end
% end
%
% %cfg.conditions = mm_ft_checkConditions(cfg.conditions,ana,cfg.condMethod);
%
% % temporary hack for plotting a single subject and 1 evVal, because the
% % data struct will have a field titled parameter and won't have the typical
% % data.evVal.sub.ses.data structure
% if ~isfield(data,cfg.parameter)
%   cfg.conditions = mm_ft_checkConditions(cfg.conditions,ana,cfg.condMethod);
% elseif isfield(data,cfg.parameter) && length(cfg.conditions{:}) > 1
%   error('Cannot have more than one condition if data is only one evVal');
% end

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

if nROIs == 1
  cfg.nCol = 1;
  cfg.nRow = 1;
end

% record all of the significant electrodes common across pairwise
% comparisons in the time ranges included in cfg.clusTimes
if cfg.plotClusSig
  sigElecsAcrossComparisons.pos = cell(length(cfg.conditions),nFreq);
  sigElecsAcrossComparisons.neg = cell(length(cfg.conditions),nFreq);
  if ~isfield(cfg,'sigElecAnyTime')
    % false produces only true values when an electrode was significantly
    % different between comparisons at all time windows; set to true to get
    % a true value for when any window is different
    cfg.sigElecAnyTime = false;
  end
  validComparisons = cell(1,length(cfg.conditions));
else
  sigElecsAcrossComparisons.pos = [];
  sigElecsAcrossComparisons.neg = [];
end

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
  cond_str = sprintf(repmat('%s_',1,length(cfg.conditions{typ})),cfg.conditions{typ}{:});
  if length(cfg.conditions{typ}) > 2
    nCond = length(cfg.conditions{typ});
  else
    nCond = 1;
  end
  
  if cfg.plotClusSig
    if length(cfg.conditions{typ}) >= 2
      condCombos = nchoosek(1:length(cfg.conditions{typ}),2);
    else
      error('Need two or more conditions to compare if plotting cluster significance tests');
    end
    
  end
  
  for f = 1:nFreq
    if cfg.plotClusSig
      validComparisons{typ} = true(1,size(condCombos,1));
      if cfg.sigElecAnyTime
        sigElecsAcrossComparisons.pos{typ,f} = false(length(data.(exper.sesStr{1}).(cfg.conditions{typ}{1}).label),size(condCombos,1));
        sigElecsAcrossComparisons.neg{typ,f} = false(length(data.(exper.sesStr{1}).(cfg.conditions{typ}{1}).label),size(condCombos,1));
      else
        sigElecsAcrossComparisons.pos{typ,f} = true(length(data.(exper.sesStr{1}).(cfg.conditions{typ}{1}).label),size(condCombos,1));
        sigElecsAcrossComparisons.neg{typ,f} = true(length(data.(exper.sesStr{1}).(cfg.conditions{typ}{1}).label),size(condCombos,1));
      end
    end
    
    figure
    
    % get all the ROI names together
    if ischar(cfg.rois) || (iscell(cfg.rois) && length(cellflat(cfg.rois)) <= 10)
      chan_str_all = cellflat(cfg.rois);
      chan_str_all = sprintf(repmat('_%s',1,length(chan_str_all)),chan_str_all{:});
    elseif iscell(cfg.rois) && length(cellflat(cfg.rois)) > 10
      chan_str_all = sprintf('_%dROIs',length(cellflat(cfg.rois)));
    else
      keyboard
    end
    
    % initialize
    h = nan(nROIs,nCond);
    dataVec = nan(nCond,nROIs,nTime);
    if cfg.plotErrorBars
      dataVar = nan(nCond,nROIs,nTime);
    end
    
    for r = 1:nROIs
      cfg.this_roi = cfg.rois{r};
      if ~iscell(cfg.this_roi)
        cfg.this_roi = {cfg.this_roi};
      end
      
      subplot(cfg.nRow,cfg.nCol,r);
      
      plot(cfg.xlim,[0 0],'k--'); % horizontal
      hold on
      
      % set the channel information
      if ismember(cfg.this_roi,ana.elecGroupsStr)
        % if it's in the predefined ROIs, get the channel numbers
        cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.this_roi)});
        % set the string for the filename
        chan_str = sprintf(repmat('%s_',1,length(cfg.this_roi)),cfg.this_roi{:});
      else
        % otherwise it should be the channel number(s) or 'all'
        cfg.channel = cfg.this_roi;
        
        % set the string for the filename
        if isfield(cfg,'cohrefchannel')
          chan_str = [cfg.cohrefchannel,'-',sprintf(repmat('%s_',1,length(cfg.this_roi)),cfg.this_roi{:})];
        else
          chan_str = sprintf(repmat('%s_',1,length(cfg.this_roi)),cfg.this_roi{:});
        end
      end
      
      for evVal = 1:length(cfg.conditions{typ})
        for t = 1:nTime
          % find the indices for the chans, freqs, and times that we want
          chansel = ismember(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).label,cfg.channel);
          
          fbeg = nearest(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).freq,cfg.freqs(f,1));
          fend = nearest(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).freq,cfg.freqs(f,2));
          freqsel = false(size(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).freq));
          freqsel(fbeg:fend) = true;
          
          tbeg = nearest(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).time,cfg.times(t,1));
          tend = nearest(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).time,cfg.times(t,2));
          timesel = false(size(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).time));
          timesel(tbeg:tend) = true;
          
          if ndims(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).(cfg.parameter)) == 4
            dataVec(evVal,r,t) = nanmean(nanmean(nanmean(nanmean(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).(cfg.parameter)(:,chansel,freqsel,timesel),4),3),2),1);
          else
            dataVec(evVal,r,t) = nanmean(nanmean(nanmean(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).(cfg.parameter)(chansel,freqsel,timesel),3),2),1);
          end
          if cfg.plotErrorBars
            dataVar(evVal,r,t) = nanste(nanmean(nanmean(nanmean(data.(exper.sesStr{sesNum}).(cfg.conditions{typ}{evVal}).(cfg.parameter)(:,chansel,freqsel,timesel),4),3),2));
          end
        end % t
        
        % set up the line color
        if ischar(cfg.graphcolor)
          % using matlab's single-character colors
          thisColor = cfg.graphcolor(evVal);
        elseif ~ischar(cfg.graphcolor) && iscell(cfg.graphcolor)
          % defined own color strings (for rgb.m)
          if ischar(cfg.graphcolor{evVal})
            thisColor = rgb(cfg.graphcolor{evVal});
          else
            error('Do not know what to do with cfg.graphcolor settings.');
          end
        elseif ~ischar(cfg.graphcolor) && ismatrix(cfg.graphcolor) && length(cfg.graphcolor(evVal,:)) == 3
          % defined own RGB triplets (e.g., using linspecer.m)
          thisColor = cfg.graphcolor(evVal,:);
        end
        % plot it
        if cfg.plotErrorBars
          this_h = shadedErrorBar(mean(cfg.times,2),squeeze(dataVec(evVal,r,:)),squeeze(dataVar(evVal,r,:)),{cfg.linestyle{evVal},'Color',thisColor,'LineWidth',cfg.linewidth},cfg.eb_transp);
          h(r,evVal) = this_h.mainLine;
        else
          % TODO: don't plot at the average time, plot at the first time and
          % add on the final time point manually?
          h(r,evVal) = plot(mean(cfg.times,2),squeeze(dataVec(evVal,r,:)),cfg.linestyle{evVal},'Color',thisColor,'LineWidth',cfg.linewidth);
        end
        
      end % evVal
      
      if cfg.setylim
        cfg.ylim = [(min(reshape(dataVec,[],1)) - 0.1) (max(reshape(dataVec,[],1)) + 0.1)];
      else
        cfg.ylim = [cfg.yminmax(1) cfg.yminmax(2)];
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
        legend(h(r,:),strrep(cfg.rename_conditions{typ},'_','-'),'Location',cfg.legendloc);
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
        
        %for evVal = 1:nCond
        for evVal = 1:size(condCombos,1)
          vs_str = sprintf('%svs%s',cfg.conditions{typ}{condCombos(evVal,1)},cfg.conditions{typ}{condCombos(evVal,2)});
          savedFile = fullfile(dirs.saveDirClusStat,sprintf('tfr_stat_clus_%s_%.1f_%.1f_%d_%d.mat',vs_str,cfg.freqs(f,1),cfg.freqs(f,2),round(cfg.clusTimes(t,1)*1000),round(cfg.clusTimes(t,2)*1000)));
          if exist(savedFile,'file')
            %fprintf('Loading %s\n',savedFile);
            load(savedFile);
            %fprintf('No stat_clus file found for %s. Checking reverse condition order.\n',vs_str);
            cond1 = cfg.conditions{typ}{condCombos(evVal,1)};
            cond2 = cfg.conditions{typ}{condCombos(evVal,2)};
            cond1ind = condCombos(evVal,1);
            cond2ind = condCombos(evVal,2);
          else
            % try the reverse order
            vs_str = sprintf('%svs%s',cfg.conditions{typ}{condCombos(evVal,2)},cfg.conditions{typ}{condCombos(evVal,1)});
            savedFile = fullfile(dirs.saveDirClusStat,sprintf('tfr_stat_clus_%s_%.1f_%.1f_%d_%d.mat',vs_str,cfg.freqs(f,1),cfg.freqs(f,2),round(cfg.clusTimes(t,1)*1000),round(cfg.clusTimes(t,2)*1000)));
            if exist(savedFile,'file')
              %fprintf('Loading %s\n',savedFile);
              load(savedFile);
              cond1 = cfg.conditions{typ}{condCombos(evVal,2)};
              cond2 = cfg.conditions{typ}{condCombos(evVal,1)};
              cond1ind = condCombos(evVal,2);
              cond2ind = condCombos(evVal,1);
            else
              validComparisons{typ}(evVal) = false;
              warning([mfilename,':FileNotFound'],'No stat_clus file found for %s (or reverse order): %s. Going to next comparison.\n',vs_str,savedFile);
              continue
            end
          end
          
          if ~isfield(stat_clus.(vs_str),'posclusters') && ~isfield(stat_clus.(vs_str),'negclusters')
            fprintf('%.1fto%.1fHz: %.2fto%.2f sec, %s:\tNo positive or negative clusters\n',...
              cfg.freqs(f,1),cfg.freqs(f,2),cfg.clusTimes(t,1),cfg.clusTimes(t,2),vs_str);
            %fprintf('%s:\tNo positive or negative clusters\n',vs_str);
            continue
          else
            % find the positive clusters
            if ~isempty(stat_clus.(vs_str).posclusters)
              sigpos = zeros(1,length(stat_clus.(vs_str).posclusters));
              for iPos = 1:length(stat_clus.(vs_str).posclusters)
                sigpos(iPos) = stat_clus.(vs_str).posclusters(iPos).prob <= cfg.clusAlpha;
              end
              sigpos = find(sigpos == 1);
            else
              sigpos = [];
            end
            if ~isempty(sigpos)
              posCLM = squeeze(stat_clus.(vs_str).posclusterslabelmat);
              sigposCLM = zeros(size(posCLM,1),size(sigpos,2));
              for iPos = 1:length(sigpos)
                sigposCLM(:,iPos) = (posCLM == sigpos(iPos));
              end
            end
            
            % find the negative clusters
            if ~isempty(stat_clus.(vs_str).negclusters)
              signeg = zeros(1,length(stat_clus.(vs_str).negclusters));
              for iNeg = 1:length(stat_clus.(vs_str).negclusters)
                signeg(iNeg) = stat_clus.(vs_str).negclusters(iNeg).prob <= cfg.clusAlpha;
              end
              signeg = find(signeg == 1);
            else
              signeg = [];
            end
            if ~isempty(signeg)
              negCLM = squeeze(stat_clus.(vs_str).negclusterslabelmat);
              signegCLM = zeros(size(negCLM,1),size(signeg,2));
              for iNeg = 1:length(signeg)
                signegCLM(:,iNeg) = (negCLM == signeg(iNeg));
              end
            end
            
            if cfg.sigElecAnyTime
              %sigElecsAcrossComparisons{typ,f}(stat_clus.(vs_str).prob <= cfg.clusAlpha,evVal) = true;
              if ~isempty(sigpos)
                sigElecsAcrossComparisons.pos{typ,f}(logical(sum(sigposCLM,2)),evVal) = true;
              end
              if ~isempty(signeg)
                sigElecsAcrossComparisons.neg{typ,f}(logical(sum(signegCLM,2)),evVal) = true;
              end
            else
              %sigElecsAcrossComparisons{typ,f}(stat_clus.(vs_str).prob > cfg.clusAlpha,evVal) = false;
              if ~isempty(sigpos)
                sigElecsAcrossComparisons.pos{typ,f}(~logical(sum(sigposCLM,2)),evVal) = false;
              end
              if ~isempty(signeg)
                sigElecsAcrossComparisons.neg{typ,f}(~logical(sum(signegCLM,2)),evVal) = false;
              end
            end
          end
          
          for r = 1:nROIs
            % choose the right roi/subplot
            cfg.this_roi = cfg.rois{r};
            subplot(cfg.nRow,cfg.nCol,r);
            
            if ~iscell(cfg.this_roi)
              cfg.this_roi = {cfg.this_roi};
            end
            
            % set the channel information
            if ismember(cfg.this_roi,ana.elecGroupsStr)
              % if it's in the predefined ROIs, get the channel numbers
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.this_roi)});
            else
              % otherwise it should be the channel number(s)
              cfg.channel = cfg.this_roi;
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
              plot([cfg.clusTimes(t,1) cfg.clusTimes(t,1)],cfg.ylim,'Color',clusLimColorL,'LineWidth',cfg.limitlinewidth); % vertical
              plot([cfg.clusTimes(t,2) cfg.clusTimes(t,2)],cfg.ylim,'Color',clusLimColorR,'LineWidth',cfg.limitlinewidth); % vertical
              hold off
            end
            
            % find the positive clusters
            if ~isempty(sigpos)
              for iPos = 1:length(sigpos)
                
                if sum(ismember(cfg.channel,data.(exper.sesStr{sesNum}).(cond1).label(sigposCLM(:,iPos) > 0))) > 0
                  %fprintf('%s:\t***Found positive clusters***\n',vs_str);
                  foundpos(r,t) = foundpos(r,t) + 1;
                  
                  % debug
                  fprintf('%.1fto%.1fHz: %.2fto%.2f sec, %5s, %s: Pos clus #%d (%d chan), p=%.3f (found=%d)\n',...
                    cfg.freqs(f,1),cfg.freqs(f,2),cfg.clusTimes(t,1),cfg.clusTimes(t,2),cell2mat(cfg.this_roi),vs_str,sigpos(iPos),...
                    sum(ismember(cfg.channel,data.(exper.sesStr{sesNum}).(cond1).label(sigposCLM(:,iPos) > 0))),...
                    stat_clus.(vs_str).posclusters(iPos).prob,foundpos(r,t));
                  
                  clus_symb = cfg.clusSymb(find(stat_clus.(vs_str).posclusters(iPos).prob < cfg.clusSize,1,'first'));
                  if length(cond1) < 7 && length(cond2) < 7
                    clus_str = sprintf('%s>%s:%s',cond1,cond2,clus_symb);
                  else
                    clus_str = sprintf('%d>%d:%s',cond1ind,cond2ind,clus_symb);
                  end
                  h_t = text(mean(cfg.clusTimes(t,:),2) - (textSpaceHorz * 2),max(reshape(dataVec(:,r,t),[],1)) + (textSpaceVert * foundpos(r,t)),...
                    clus_str);
                  set(h_t,'FontSize',cfg.textFontSize);
                end
              end % iPos
            end
            
            % find the negative clusters
            if ~isempty(signeg)
              for iNeg = 1:length(signeg)
                
                if sum(ismember(cfg.channel,data.(exper.sesStr{sesNum}).(cond2).label(signegCLM(:,iNeg) > 0))) > 0
                  %fprintf('%s:\t***Found negative clusters***\n',vs_str);
                  foundneg(r,t) = foundneg(r,t) + 1;
                  
                  % debug
                  fprintf('%.1fto%.1fHz: %.2fto%.2f sec, %5s, %s: Neg clus #%d (%d chan), p=%.3f (found=%d)\n',...
                    cfg.freqs(f,1),cfg.freqs(f,2),cfg.clusTimes(t,1),cfg.clusTimes(t,2),cell2mat(cfg.this_roi),vs_str,signeg(iNeg),...
                    sum(ismember(cfg.channel,data.(exper.sesStr{sesNum}).(cond2).label(signegCLM(:,iNeg) > 0))),...
                    stat_clus.(vs_str).negclusters(iNeg).prob,foundneg(r,t));
                  
                  clus_symb = cfg.clusSymb(find(stat_clus.(vs_str).negclusters(iNeg).prob < cfg.clusSize,1,'first'));
                  if length(cond1) < 7 && length(cond2) < 7
                    clus_str = sprintf('%s<%s:%s',cond1,cond2,clus_symb);
                  else
                    clus_str = sprintf('%d<%d:%s',cond1ind,cond2ind,clus_symb);
                  end
                  h_t = text(mean(cfg.clusTimes(t,:),2) - (textSpaceHorz * 2),min(reshape(dataVec(:,r,t),[],1)) - (textSpaceVert * foundneg(r,t)),...
                    clus_str);
                  set(h_t,'FontSize',cfg.textFontSize);
                end
              end % iNeg
            end
            
            %hold off
          end % r
        end % evVal
      end % t
      
      sigElecsAcrossComparisons.pos{typ,f} = sigElecsAcrossComparisons.pos{typ,f}(:,validComparisons{typ});
      sigElecsAcrossComparisons.neg{typ,f} = sigElecsAcrossComparisons.neg{typ,f}(:,validComparisons{typ});
    end % if plotClusSig
    
    if ~isfield(files,'figFontName')
      files.figFontName = 'Helvetica';
    end
    
    if files.saveFigs
      %publishfig(gcf,~cfg.plotTitle,cfg.ticFontSize,cfg.labelFontSize,files.figFontName);
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
      
      cfg.figfilename = sprintf('tfr_%s_ga_%s%s%d_%d_%d_%d%s%s%s%s%s',cfg.type,type_str,cond_str,round(cfg.freqs(f,1)),round(cfg.freqs(f,2)),round(cfg.xlim(1)*1000),round(cfg.xlim(2)*1000),chan_str_all,alpha_str,cfg.legend_str,cfg.title_str,cfg.eb_str);
      
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
      
      if cfg.nRow > 2
        % make it taller if necessary
        cfg.pos = get(gcf, 'Position');
        % get the height x width ratio
        hwRatio = cfg.pos(3) / cfg.pos(4);
        % % square figure
        % cfg.figSize = [ceil(min(cfg.screenXY) * 0.85) ceil(min(cfg.screenXY) * 0.85)];
        % maintain figure height x width ratio
        cfg.figSize = [ceil(min(cfg.screenXY) * 0.85) ceil(min(cfg.screenXY) * 0.85 * hwRatio)];
        % resize the figure window
        set(gcf, 'Units', 'pixels', 'Position', [ceil(cfg.pos(1) * 0.6), cfg.pos(2), ceil(cfg.figSize(2) / floor(cfg.nRow/2)), cfg.figSize(1)]);
      end
      %print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsTFR,cfg.figfilename));
      screen2file(fullfile(dirs.saveDirFigsTFR,cfg.figfilename),files);
    end
    
    if ~files.saveFigs
      publishfig(gcf,~cfg.plotTitle,cfg.ticFontSize,cfg.labelFontSize,files.figFontName);
      if exist('tightfig','file')
        tightfig(gcf);
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
      if cfg.nRow > 2
        set(gcf, 'Units', 'pixels', 'Position', [ceil(cfg.pos(1) * 0.6), cfg.pos(2), ceil(cfg.figSize(2) / floor(cfg.nRow/2)), cfg.figSize(1)]);
      else
        set(gcf, 'Units', 'pixels', 'Position', [ceil(cfg.pos(1) * 0.6), cfg.pos(2), cfg.figSize(2), cfg.figSize(1)]);
      end
    else
      close all
    end
    
  end % frq
end % typ

end
