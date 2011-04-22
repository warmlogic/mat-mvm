function mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs)
%mm_ft_clusterplotTFR Plot (and save) significant clusters
%   

close all

% p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
cfg_ft.highlightsymbolseries = ['*','x','+','o','.'];
cfg_ft.highlightcolorpos = [0.5 0 1];
cfg_ft.highlightcolorneg = [0 0.5 0];
cfg_ft.elec = ana.elec;
cfg_ft.contournum = 0;
cfg_ft.emarker = '.';
cfg_ft.zparam = 'stat';
if ~isfield(cfg_ft,'alpha')
  cfg_ft.alpha  = 0.05;
end
if ~isfield(cfg_ft,'zlim')
  cfg_ft.zlim = [-5 5];
end

% make sure cfg_plot.conditions is set correctly
if ~isfield(cfg_plot,'condMethod')
  if ~iscell(cfg_plot.conditions) && (strcmp(cfg_plot.conditions,'all') || strcmp(cfg_plot.conditions,'all_across_types') || strcmp(cfg_plot.conditions,'all_within_types'))
    cfg_plot.condMethod = 'pairwise';
  elseif iscell(cfg_plot.conditions) && ~iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
    cfg_plot.condMethod = 'pairwise';
  elseif iscell(cfg_plot.conditions) && iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions{1}) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
    cfg_plot.condMethod = 'pairwise';
  else
    cfg_plot.condMethod = [];
  end
end
cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);

% set the directory to load the file from
dirs.saveDirClusStat = fullfile(dirs.saveDir,sprintf('tfr_stat_clus_%d_%d',cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000));

for cnd = 1:length(cfg_plot.conditions)
  % set the number of conditions that we're testing
  cfg_plot.numConds = size(cfg_plot.conditions{cnd},2);
  vs_str = sprintf('%s%s',cfg_plot.conditions{cnd}{1},sprintf(repmat('vs%s',1,cfg_plot.numConds-1),cfg_plot.conditions{cnd}{2:end}));
  
  fprintf('%s, %d--%d ms, %.1f--%.1f Hz\n',vs_str,cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000,cfg_ft.frequency(1),cfg_ft.frequency(2));
  
  savedFile = fullfile(dirs.saveDirClusStat,sprintf('tfr_stat_clus_%s_%d_%d_%d_%d.mat',vs_str,cfg_ft.frequency(1),cfg_ft.frequency(2),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000));
  if exist(savedFile,'file')
    fprintf('Loading %s\n',savedFile);
    load(savedFile);
  else
    fprintf('No stat_clus file found for %s. Going to next comparison.\n',vs_str);
    continue
  end
  
  if ~isfield(stat_clus.(vs_str),'posclusters') && ~isfield(stat_clus.(vs_str),'negclusters')
    fprintf('%s:\tNo positive or negative clusters\n',vs_str);
    continue
  end
  
  if isfield(stat_clus.(vs_str),'posclusters') || isfield(stat_clus.(vs_str),'negclusters')
    if ~isempty(stat_clus.(vs_str).posclusters)
      %for i = 1:length(stat_clus.(vs_str).posclusters)
      %  fprintf('%s, Pos %d = %.5f\n',vs_str,i,stat_clus.(vs_str).posclusters(i).prob);
      %end
      fprintf('%s\tSmallest Pos = %.5f\n',vs_str,stat_clus.(vs_str).posclusters(1).prob);
    end
    if ~isempty(stat_clus.(vs_str).negclusters)
      %for i = 1:length(stat_clus.(vs_str).negclusters)
      %  fprintf('%s, Neg %d = %.5f\n',vs_str,i,stat_clus.(vs_str).negclusters(i).prob);
      %end
      fprintf('%s\tSmallest Neg = %.5f\n',vs_str,stat_clus.(vs_str).negclusters(1).prob);
    end
    
    %if ~isempty(stat_clus.(vs_str).posclusters) || ~isempty(stat_clus.(vs_str).negclusters)
    
    if (~isempty(stat_clus.(vs_str).posclusters) && ~isempty(find(stat_clus.(vs_str).posclusters(1).prob < cfg_ft.alpha,1))) || (~isempty(stat_clus.(vs_str).negclusters) && ~isempty(find(stat_clus.(vs_str).negclusters(1).prob < cfg_ft.alpha,1)))
      fprintf('%s:\t***Found positive or negative clusters***\n',vs_str);
      ft_clusterplot(cfg_ft,stat_clus.(vs_str));
      
      if files.saveFigs
        fignums = findobj('Type','figure');
        for f = 1:length(fignums)
          figure(f)
          
          cfg_plot.figfilename = sprintf('tfr_clus_ga_%s_%d_%d_%d_%d_fig%d.%s',vs_str,round(cfg_ft.frequency(1)),round(cfg_ft.frequency(2)),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000,f,files.figFileExt);
          
          dirs.saveDirFigsClus = fullfile(dirs.saveDirFigs,sprintf('tfr_stat_clus_%d_%d',cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000),vs_str);
          %dirs.saveDirFigsClus = fullfile(dirs.saveDirFigs,'tfr_stat_clus',vs_str);
          if ~exist(dirs.saveDirFigsClus,'dir')
            mkdir(dirs.saveDirFigsClus)
          end
          print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigsClus,cfg_plot.figfilename));
        end
        close all
      else
        keyboard
        close all
      end % if
      
    elseif (~isempty(stat_clus.(vs_str).posclusters) && isempty(find(stat_clus.(vs_str).posclusters(1).prob < cfg_ft.alpha,1))) || (~isempty(stat_clus.(vs_str).negclusters) && isempty(find(stat_clus.(vs_str).negclusters(1).prob < cfg_ft.alpha,1)))
      fprintf('%s:\tNo significant positive or negative clusters\n',vs_str);
    elseif isempty(stat_clus.(vs_str).posclusters) && isempty(stat_clus.(vs_str).negclusters)
      fprintf('%s:\tNo positive or negative clusters\n',vs_str);
    end
  end % if isfield
  
end % for cnd

end
