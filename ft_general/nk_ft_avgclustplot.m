function outdata = nk_ft_avgclustplot(stat_clus,cfg_plot,cfg_ft,dirs,files,savefigs)

if ~isempty(savefigs)
  files.saveFigs = 1;
end

if ~isfield(cfg_ft,'alpha')
  cfg_ft.alpha = 'all';
end
if isequal(cfg_ft.alpha,'all')
  cfg_ft.alpha = 1.0;
end

if ~isfield(cfg_plot,'isfreq')
  cfg_plot.isfreq = true;
end

vs_str = fieldnames(stat_clus);%[cfg.conds{1} 'vs' cfg.conds{2}];
vs_str = vs_str{1};

if ~isfield(stat_clus.(vs_str),'posclusters') && ~isfield(stat_clus.(vs_str),'negclusters')
  fprintf('%s:\tNo positive or negative clusters found.\n',vs_str);
  return
end

load('dummy1dclus.mat');

% for denoting directionality of found clusters
plot_clus_str = {};

% find how many clusters were below cfg_ft.alpha
if isfield(stat_clus.(vs_str),'posclusters') || isfield(stat_clus.(vs_str),'negclusters')
  if ~isempty(stat_clus.(vs_str).posclusters)
    %for i = 1:length(stat_clus.(vs_str).posclusters)
    %  fprintf('%s, Pos (%d of %d) p=%.5f\n',vs_str,i,length(stat_clus.(vs_str).posclusters),stat_clus.(vs_str).posclusters(i).prob);
    %end
    fprintf('%s\tSmallest Pos: p=%.5f\n',vs_str,stat_clus.(vs_str).posclusters(1).prob);
  end
  if ~isempty(stat_clus.(vs_str).negclusters)
    %for i = 1:length(stat_clus.(vs_str).negclusters)
    %  fprintf('%s, Neg (%d of %d) p=%.5f\n',vs_str,i,length(stat_clus.(vs_str).negclusters),stat_clus.(vs_str).negclusters(i).prob);
    %end
    fprintf('%s\tSmallest Neg: p=%.5f\n',vs_str,stat_clus.(vs_str).negclusters(1).prob);
  end
  
  if ~isempty(stat_clus.(vs_str).posclusters) || ~isempty(stat_clus.(vs_str).negclusters)
    sigpos = [];
    if ~isempty(stat_clus.(vs_str).posclusters)
      for iPos = 1:length(stat_clus.(vs_str).posclusters)
        sigpos(iPos) = stat_clus.(vs_str).posclusters(iPos).prob <= cfg_ft.alpha;
      end
      sigpos = find(sigpos == 1);
    end
    signeg = [];
    if ~isempty(stat_clus.(vs_str).negclusters)
      for iNeg = 1:length(stat_clus.(vs_str).negclusters)
        signeg(iNeg) = stat_clus.(vs_str).negclusters(iNeg).prob <= cfg_ft.alpha;
      end
      signeg = find(signeg == 1);
    end
    Nsigpos = length(sigpos);
    Nsigneg = length(signeg);
    Nsigall = Nsigpos + Nsigneg;
    
    clus_str = '';
    if Nsigpos > 0
      clus_str = cat(2,clus_str,'positive');
      plot_clus_str{1} = 'pos';
    end
    if Nsigneg > 0 && isempty(clus_str)
      clus_str = cat(2,clus_str,'negative');
      plot_clus_str{1} = 'neg';
    elseif Nsigneg > 0 && ~isempty(clus_str)
      clus_str = cat(2,clus_str,' and negative');
      plot_clus_str{2} = 'neg';
    end
    
    if Nsigall > 0
      if Nsigall == 1
        clus_str = cat(2,clus_str,' cluster');
      elseif Nsigall > 1
        clus_str = cat(2,clus_str,' clusters');
      end
      fprintf('%s:\t***Found significant %s at p<%.3f***\n',vs_str,clus_str,cfg_ft.alpha);
    else
      warning('No clusters found at p<%.5f',cfg_ft.alpha);
      return
    end
  end
end


for cl = 1:length(plot_clus_str)
  for Nsig = 1:eval(sprintf('Nsig%s',plot_clus_str{cl}))
    cfg_ft.clusnum = Nsig;
    
    ind = find(stat_clus.(vs_str).(sprintf('%sclusterslabelmat',plot_clus_str{cl}))==cfg_ft.clusnum);
    [x,y,z] = ind2sub(size(stat_clus.(vs_str).(sprintf('%sclusterslabelmat',plot_clus_str{cl}))),ind);
    elecs = stat_clus.(vs_str).label(unique(x));
    elecnums = unique(x);
    
    tvals = zeros(length(stat_clus.(vs_str).label),1,1);
    for i = 1:length(elecnums)
      ind = find(x==elecnums(i));
      tempdata = [];
      for j = 1:length(ind)
        tempdata(j) = stat_clus.(vs_str).stat(x(ind(j)),y(ind(j)),z(ind(j)));
      end
      tvals(elecnums(i)) = sum(tempdata);
    end
    
    %new_clus = stat_clus.(vs_str);
    new_clus.posclusters = [];
    new_clus.negclusters = [];
    new_clus.stat = tvals;
    new_clus.prob = tvals;
    new_clus.label = stat_clus.(vs_str).label;
    lmat = zeros(size(tvals));
    lmat(elecnums) = 1;
    new_clus.(sprintf('%sclusterslabelmat',plot_clus_str{cl})) = lmat;
    new_clus.(sprintf('%sclusters',plot_clus_str{cl})) = stat_clus.(vs_str).(sprintf('%sclusters',plot_clus_str{cl}));
    cfg_ft.highlightcolorpos = [1 1 1];
    cfg_ft.highlightcolorneg = [1 1 1];
    %cfg_ft.highlightseries = {'numbers','numbers','numbers','numbers','numbers'};
    ft_clusterplot(cfg_ft,new_clus);
    colorbar;
    title([regexprep(vs_str,'_','') ', Cluster ' num2str(cfg_ft.clusnum) ', p=' sprintf('%.5f',stat_clus.(vs_str).(sprintf('%sclusters',plot_clus_str{cl}))(cfg_ft.clusnum).prob)],'fontsize',20);
    outdata = tvals;
    
    if ~isfield(cfg_plot,'dirStr')
      cfg_plot.dirStr = '';
    end
    p = stat_clus.(vs_str).(sprintf('%sclusters',plot_clus_str{cl}))(cfg_ft.clusnum).prob;
    if files.saveFigs
      %     fignums = findobj('Type','figure');
      %     for f = 1:length(fignums)
      %         figure(f)
      f = cfg_ft.clusnum;
      p_str = strrep(sprintf('%.3f',p),'.','p');
      if cfg_plot.isfreq
        cfg_plot.figfilename = sprintf('tfr_clus_avgclus_%s_%d_%d_%s%d_%d_%d_%s',vs_str,round(cfg_plot.frequency(1)),round(cfg_plot.frequency(2)),plot_clus_str{cl},f,round(cfg_plot.latency(1)*1000),round(cfg_plot.latency(2)*1000),p_str);
        dirs.saveDirFigsClus = fullfile(dirs.saveDirFigs,sprintf('tfr_stat_clus_%d_%d%s',round(cfg_plot.latency(1)*1000),round(cfg_plot.latency(2)*1000),cfg_plot.dirStr),vs_str);
      else
        cfg_plot.figfilename = sprintf('tla_clus_avgclus_%s_%s%d_%d_%d_%s',vs_str,plot_clus_str{cl},f,round(cfg_plot.latency(1)*1000),round(cfg_plot.latency(2)*1000),p_str);
        dirs.saveDirFigsClus = fullfile(dirs.saveDirFigs,sprintf('tla_stat_clus_%d_%d%s',round(cfg_plot.latency(1)*1000),round(cfg_plot.latency(2)*1000),cfg_plot.dirStr),vs_str);
      end
      if ~exist(dirs.saveDirFigsClus,'dir')
        mkdir(dirs.saveDirFigsClus)
      end
      
      %while exist([fullfile(dirs.saveDirFigsClus,cfg_plot.figfilename),'.',files.figPrintFormat],'file')
      %  f=f+1;
      %  cfg_plot.figfilename = sprintf('tfr_clus_avgclus_%s_%d_%d_%d_%d_fig%d',vs_str,round(cfg_plot.frequency(1)),round(cfg_plot.frequency(2)),round(cfg_plot.latency(1)*1000),round(cfg_plot.latency(2)*1000),f);
      %end
      
      if strcmp(files.figPrintFormat,'fig')
        saveas(gcf,fullfile(dirs.saveDirFigsClus,[cfg_plot.figfilename,'.',files.figPrintFormat]),'fig');
      else
        if strcmp(files.figPrintFormat(1:2),'-d')
          files.figPrintFormat = files.figPrintFormat(3:end);
        end
        if ~isfield(files,'figPrintRes')
          files.figPrintRes = 150;
        end
        set(gcf,'InvertHardCopy','off');
        set(gcf,'Color','White');
        print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsClus,[cfg_plot.figfilename,'.',files.figPrintFormat]));
      end
    end % if
    
  end
end

end