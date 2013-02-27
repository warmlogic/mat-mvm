function [outdata] = nk_ft_avgpowerbytime(freqdata,stat_clus,cfg_plot,cfg,dirs,files,savefigs)

% create a 2d power as a function of time plot for selected electrodes,
% frequencies, and conditions
%
% input
%   freqdata = power data for individual subjects
%   cfg = plot specifications
%       cfg.time = time points to display
%       cfg.freq = frequencies average over
%       cfg.elecs = cell array of electrodes to average over
%       cfg.conds = cell array of conditions names to display
%
% output
%   outdata = avgerage data plotted
if isempty(savefigs)
  files.saveFigs = 1;
end
if ~isfield(cfg,'transp')
  cfg.transp = 0;
end

if ~isfield(cfg,'alpha')
  cfg.alpha = 'all';
end
if isequal(cfg.alpha,'all')
  cfg.alpha = 1.0;
end

vs_str = fieldnames(stat_clus);
vs_str = vs_str{1};

if ~isfield(stat_clus.(vs_str),'posclusters') && ~isfield(stat_clus.(vs_str),'negclusters')
  fprintf('%s:\tNo positive or negative clusters found.\n',vs_str);
  return
end

% for denoting directionality of found clusters
plot_clus_str = {};

% find how many clusters were below cfg.alpha
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
        sigpos(iPos) = stat_clus.(vs_str).posclusters(iPos).prob <= cfg.alpha;
      end
      sigpos = find(sigpos == 1);
    end
    signeg = [];
    if ~isempty(stat_clus.(vs_str).negclusters)
      for iNeg = 1:length(stat_clus.(vs_str).negclusters)
        signeg(iNeg) = stat_clus.(vs_str).negclusters(iNeg).prob <= cfg.alpha;
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
      fprintf('%s:\t***Found significant %s at p<%.3f***\n',vs_str,clus_str,cfg.alpha);
    else
      warning('No clusters found at p<%.5f',cfg.alpha);
      return
    end
  end
end

for cl = 1:length(plot_clus_str)
  for Nsig = 1:eval(sprintf('Nsig%s',plot_clus_str{cl}))
    cfg.clusnum = Nsig;
    
    ind = find(stat_clus.(vs_str).(sprintf('%sclusterslabelmat',plot_clus_str{cl}))==cfg.clusnum);
    [x,y,z] = ind2sub(size(stat_clus.(vs_str).(sprintf('%sclusterslabelmat',plot_clus_str{cl}))),ind);
    elecs = stat_clus.(vs_str).label(unique(x));
    sigt = stat_clus.(vs_str).time(unique(z));
    if isfield(stat_clus.(vs_str),'freq')
      sigf = stat_clus.(vs_str).freq(unique(y));
      cfg.freq = sigf;%stat_clus.(vs_str).cfg.frequency;%cfg_ana.frequencies;%[8 12];
    else
      sigf = stat_clus.(vs_str).cfg.frequency;
      cfg.freq = stat_clus.(vs_str).cfg.frequency;
    end
    
    cfg.elecs = {};
    for i = 1:length(elecs)
      cfg.elecs = cat(2,cfg.elecs,elecs{i});
    end
    
    cfg.sigt = sigt;
    
    fprintf('%s: p=%f, t=%f to %f\n',vs_str,stat_clus.(vs_str).(sprintf('%sclusters',plot_clus_str{cl}))(cfg.clusnum).prob,cfg.sigt);
    
    %make condition average data matrix across subjects
    conddata = [];
    for icond = 1:length(cfg.conds)
      inconds = fieldnames(freqdata);
      condidx = strcmp(inconds,cfg.conds{icond});
      if sum(condidx) ~= 1
        error('nonexistant or nonunique condition name');
      else
        cond = inconds{condidx};
      end
      
      %average condition data across subjects
      tempdata = [];
      for isub = 1:length(freqdata.(cond).sub)
        data = ft_selectdata(freqdata.(cond).sub(isub).ses(1).data, 'foilim',cfg.freq,'avgoverfreq','yes');
        data = ft_selectdata(data, 'channel', cfg.elecs,'avgoverchan','yes');
        data = ft_selectdata(data, 'toilim',cfg.time,'avgovertime','no');
        tempdata(:,isub) = squeeze(data.powspctrm);
      end
      conddata(:,icond) = mean(tempdata,2);
      conddata_var(:,icond) = nanste(tempdata,2)';
    end
    
    outdata.cfg = cfg;
    outdata.data = conddata;
    outdata.conds = cfg.conds;
    outdata.time = data.time;
    outdata.var = conddata_var;
    
    
    figure('color','white');
    %plot(outdata.time,outdata.data,'linewidth',5);
    if ~isfield(cfg_plot,'colors')
      colors = get(gca,'ColorOrd');
    else
      colors = cfg_plot.colors;
    end
    for i = 1:size(outdata.data,2)
      h = shadedErrorBar(outdata.time,outdata.data(:,i),outdata.var(:,i),{'color',colors(i,:),'linewidth',5},cfg.transp);
      lh(i) = h.mainLine;
      hold on
    end
    
    
    set(gca,'fontsize',22);
    %xlim([-.5,max(data.time)]);
    xlim(cfg.time);
    xlabel('Time (s)');
    ylabel('Power');
    %outdata.conds{3} = 'pB';outdata.conds{2} = 'T';
    %legend(lh,regexprep(outdata.conds,'_',''),'location','best','fontsize',10);
    title(sprintf('%s',['Cluster ' num2str(cfg.clusnum) ', AvgPwr' num2str(cfg.freq(1)) 'to' num2str(cfg.freq(end)) 'Hz']));
    if isfield(cfg,'sigt')
      plot(repmat(min(cfg.sigt),2,1),ylim,'--k','linewidth',3);
      hold on
      plot(repmat(max(cfg.sigt),2,1),ylim,'--k','linewidth',3);
    end
    
    plot(xlim,[0 0],'k--'); % horizontal
    plot([0 0],ylim,'k--'); % vertical
    
    box off
    
    if ~isfield(cfg_plot,'dirStr')
      cfg_plot.dirStr = '';
    end
    p = stat_clus.(vs_str).(sprintf('%sclusters',plot_clus_str{cl}))(cfg.clusnum).prob;
    if files.saveFigs
      %fignums = findobj('Type','figure');
      %for f = 1:length(fignums)
      %figure(f)
      f=cfg.clusnum;
      p_str = strrep(sprintf('%.3f',p),'.','p');
      cfg_plot.figfilename = sprintf('tfr_clus_avgfreq_%s_%d_%d_%s%d_%d_%d_%s',vs_str,round(sigf(1)),round(sigf(2)),plot_clus_str{cl},f,round(sigt(1)*1000),round(sigt(end)*1000),p_str);
      
      dirs.saveDirFigsClus = fullfile(dirs.saveDirFigs,sprintf('tfr_stat_clus_%d_%d%s',round(cfg_plot.latency(1)*1000),round(cfg_plot.latency(2)*1000),cfg_plot.dirStr),vs_str);
      if ~exist(dirs.saveDirFigsClus,'dir')
        mkdir(dirs.saveDirFigsClus)
      end
      
      %         while exist([fullfile(dirs.saveDirFigsClus,cfg_plot.figfilename) '.' files.figPrintFormat],'file')
      %             f=f+1;
      %             cfg_plot.figfilename = sprintf('tfr_clus_avgclus_%s_%d_%d_%d_%d_fig%d',vs_str,round(cfg_plot.frequencies(1)),round(cfg_plot.frequencies(2)),round(cfg_plot.latencies(1)*1000),round(cfg_plot.latencies(2)*1000),f);
      %         end
      if strcmp(files.figPrintFormat,'fig')
        saveas(gcf,fullfile(dirs.saveDirFigsClus,[cfg_plot.figfilename,'.',files.figPrintFormat]),'fig');
      else
        if strcmp(files.figPrintFormat(1:2),'-d')
          files.figPrintFormat = files.figPrintFormat(3:end);
        end
        if ~isfield(files,'figPrintRes')
          files.figPrintRes = 300;
        end
        print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigsClus,[cfg_plot.figfilename,'.',files.figPrintFormat]));
      end
    end % if
    
  end % Nsig
end % plot_clus_str