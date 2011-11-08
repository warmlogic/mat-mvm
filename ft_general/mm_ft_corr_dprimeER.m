function mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data)
%MM_FT_CORR_DPRIMEER Correlate voltage differences with dprime
%

if ~isfield(cfg_ana,'parameter')
  error('Must specify cfg_ana.parameter, denoting the data to test (e.g., ''avg'' or ''individual'')');
end

% get the label info for this data struct
if isfield(data.(exper.eventValues{1}).sub(1).ses(1).data,'label');
  lab = data.(exper.eventValues{1}).sub(1).ses(1).data.label;
else
  error('label information not found in data struct');
end

% make sure cfg_ana.conditions is set correctly
if ~isfield(cfg_ana,'condMethod') || isempty(cfg_ana.condMethod)
  if ~iscell(cfg_ana.conditions) && (strcmp(cfg_ana.conditions,'all') || strcmp(cfg_ana.conditions,'all_across_types') || strcmp(cfg_ana.conditions,'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  elseif iscell(cfg_ana.conditions) && ~iscell(cfg_ana.conditions{1}) && length(cfg_ana.conditions) == 1 && (strcmp(cfg_ana.conditions{1},'all') || strcmp(cfg_ana.conditions{1},'all_across_types') || strcmp(cfg_ana.conditions{1},'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  elseif iscell(cfg_ana.conditions) && iscell(cfg_ana.conditions{1}) && length(cfg_ana.conditions{1}) == 1 && (strcmp(cfg_ana.conditions{1},'all') || strcmp(cfg_ana.conditions{1},'all_across_types') || strcmp(cfg_ana.conditions{1},'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  else
    cfg_ana.condMethod = [];
  end
end
cfg_ana.conditions = mm_ft_checkConditions(cfg_ana.conditions,ana,cfg_ana.condMethod);
% make sure conditions are set up for the for loop
if ~isfield(cfg_ana,'types')
  cfg_ana.types = repmat({''},size(cfg_ana.conditions));
end

% set up how the lines will look
cfg_plot = [];
cfg_plot.linespec = 'ko';
cfg_plot.linewidth = 1.5;
cfg_plot.marksize = 8;
cfg_plot.markcolor = 'w';
cfg_plot.linefit = 'r-';

% set the channel information
if ~isfield(cfg_ana,'roi')
  error('Must specify either ROI names or channel names in cfg_ana.roi');
elseif isfield(cfg_ana,'roi')
  if ismember(cfg_ana.roi,ana.elecGroupsStr)
    % if it's in the predefined ROIs, get the channel numbers
    cfg_ana.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi)});
    % find the channel indices for averaging
    cfg_ana.chansel = ismember(lab,cfg_ana.channel);
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:});
  else
    % otherwise it should be the channel number(s) or 'all'
    if ~iscell(cfg_ana.roi)
      cfg_ana.roi = {cfg_ana.roi};
    end
    cfg_ana.channel = cfg_ana.roi;
    
    % find the channel indices for averaging
    if strcmp(cfg_ana.channel,'all')
      cfg_ana.chansel = ismember(lab,ft_channelselection(cfg_ana.channel,lab));
    else
      cfg_ana.chansel = ismember(lab,cfg_ana.channel);
    end
    % set the string for the filename
    cfg_plot.chan_str = sprintf(repmat('%s_',1,length(cfg_ana.channel)),cfg_ana.channel{:});
  end
end

% exclude the bad subjects from the subject count
cfg_ana.numSub = length(exper.subjects) - sum(exper.badSub);

% get the conditions for this type
for typ = 1:length(cfg_ana.conditions)
  % which d' values are we dealing with (e.g., item or source)
  for dp = 1:length(cfg_ana.dpTypes)
    
    % initialize for storing the voltage values
    for evVal = 1:length(cfg_ana.conditions{typ}{dp})
      cfg_ana.values.(cfg_ana.conditions{typ}{dp}{evVal}) = [];
    end
    % get the voltage for the correct region
    for evVal = 1:length(cfg_ana.conditions{typ}{dp})
      ev = cfg_ana.conditions{typ}{dp}{evVal};
      cfg_ana.values.(ev) = nan(cfg_ana.numSub,length(exper.sessions));
      goodSubInd = 0;
      for sub = 1:length(exper.subjects)
        ses = 1;
        %for ses = 1:length(exper.sessions)
        if exper.badSub(sub,ses)
          fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
          continue
        else
          goodSubInd = goodSubInd + 1;
          cfg_ana.timesel.(ev) = find(data.(ev).sub(sub).ses(ses).data.time >= cfg_ana.latency(1) & data.(ev).sub(sub).ses(ses).data.time <= cfg_ana.latency(2));
          cfg_ana.values.(ev)(goodSubInd,ses) = mean(mean(data.(ev).sub(sub).ses(ses).data.(cfg_ana.parameter)(cfg_ana.chansel,cfg_ana.timesel.(ev)),1),2);
        end
        %end % ses
      end % sub
      
      %cfg_ana.timesel.(ev) = find(data.(ev).sub(sub).ses(ses).data.time >= cfg_ana.latency(1) & data.(ev).sub(sub).ses(ses).data.time <= cfg_ana.latency(2));
      %cfg_ana.timesel.(cfg_ana.conditions{typ}{dp}{evVal}) = find(data.(cfg_ana.conditions{typ}{dp}{evVal}).time >= cfg_ana.latency(1) & data.(cfg_ana.conditions{typ}{dp}{evVal}).time <= cfg_ana.latency(2));
      % data does not include bad subjects
      %cfg_ana.values.(cfg_ana.conditions{typ}{dp}{evVal}) = mean(mean(data.(cfg_ana.conditions{typ}{dp}{evVal}).(cfg_ana.parameter)(:,cfg_ana.chansel,cfg_ana.timesel.(cfg_ana.conditions{typ}{dp}{evVal})),2),3)';
    end
    
    % d' values on the x-axis
    if ~isempty(cfg_ana.types{typ})
      x1 = cfg_ana.(cfg_ana.types{typ}).d_item(~exper.badSub);
    else
      x1 = cfg_ana.d_item(~exper.badSub);
    end
    % voltage values on the y-axis
    y1 = (cfg_ana.values.(cfg_ana.conditions{typ}{dp}{1}) - cfg_ana.values.(cfg_ana.conditions{typ}{dp}{2}));
    % run the correlation
    [rho p] = corr([x1' y1]);
    % plot it
    figure
    [m,b] = mm_ft_linefit(x1',y1,1,cfg_plot.linefit);
    hold on
    plot(x1,y1,cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
    hold off
    title(sprintf('%s,%s',cfg_ana.dpTypes{dp},sprintf(repmat(' %s',1,length(cfg_ana.roi)),cfg_ana.roi{:})));
    xlabel(sprintf('%s d''',cfg_ana.dpTypes{dp}));
    ylabel(sprintf('Voltage (%s): %s - %s','\muV',cfg_ana.conditions{typ}{dp}{1},cfg_ana.conditions{typ}{dp}{2}));
    %axis([floor(min(x1)) ceil(max(x1)) floor(min(y1)) ceil(max(y1))])
    axis([0 ceil(max(x1)) -4 4])
    axis square
    if ~isfield(files,'figFontName')
      files.figFontName = 'Helvetica';
    end
    publishfig(gcf,0,[],[],files.figFontName);
    if ~isempty(cfg_ana.types{typ})
      fprintf('%s: ',cfg_ana.types{typ});
    end
    fprintf('%s, %s d'', %s - %s:\tPearson''s r=%.4f (r^2=%.4f), p=%.4f, m=%.3f, b=%.3f\n',sprintf(repmat('%s ',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ana.dpTypes{dp},cfg_ana.conditions{typ}{dp}{1},cfg_ana.conditions{typ}{dp}{2},rho(2,1),rho(2,1)^2,p(2,1),m,b);
    if files.saveFigs
      if ~isempty(cfg_ana.types{typ})
        cfg_plot.figfilename = sprintf('tla_corr_ga_%s_%s_%s%s%d_%d',cfg_ana.types{typ},cfg_ana.dpTypes{dp},sprintf(repmat('%s_',1,length(cfg_ana.conditions{typ}{dp})),cfg_ana.conditions{typ}{dp}{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ana.latency(1)*1000,cfg_ana.latency(2)*1000);
      else
        cfg_plot.figfilename = sprintf('tla_corr_ga_%s_%s%s%d_%d',cfg_ana.dpTypes{dp},sprintf(repmat('%s_',1,length(cfg_ana.conditions{typ}{dp})),cfg_ana.conditions{typ}{dp}{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ana.latency(1)*1000,cfg_ana.latency(2)*1000);
      end
      dirs.saveDirFigsCorr = fullfile(dirs.saveDirFigs,'tla_corr');
      if ~exist(dirs.saveDirFigsCorr,'dir')
        mkdir(dirs.saveDirFigsCorr)
      end
      
      if strcmp(files.figPrintFormat(1:2),'-d')
        files.figPrintFormat = files.figPrintFormat(3:end);
      end
      if ~isfield(files,'figPrintRes')
        files.figPrintRes = 150;
      end
      print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),fullfile(dirs.saveDirFigs,dirs.saveDirFigsCorr,cfg_plot.figfilename));
    end
  end % dp
end % typ

end
