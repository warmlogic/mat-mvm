function mm_ft_subjplotTFR(cfg_ft,cfg_plot,ana,exper,data)
%MM_FT_SUBPLOTTFR make a subplot of the individual subject data
%   

% FIXME: implement support for multiple sessions
ses = 1;

vertTextLoc = 0.18;

% make sure some fields are set
if ~isfield(cfg_ft,'zparam')
  error('Must define cfg_ft.zparam for ft_singleplotTFR');
end
if ~isfield(cfg_plot,'excludeBadSub')
  cfg_plot.excludeBadSub = 0;
end
if ~isfield(cfg_plot,'numCols')
  cfg_plot.numCols = 5;
end

% make sure conditions are set correctly
if ~isfield(cfg_plot,'condMethod')
  if ~iscell(cfg_plot.conditions) && (strcmp(cfg_plot.conditions,'all') || strcmp(cfg_plot.conditions,'all_across_types') || strcmp(cfg_plot.conditions,'all_within_types'))
    cfg_plot.condMethod = 'single';
  elseif iscell(cfg_plot.conditions) && ~iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
    cfg_plot.condMethod = 'single';
  elseif iscell(cfg_plot.conditions) && iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions{1}) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
    cfg_plot.condMethod = 'single';
  else
    cfg_plot.condMethod = [];
  end
end
cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);
% make sure conditions are set up for the for loop
if ~isfield(cfg_plot,'types')
  cfg_plot.types = repmat({''},size(cfg_plot.conditions));
end

% set the channel information
if ~isfield(cfg_plot,'roi')
  error('Must specify either ROI names or channel names in cfg_plot.roi');
elseif isfield(cfg_plot,'roi')
  if ismember(cfg_plot.roi,ana.elecGroupsStr)
    % if it's in the predefined ROIs, get the channel numbers
    cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
  else
    % otherwise it should be the channel number(s) or 'all'
    if ~iscell(cfg_plot.roi)
      cfg_plot.roi = {cfg_plot.roi};
    end
    cfg_ft.channel = cfg_plot.roi;
  end
end

% set the number of columns and rows
if cfg_plot.numCols > length(exper.subjects) + 1
  cfg_plot.numCols = length(exper.subjects) + 1;
end
if cfg_plot.excludeBadSub == 1
  cfg_plot.numRows = ceil((length(exper.subjects) - sum(exper.badSub)) / cfg_plot.numCols);
else
  cfg_plot.numRows = ceil((length(exper.subjects)) / cfg_plot.numCols);
end
if (cfg_plot.numCols * cfg_plot.numRows) == length(exper.subjects)
  cfg_plot.numRows = cfg_plot.numRows + 1;
end

for typ = 1:length(cfg_plot.conditions)
  for evVal = 1:length(cfg_plot.conditions{typ})
    figure
    for sub = 1:length(exper.subjects)
      if cfg_plot.excludeBadSub == 1 && exper.badSub(sub,ses)
        fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
        continue
      else
        subplot(cfg_plot.numRows,cfg_plot.numCols,sub);
        
        if isfield(data.(cfg_plot.conditions{typ}{evVal}).sub(sub).ses(ses).data,cfg_ft.zparam)
          ft_singleplotTFR(cfg_ft,data.(cfg_plot.conditions{typ}{evVal}).sub(sub).ses(ses).data);
        end
        
        if cfg_plot.excludeBadSub == 0 && exper.badSub(sub,ses)
          fprintf('Found bad subject: %s\n',exper.subjects{sub});
          title(sprintf('*BAD*: %s;%s',exper.subjects{sub},num2str(exper.nTrials.(cfg_plot.conditions{typ}{evVal})(sub))));
        else
          title(sprintf('%s;%s',exper.subjects{sub},num2str(exper.nTrials.(cfg_plot.conditions{typ}{evVal})(sub))));
        end
      end
    end
    
    % give some info
    subplot(cfg_plot.numRows,cfg_plot.numCols,length(exper.subjects) + 1);
    ycoord = 1.0;
    if ~isempty(cfg_plot.types{typ})
      text(0.5,ycoord,sprintf('%s',cfg_plot.types{typ}),'color','k');
      ycoord = ycoord - vertTextLoc;
      set(gcf,'Name',sprintf('%s, %s,%s',cfg_plot.types{typ},cfg_plot.conditions{typ}{evVal},sprintf(repmat(' %s',1,length(cfg_plot.roi)),cfg_plot.roi{:})))
    else
      set(gcf,'Name',sprintf('%s,%s',cfg_plot.conditions{typ}{evVal},sprintf(repmat(' %s',1,length(cfg_plot.roi)),cfg_plot.roi{:})))
    end
    text(0.5,ycoord,sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}),'color','k');
    ycoord = ycoord - vertTextLoc;
    text(0.5,ycoord,cfg_plot.conditions{typ}{evVal},'color','k');
    axis off
  end
end

end
