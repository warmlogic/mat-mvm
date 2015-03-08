function mm_ft_subjplotER(cfg_plot,ana,exper,data,ses)
%MM_FT_SUBPLOTER make a subplot of the individual subject data
%   

% FIXME: implement support for multiple sessions
% ses = 1;

vertTextLoc = 0.18;

if ~isfield(cfg_plot,'parameter')
  cfg_plot.parameter = 'avg';
end

% make sure some fields are set
if ~isfield(cfg_plot,'excludeBadSub')
  cfg_plot.excludeBadSub = 0;
end
if ~isfield(cfg_plot,'numCols')
  cfg_plot.numCols = 5;
end
if ~isfield(cfg_plot,'titleTrialCount')
  cfg_plot.titleTrialCount = true;
end
if ~isfield(cfg_plot,'graphcolor')
  cfg_plot.graphcolor = 'rbkgcmyrbkgcmyrbkgcmy';
end

% % make sure conditions are set correctly
% if ~isfield(cfg_plot,'condMethod')
%   if ~iscell(cfg_plot.conditions) && (strcmp(cfg_plot.conditions,'all') || strcmp(cfg_plot.conditions,'all_across_types') || strcmp(cfg_plot.conditions,'all_within_types'))
%     cfg_plot.condMethod = 'single';
%   elseif iscell(cfg_plot.conditions) && ~iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
%     cfg_plot.condMethod = 'single';
%   elseif iscell(cfg_plot.conditions) && iscell(cfg_plot.conditions{1}) && length(cfg_plot.conditions{1}) == 1 && (strcmp(cfg_plot.conditions{1},'all') || strcmp(cfg_plot.conditions{1},'all_across_types') || strcmp(cfg_plot.conditions{1},'all_within_types'))
%     cfg_plot.condMethod = 'single';
%   else
%     cfg_plot.condMethod = [];
%   end
% end
% cfg_plot.conditions = mm_ft_checkConditions(cfg_plot.conditions,ana,cfg_plot.condMethod);

% make sure conditions are set up for the for loop
if ~isfield(cfg_plot,'types')
  cfg_plot.types = repmat({''},size(cfg_plot.conditions));
end

% get the label info for this data struct
if isfield(data.(exper.sesStr{ses}).(cfg_plot.conditions{1}{1}{1}).sub(1).data,'label')
  lab = data.(exper.sesStr{ses}).(cfg_plot.conditions{1}{1}{1}).sub(1).data.label;
else
  error('label information not found in data struct\n');
end

% set the channel information
if ~isfield(cfg_plot,'roi')
  error('Must specify either ROI names or channel names in cfg_plot.roi\n');
elseif isfield(cfg_plot,'roi')
  if ismember(cfg_plot.roi,ana.elecGroupsStr)
    % if it's in the predefined ROIs, get the channel numbers
    cfg_plot.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
    % find the channel indices for averaging
    cfg_plot.chansel = ismember(lab,cfg_plot.channel);
  else
    % otherwise it should be the channel number(s) or 'all'
    if ~iscell(cfg_plot.roi)
      cfg_plot.roi = {cfg_plot.roi};
    end
    cfg_plot.channel = cfg_plot.roi;
    
    % find the channel indices for averaging
    if strcmp(cfg_plot.channel,'all')
      cfg_plot.chansel = ismember(lab,ft_channelselection(cfg_plot.channel,lab));
    else
      cfg_plot.chansel = ismember(lab,cfg_plot.channel);
    end
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

% each entry in cfg_plot.conditions is a cell containing one type of events
for typ = 1:length(cfg_plot.conditions{ses})
  figure
  for sub = 1:length(exper.subjects)
    if cfg_plot.excludeBadSub == 1 && exper.badSub(sub,ses)
      fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
      continue
    else
      subplot(cfg_plot.numRows,cfg_plot.numCols,sub);
      plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal
      hold on
      plot([0 0],[cfg_plot.ylim(1) cfg_plot.ylim(2)],'k--'); % vertical
      
      evStr = [];
      % go through each event value for this type
      for evVal = 1:length(cfg_plot.conditions{ses}{typ})
        if cfg_plot.titleTrialCount
          evStr = cat(2,evStr,sprintf(' %s:%d',strrep(cfg_plot.conditions{ses}{typ}{evVal},'_','-'),exper.nTrials.(exper.sesStr{ses}).(cfg_plot.conditions{ses}{typ}{evVal})(sub)));
        end
        if isfield(data.(exper.sesStr{ses}).(cfg_plot.conditions{ses}{typ}{evVal}).sub(sub).data,cfg_plot.parameter)
          plot(data.(exper.sesStr{ses}).(cfg_plot.conditions{ses}{typ}{evVal}).sub(sub).data.time,mean(data.(exper.sesStr{ses}).(cfg_plot.conditions{ses}{typ}{evVal}).sub(sub).data.(cfg_plot.parameter)(cfg_plot.chansel,:),1),cfg_plot.graphcolor(evVal));
        end
      end
      
      titleStr = sprintf('%s',exper.subjects{sub});
      
      if cfg_plot.excludeBadSub == 0 && exper.badSub(sub,ses)
        fprintf('Found bad subject: %s\n',exper.subjects{sub});
        titleStr = sprintf('*BAD*: %s',titleStr);
      end
      if cfg_plot.titleTrialCount
        titleStr = sprintf('%s;%s',titleStr,evStr);
      end
      title(titleStr);
      
      axis([cfg_plot.xlim(1) cfg_plot.xlim(2) cfg_plot.ylim(1) cfg_plot.ylim(2)]);
      hold off
    end
  end
  
  % give some info
  subplot(cfg_plot.numRows,cfg_plot.numCols,length(exper.subjects) + 1);
  ycoord = 1.0;
  if ~isempty(cfg_plot.types{typ})
    text(0.5,ycoord,sprintf('%s',cfg_plot.types{typ}),'color','k');
    ycoord = ycoord - vertTextLoc;
  end
  text(0.5,ycoord,sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}),'color','k');
  ycoord = ycoord - vertTextLoc;
  for evVal = 1:length(cfg_plot.conditions{ses}{typ})
    text(0.5,ycoord,strrep(cfg_plot.conditions{ses}{typ}{evVal},'_','-'),'color',cfg_plot.graphcolor(evVal));
    ycoord = ycoord - vertTextLoc;
  end
  axis off
end

end
