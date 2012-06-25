function mm_printDataToText(cfg,exper,ana,dirs,data)
%mm_printDataToText - print data (e.g., voltages) to a file
%
% See also:
%   MM_FT_CHECKCONDITIONS

if ~isfield(cfg,'parameter')
  error('Must specify cfg.parameter, denoting the data to test (e.g., ''avg'' or ''individual'')');
end

if ~isfield(cfg,'direction')
  error('Must specify cfg.direction, denoting whether the data should be in ''columns'' or ''rows''.');
end

allROIs = cellflat(cfg.rois);
allROIs_str = sprintf(repmat('_%s',1,length(allROIs)),allROIs{:});
allConds = unique(cellflat(cfg.condByROI));
allConds_str = sprintf(repmat('_%s',1,length(allConds)),allConds{:});

if ~isfield(cfg,'excludeBadSub')
  cfg.excludeBadSub = 1;
elseif isfield(cfg,'excludeBadSub') && cfg.excludeBadSub ~= 1
  fprintf('Must exclude bad subjects. Setting cfg.excludeBadSub = 1;');
  cfg.excludeBadSub = 1;
end

% get the label info for this data struct
if isfield(data.(allConds{1}).sub(1).ses(1).data,'label');
  lab = data.(allConds{1}).sub(1).ses(1).data.label;
else
  error('label information not found in data struct');
end

% find out which were the good subjects
cfg.goodSub = exper.subjects(~exper.badSub);

if strcmp(cfg.direction,'columns')
  f_suf = '_col';
elseif strcmp(cfg.direction,'rows')
  f_suf = '_row';
end

outfile = fullfile(dirs.saveDirProc,sprintf('%s%s%s%s.txt',exper.name,allConds_str,allROIs_str,f_suf));
% open the file
fid = fopen(outfile,'w+');

if strcmp(cfg.direction,'columns')
  % print the main header
  fprintf(fid,'Experiment\tSubject\tSubject Index\tROI_Cond_Latency\n');
  
  fprintf(fid,'\t\t');
  for r = 1:length(cfg.rois)
    cfg.roi = cfg.rois{r};
    
    % set the channel information
    if ismember(cfg.roi,ana.elecGroupsStr)
      % set the string
      cfg.chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
    else
      % otherwise it should be the channel number(s) or 'all'
      if ~iscell(cfg.roi)
        cfg.roi = {cfg.roi};
      end
      % set the string
      cfg.chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
    end
    
    cfg.conditions = cfg.condByROI{r};
    cfg.latency = cfg.latencies(r,:);
    for evVal = 1:length(cfg.conditions)
      fprintf(fid,'\t%s%s_%d_%d',cfg.chan_str,cfg.conditions{evVal},cfg.latency(1)*1000,cfg.latency(2)*1000);
    end
  end
  fprintf(fid,'\n');
  
  ses = 1;
  goodSubjInd = 0;
  for sub = 1:length(exper.subjects)
    if ismember(exper.subjects{sub},cfg.goodSub)
      goodSubjInd = goodSubjInd + 1;
      % print the experiment name
      fprintf(fid,'%s\t',exper.name);
      fprintf(fid,'%s\t',exper.subjects{sub});
      fprintf(fid,'%d',goodSubjInd);
      
      for r = 1:length(cfg.rois)
        cfg.roi = cfg.rois{r};
        
        % set the channel information
        if ismember(cfg.roi,ana.elecGroupsStr)
          % if it's in the predefined ROIs, get the channel numbers
          cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
          % find the channel indices for averaging
          cfg.chansel = ismember(lab,cfg.channel);
          % % set the string for the filename
          % cfg.chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
        else
          % otherwise it should be the channel number(s) or 'all'
          if ~iscell(cfg.roi)
            cfg.roi = {cfg.roi};
          end
          cfg.channel = cfg.roi;
          
          % find the channel indices for averaging
          if strcmp(cfg.channel,'all')
            cfg.chansel = ismember(lab,ft_channelselection(cfg.channel,lab));
          else
            cfg.chansel = ismember(lab,cfg.channel);
          end
          % % set the string for the filename
          % cfg.chan_str = sprintf(repmat('%s_',1,length(cfg.channel)),cfg.channel{:});
        end
        
        cfg.conditions = cfg.condByROI{r};
        cfg.latency = cfg.latencies(r,:);
        
        for evVal = 1:length(cfg.conditions)
          ev = cfg.conditions{evVal};
          
          timesel = find(data.(ev).sub(sub).ses(ses).data.time >= cfg.latency(1) & data.(ev).sub(sub).ses(ses).data.time <= cfg.latency(2));
          voltage = mean(mean(data.(ev).sub(sub).ses(ses).data.(cfg.parameter)(cfg.chansel,timesel),1),2);
          fprintf(fid,'\t%.4f',voltage);
          
        end % c
      end % r
      fprintf(fid,'\n');
      
    else
      fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
      continue
    end % if good
  end % sub
  
elseif strcmp(cfg.direction,'rows')
  fprintf(fid,'Experiment\t%s\n',exper.name);
  
  % exclude the bad subjects from the subject count
  cfg.numSub = length(exper.subjects) - sum(exper.badSub);
  
  fprintf(fid,'Cond_ROI_Latency%s\n',sprintf(repmat('\t%s',1,length(cfg.goodSub)),cfg.goodSub{:}));
  
  for r = 1:length(cfg.rois)
    cfg.roi = cfg.rois{r};
    cfg.latency = cfg.latencies(r,:);
    cfg.conditions = cfg.condByROI{r};
    
    % set the channel information
    if ~isfield(cfg,'roi')
      error('Must specify either ROI names or channel names in cfg.roi');
    elseif isfield(cfg,'roi')
      if ismember(cfg.roi,ana.elecGroupsStr)
        % if it's in the predefined ROIs, get the channel numbers
        cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
        % find the channel indices for averaging
        cfg.chansel = ismember(lab,cfg.channel);
        % set the string for the filename
        cfg.chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
      else
        % otherwise it should be the channel number(s) or 'all'
        if ~iscell(cfg.roi)
          cfg.roi = {cfg.roi};
        end
        cfg.channel = cfg.roi;
        
        % find the channel indices for averaging
        if strcmp(cfg.channel,'all')
          cfg.chansel = ismember(lab,ft_channelselection(cfg.channel,lab));
        else
          cfg.chansel = ismember(lab,cfg.channel);
        end
        % set the string for the filename
        cfg.chan_str = sprintf(repmat('%s_',1,length(cfg.channel)),cfg.channel{:});
      end
    end
    
    
    % get times, data, SEM
    for evVal = 1:length(cfg.conditions)
      ev = cfg.conditions{evVal};
      cfg.values.(ev) = nan(cfg.numSub,length(exper.sessions));
      goodSubInd = 0;
      for sub = 1:length(exper.subjects)
        ses = 1;
        %for ses = 1:length(exper.sessions)
        if exper.badSub(sub,ses)
          fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
          continue
        else
          goodSubInd = goodSubInd + 1;
          cfg.timesel.(ev) = find(data.(ev).sub(sub).ses(ses).data.time >= cfg.latency(1) & data.(ev).sub(sub).ses(ses).data.time <= cfg.latency(2));
          cfg.values.(ev)(goodSubInd,ses) = mean(mean(data.(ev).sub(sub).ses(ses).data.(cfg.parameter)(cfg.chansel,cfg.timesel.(ev)),1),2);
        end
        %end % ses
      end % sub
      %cfg.sem.(ev) = std(cfg.values.(ev))/sqrt(length(cfg.values.(ev)));
      
      fprintf(fid,'%s%s_%d_%d%s\n',cfg.chan_str,ev,cfg.latency(1)*1000,cfg.latency(2)*1000,sprintf(repmat('\t%.4f',1,length(cfg.values.(ev))),cfg.values.(ev)));
      
    end % evVal
    
  end % r
end % direction

fprintf('Saving %s\n',fullfile(dirs.saveDirProc,sprintf('%s%s%s%s.txt',exper.name,allConds_str,allROIs_str,f_suf)));
fclose(fid);


end
