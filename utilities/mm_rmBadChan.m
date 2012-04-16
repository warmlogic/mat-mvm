function [data,ana] = mm_rmBadChan(cfg,exper,ana,data)
%MM_RMBADCHAN - Remove bad channels from individual subject data
%
% [data,ana] = mm_rmBadChan(cfg,exper,ana,data)
%
% Input:
%   cfg     = configuration struct
%   exper   = used for badChan field
%   ana     = used for eventValues and ana.elec.label fields
%   data    = ERP or TF data struct
%
% Output
%   data    = data structs with bad channels removed (deleted)
%
% NB: Currently this runs FieldTrip's ft_selectdata_old because there is a
%     bug in ft_selectdata_new (dimord field gets removed)
%     http://bugzilla.fcdonders.nl/show_bug.cgi?id=1409
%
%
% See also: MM_GETBADCHAN
%

if ~isfield(exper,'badChan')
  error('No bad channel information in exper struct. Run mm_getBadChan.');
end

if ~isfield(cfg,'rmBadChan')
  cfg.rmBadChan = 'remove';
end

if ~strcmp(cfg.rmBadChan,'no') && ~strcmp(cfg.rmBadChan,'remove') && ~strcmp(cfg.rmBadChan,'nan')
  error('cfg.rmBadChan must be set to ''no'', ''remove'', or ''nan''.');
end

if ~isfield(cfg,'printRoi')
  cfg.printRoi = [];
end

if ~isfield(cfg,'badChanThresh')
  cfg.badChanThresh = 0.2;
end

if ~isfield(cfg,'printBadChan')
  cfg.printBadChan = true;
end

% find the channels that are common across subjects; start with full set
cfg_com = [];
cfg_com.channel = ana.elec.label;

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    
    % this subject's channels
    subChan = eval(sprintf('{''all'' ''-Fid*''%s};',sprintf(repmat(' ''-E%d''',1,length(exper.badChan{sub,ses})),exper.badChan{sub,ses})));
    
    % whittle down to common channels
    cfg_com.channel = ft_channelselection(subChan,cfg_com.channel);
    
    cfg_sd = [];
    cfg_sd.channel = ft_channelselection(subChan,data.(ana.eventValues{1}{1}).sub(sub).ses(ses).data.label);
    
    if ~isempty(exper.badChan{sub,ses}) && sum(~ismember(data.(ana.eventValues{1}{1}).sub(sub).ses(ses).data.label,cfg_sd.channel)) > 0
      for typ = 1:length(ana.eventValues)
        for evVal = 1:length(ana.eventValues{typ})
          if isfield(data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'label')
            
            if strcmp(cfg.rmBadChan,'remove')
              % remove the bad channels
              fprintf('%s, %s, %s: Removing %d bad channel(s).\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal},length(exper.badChan{sub,ses}));
              % % ft_selectdata_new
              % data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_selectdata(cfg_sd,data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
              
              % ft_selectdata_old
              data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_selectdata(data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'channel',cfg_sd.channel);
            elseif strcmp(cfg.rmBadChan,'nan')
              % setting bad channels to nan - doesn't work with GA
              fprintf('%s, %s, %s: Setting %s bad channel(s) to NaNs.\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal},length(exper.badChan{sub,ses}));
              data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(exper.badChan{sub,ses},:) = NaN;
            end
            
          else
            warning([mfilename,':noData'],'%s %s %s: No data in struct.\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
          end
        end % evVal
      end % typ
    elseif ~isempty(exper.badChan{sub,ses}) && sum(~ismember(data.(ana.eventValues{1}{1}).sub(sub).ses(ses).data.label,cfg_sd.channel)) == 0
      fprintf('%s %s: Bad channels already removed.\n',exper.subjects{sub},exper.sesStr{ses});
    elseif isempty(exper.badChan{sub,ses})
      fprintf('%s %s: No bad channels.\n',exper.subjects{sub},exper.sesStr{ses});
    end % if
    continue
  end % ses
end % sub

if cfg.printBadChan
  fprintf('\n');
  % first list all the channels
  chansel = 'all';
  channels = ft_channelselection(eval(sprintf('{''%s'', ''-Fid*''}',chansel)),ana.elec.label);
  roi_str = sprintf('\t%s (%d)',chansel,length(channels));
  % then each roi
  for r = 1:length(cfg.printRoi)
    if ~strcmp(cfg.printRoi{r},'all')
      chansel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.printRoi{r})});
      channels = ft_channelselection(chansel,ana.elec.label);
      
      if length(cfg.printRoi{r}) > 1
        roi_substr = cfg.printRoi{r}{1};
        for i = 2:length(cfg.printRoi{r})
          roi_substr = sprintf('%s+%s',roi_substr,cfg.printRoi{r}{i});
        end
      else
        roi_substr = cfg.printRoi{r}{1};
      end
      
      roi_str = sprintf('%s\t%s (%d)',roi_str,roi_substr,length(channels));
      %chansel = ismember(ana.elec.label,channels);
    end
  end
  
  % find common channels across all subjects
  %channels = ft_channelselection(cfg_com.channel,ana.elec.label);
  if 1 - (length(cfg_com.channel) / length(ana.elec.label)) >= cfg.badChanThresh
    low_str = '<';
  else
    low_str = '';
  end
  common_str = sprintf('\t%d%s',length(cfg_com.channel),low_str);
  roi_split = regexp(roi_str,'\t','split');
  for r = 1:length(cfg.printRoi)
    if length(roi_split{r+1}) > 7
      tabchar = '\t';
    else
      tabchar = '';
    end
    if ~strcmp(cfg.printRoi{r},'all')
      chansel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.printRoi{r})});
      channels = ft_channelselection(chansel,cfg_com.channel);
      if 1 - (length(channels) / length(chansel)) >= cfg.badChanThresh
        low_str = '<';
      else
        low_str = '';
      end
      common_str = sprintf('%s%s\t%d',common_str,sprintf(tabchar),length(channels),low_str);
    end
  end
  
  if length(exper.subjects{1}) > 7
    tabchar_sub = '\t';
  else
    tabchar_sub = '';
  end
  
  fprintf('ROI%s%s\n',sprintf(tabchar_sub),roi_str);
  fprintf('Common%s%s\n',sprintf(tabchar_sub),common_str);
  %fprintf('\n');
  
  % print out the channel counts for each subject
  %fprintf('Subject%s%s\n',sprintf(tabchar_sub),sprintf(repmat('\t%s',1,length(eventValues)),eventValues{:}));
  for ses = 1:length(exper.sessions)
    for sub = 1:length(exper.subjects)
      % this subject's channels
      subChan = eval(sprintf('{''all'' ''-Fid*''%s};',sprintf(repmat(' ''-E%d''',1,length(exper.badChan{sub,ses})),exper.badChan{sub,ses})));
      %subChan = data.(ana.eventValues{1}{1}).sub(sub).ses(ses).data.label;
      subChan = ft_channelselection(subChan,ana.elec.label);
      sub_str = sprintf('%s%s%d',exper.subjects{sub},sprintf(tabchar_sub),length(subChan));
      for r = 1:length(cfg.printRoi)
        if length(roi_split{r+1}) > 7
          tabchar = '\t';
        else
          tabchar = '';
        end
        if ~strcmp(cfg.printRoi{r},'all')
          chansel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.printRoi{r})});
          channels = ft_channelselection(chansel,subChan);
          
          if 1 - (length(channels) / length(chansel)) >= cfg.badChanThresh
            low_str = '<';
          else
            low_str = '';
          end
          
          sub_str = sprintf('%s%s\t%d%s',sub_str,sprintf(tabchar),length(channels),low_str);
        end
      end
      fprintf('%s\n',sub_str);
    end
  end
  fprintf('< inducates # of remaining channels is <=%.1f%% of full ROI.\n',cfg.badChanThresh*100);
end

end
