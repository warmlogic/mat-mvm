function [data,badChan,badEv] = mm_ft_artifact(dataroot,subject,sesName,eventValue,ana,exper,elecfile,data)
%MM_FT_ARTIFACT reject artifacts
% [data,badChan] = mm_ft_artifact(dataroot,subject,sesName,eventValue,ana,exper,elecfile,data)
%
% ana.artifact.type details are described in: SEG2FT, CREATE_FT_STRUCT
%
% See also: SEG2FT, CREATE_FT_STRUCT

% % HCGSN 129 Eye channels
% EOGV_upper = [25 8]; % left, right
% EOGV_lower = [127 126]; % left, right
% EOGH_left = 128;
% EOGH_right = 125;
% eog = {[25 127], [8 126]}; % left, right

badEv = [];

%% set the artifact processing parameters

% figure out which artifact options we're using
if ismember('nsAuto',ana.artifact.type)
  rejArt_nsAuto = true;
else
  rejArt_nsAuto = false;
end
if ismember('zeroVar',ana.artifact.type)
  rejArt_zeroVar = true;
else
  rejArt_zeroVar = false;
end
if ismember('badChanManual',ana.artifact.type)
  rejArt_badChanManual = true;
else
  rejArt_badChanManual = false;
end
if ismember('badChanEP',ana.artifact.type)
  rejArt_badChanEP = true;
else
  rejArt_badChanEP = false;
end
if ismember('rmBadChan',ana.artifact.type)
  rejArt_rmBadChan = true;
else
  rejArt_rmBadChan = false;
end
if ismember('preRejManual',ana.artifact.type)
  rejArt_preRejManual = true;
else
  rejArt_preRejManual = false;
end
% if ismember('trialNum',ana.artifact.type)
%   if ~isfield(ana.artifact,'trialNum')
%     error('Must define which trials to reject due to artifacts.');
%   elseif isfield(ana.artifact,'trialNum') && ~isfield(ana.artifact,eventValue,cell2mat(eventValue))
%     error('Must define which %s trials to reject due to artifacts.',cell2mat(eventValue));
%   end
%   rejArt_trialNum = true;
% else
%   rejArt_trialNum = false;
% end
if ismember('ftAuto',ana.artifact.type)
  rejArt_ftAuto = true;
else
  rejArt_ftAuto = false;
end
if ismember('ftManual',ana.artifact.type)
  rejArt_ftManual = true;
else
  rejArt_ftManual = false;
end
if ismember('ftICA',ana.artifact.type)
  rejArt_ftICA = true;
  % manual artifact rejection is required before ICA
  if ~rejArt_ftManual
    rejArt_ftManual = true;
  end
else
  rejArt_ftICA = false;
end

if rejArt_ftAuto || rejArt_ftManual || rejArt_ftICA
  ft_artTypes = {'clip','ecg','eog','jump','manual','muscle','threshold','tms','zvalue','visual'};
end

% if rejArt_nsAuto && rejArt_zeroVar
%   error('Cannot reject both NS artifacts (''nsAuto'') and trials with zero variance (''zeroVar''). Choose one or the other.')
% end

if rejArt_preRejManual && ~rejArt_nsAuto && ~rejArt_zeroVar
  error('To manually inspect prerejected artifacts (''preRejManual''), you must also use either ''nsAuto'' or ''zeroVar''. Otherwise, just use ''ftManual'')');
end

%% check on predefined trial numbers, NS, and zero variance artifacts;
% manual inspection option will show which have been rejected

foundArt = false;

% % reject using trial numbers
% if rejArt_trialNum
%   if ~isempty(find(ana.artifact.trialNum.(cell2mat(eventValue)),1))
%     foundArt = true;
%     fprintf('Automatically rejecting %d predefined artifacts for%s.\n',length(ana.artifact.trialNum.(cell2mat(eventValue))),sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
%     cfg = [];
%     % mark the trials that have artifacts as such; select the entire sample
%     % range for the bad events
%     cfg.artfctdef.visual.artifact = data.sampleinfo(ana.artifact.trialNum.(eventValue),:);
%   else
%     fprintf('No predefined artifacts found for%s.\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
%   end
% end

if rejArt_badChanManual
  foundBadChan = false;
  fid = fopen(fullfile(dataroot,[exper.name,'_badChan.txt']));
  
  if fid == -1
    error('Could not open %s. Make sure you do not have it open in another application.',fullfile(dataroot,[exper.name,'_badChan.txt']));
  else
    [badChanInfo] = textscan(fid,'%s%s%s','delimiter','\t');
    fclose(fid);
    for i = 1:length(badChanInfo{1})
      if strcmp(badChanInfo{1}{i},subject) && strcmp(badChanInfo{2}{i},sesName)
        badChanInfo{3}{i} = strrep(strrep(badChanInfo{3}{i},'[',''),']','');
        badChanManual = eval(sprintf('[%s];',badChanInfo{3}{i}));
        foundBadChan = true;
        break
      else
        badChanManual = [];
      end
    end
    if ~foundBadChan
      error('No manual bad channel information found for %s %s.',subject,sesName);
    end
  end
end

if rejArt_badChanEP
  foundBadChan = false;
  foundNoBadChan = false;
  
  epGBC_str = 'Global bad Channels: ';
  ep_prefix = 'Artifact_Correction_Log';
  
  %   datasetFile = ft_findcfg(data.cfg,'dataset');
  %   if ~isempty(datasetFile)
  %     datasetSep = strfind(datasetFile,filesep);
  %     if ~isempty(datasetSep)
  %       datasetFile = datasetFile(datasetSep(end)+1:end);
  %     end
  %     datasetSep = strfind(datasetFile,'.');
  %     datasetFile = datasetFile(1:datasetSep(end)-1);
  %     datasetFile = strrep(datasetFile,'_e','');
  %   else
  %     datasetFile = subject;
  %   end
  %   epArtFile = dir(fullfile(dataroot,'ep_art',sesName,sprintf('%s %s*.txt',ep_prefix,datasetFile)));
  epArtFile = dir(fullfile(dataroot,'ep_art',sesName,sprintf('%s %s*.txt',ep_prefix,subject)));
  
  if ~isempty(epArtFile)
    if length(epArtFile) == 1
      epArtFile = epArtFile.name;
      fid = fopen(fullfile(dataroot,'ep_art',sesName,epArtFile),'r');
      if fid == -1
        error('Could not open the file for %s %s. Make sure you do not have it open in another application.',subject,sesName);
      else
        
        while ~foundBadChan
          % get the next line
          tline = fgetl(fid);
          if strncmpi(tline,epGBC_str,length(epGBC_str))
            foundBadChan = true;
            break
          end
        end
        fclose(fid);
        
        badChanInfo = strrep(tline,epGBC_str,'');
        if strcmp(badChanInfo(end),sprintf('\t'))
          badChanInfo = badChanInfo(1:end-1);
        end
        if strcmp(badChanInfo,'None')
          foundNoBadChan = true;
        end
        if foundBadChan
          if ~foundNoBadChan
            %badChanEP = regexp(badChanEP,'\t','split');
            badChanEP = eval(sprintf('[%s];',badChanInfo));
          elseif foundNoBadChan
            badChanEP = [];
          end
        elseif ~foundBadChan
          error('No bad channel information found.');
        end
      end
    elseif length(epArtFile) > 1
      error('More than one EP Articifact Correction Log found for %s %s, probably multiple sessions in the same directory. Make sure each subject has one Artifact Log file in individual session folders, or just put the info in a text file.',subject,sesName);
    end
  elseif isempty(epArtFile)
    error('No EP Articifact Correction Log found for %s %s.',subject,sesName);
  end
end

% collect the bad channels
if rejArt_badChanManual || rejArt_badChanEP
  if ~exist('badChanManual','var')
    badChanManual = [];
  end
  if ~exist('badChanEP','var')
    badChanEP = [];
  end
  badChan = unique(cat(2,badChanManual,badChanEP));
else
  badChan = [];
end

if rejArt_rmBadChan && ~isempty(badChan)
  cfg_rv = [];
  cfg_rv.channel = eval(sprintf('{''all''%s};',sprintf(repmat(' ''-E%d''',1,length(badChan)),badChan)));
  cfg_rv.keepchannel = 'no';
  %cfg_rv.keepchannel = 'nan';
  
  data = ft_rejectvisual(cfg_rv,data);
elseif rejArt_rmBadChan && isempty(badChan)
  fprintf('No bad channels to reject!\n');
end

if rejArt_nsAuto
  % make sure the file with NS artifact info exists
  summaryFile = dir(fullfile(dataroot,sesName,'ns_bci',[subject,'*.bci']));
  if ~isempty(summaryFile)
    summaryFile = fullfile(dataroot,sesName,'ns_bci',summaryFile.name);
    
    fid = fopen(summaryFile,'r');
    % get the header
    fgetl(fid);
    % get the next line
    tline = fgetl(fid);
    fclose(fid);
    % get only the second half of the line
    tline = tline(ceil(length(tline)/2):end);
    % figure out how many channels are in the summary file
    if ~isempty(strfind(tline,'129')) && length(strfind(tline,'129')) == 1
      nChan_summary = 129;
    elseif ~isempty(strfind(tline,'128')) && length(strfind(tline,'128')) == 1
      nChan_summary = 128;
    elseif ~isempty(strfind(tline,'257')) && length(strfind(tline,'257')) == 1
      nChan_summary = 257;
    elseif ~isempty(strfind(tline,'256')) && length(strfind(tline,'256')) == 1
      nChan_summary = 256;
    else
      nChan_summary = nChan_elecfile;
    end
    
    % read in the NS artifact file
    format_str = ['%s%d%d%s',repmat('%d',[1,nChan_summary*2]),'%s'];
    fid = fopen(summaryFile,'r');
    sesSummary = textscan(fid,format_str,'Headerlines',1,'delimiter','\t');
    fclose(fid);
  else
    error('Cannot find %s*.%s file in %s. Use the File Export tool to export Metadata > Segment Information.',subject,'bci',fullfile(dataroot,'ns_bci'));
  end
  
  % find the bad trials for this event
  %thisEv = find(ismember(sesSummary{1},eventValue));
  badEv = strcmp(sesSummary{4},'bad');
  % make sure we only have data for the bad events for this event value
  badEv = badEv(ismember(sesSummary{1},eventValue));
  %goodEv = strcmp(sesSummary{4},'good');
  %cfg.trials = logical(goodEv(min(thisEv):max(thisEv)));
  
  % trials in raw files aren't stored in the same order as they are
  % listed in the bci file. raw files seem to have trials in the order
  % they occurred in the experiment while bci files have trials in the
  % order that they're listed in the segmentation tool, then sorted by
  % time. the solution is to sort both, then put the bci good/bad status
  % in the same order as the raw file. this shouldn't affect files like
  % egis where the order is already the same.
  
  % put the trialinfo events in sort order
  [tiEvSort,tiEvInd] = sort(eventValue(data.trialinfo));
  
  % put the bci events and statuses in sort order
  [bciEvSort,bciEvInd] = sort(sesSummary{1}(ismember(sesSummary{1},eventValue)));
  badEvSort = badEv(bciEvInd);
  
  % set indices for trialinfo events
  tiOrigInd(tiEvInd) = 1:length(tiEvInd);
  
  % put badEv in order of trialinfo events
  badEv = badEvSort(tiOrigInd);
  
  % % testing
  % bciEvSort(tiOrigInd)
  
  % old method: remove trials with artifacts from the trl matrix
  %cfg.trl(logical(badEv(min(thisEv):max(thisEv))),:) = [];
  
  %   % another old method
  %   %
  %   % I don't remember why I implemented the sort(), but it doesn't work
  %   % right for me so I'm going back to the old method. It doesn't work
  %   % because badEv is now in a different order from data.trialinfo and
  %   % data.sampleinfo.
  %
  %   % sort the events by StartTime offset in the bci file
  %   [startTimeSorted,startTimeInd] = sort(sesSummary{2});
  %   badSorted = sesSummary{4}(startTimeInd);
  %   eventValueSorted = sesSummary{1}(startTimeInd);
  %   % find the bad trials for this event
  %   badEv = strcmp(badSorted,'bad');
  %   % make sure we only have data for the bad events for this event value
  %   badEv = badEv(ismember(eventValueSorted,eventValue));

  if ~isempty(find(badEv,1))
    foundArt = true;
    fprintf('Automatically rejecting %d NS artifacts for%s.\n',sum(badEv),sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
    cfg = [];
    % mark the trials that have artifacts as such; select the entire sample
    % range for the bad events
    cfg.artfctdef.visual.artifact = data.sampleinfo(badEv,:);
    
    % this doesn't work when we're passing in multiple event values
    %cfg.artfctdef.visual.artifact = data.sampleinfo(logical(badEv(min(thisEv):max(thisEv))),:);
  else
    fprintf('No NS artifacts found for%s.\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
  end
end

% if we want to reject the trials with zero variance in their voltage
if rejArt_zeroVar
  badEv = zeros(size(data.trial'));
  for i = 1:size(data.trial,2)
    % variance across time on all the channels averaged together
    if var(mean(data.trial{i},1),0,2) == 0
      badEv(i) = 1;
    end
  end
  
  if ~isempty(find(badEv,1))
    foundArt = true;
    fprintf('Automatically rejecting %d trials with zero variance for%s.\n',sum(badEv),sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
    
    if ~exist('cfg','var')
      cfg = [];
      cfg.artfctdef.visual.artifact = data.sampleinfo(logical(badEv),:);
    else
      if rejArt_nsAuto && foundArt
        % if running both nsAuto and zeroVar
        if isfield(cfg,'artfctdef')
          if isfield(cfg.artfctdef,'visual')
            if isfield(cfg.artfctdef.visual,'artifact')
              if ~isempty(cfg.artfctdef.visual.artifact)
                cfg.artfctdef.visual.artifact = unique(cat(1,cfg.artfctdef.visual.artifact,data.sampleinfo(logical(badEv),:)),'rows');
              end
            end
          end
        end
      end
    end
  else
    fprintf('No zero variance trials found for%s.\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
  end
end

% if we want to inspect manually
if (rejArt_nsAuto || rejArt_zeroVar) && rejArt_preRejManual
  foundArt = true;
  %cfg.viewmode = 'butterfly';
  cfg.viewmode = 'vertical';
  cfg.continuous = 'no';
  cfg.elecfile = elecfile;
  cfg.plotlabels = 'some';
  cfg.ylim = [-10 10];
  
  fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
  fprintf('\n\nManual artifact rejection (NS artifacts are marked):\n');
  fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
  fprintf('\tUse arrows to move to next trial.\n');
  fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
  fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
  fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n\n');
  
  cfg = ft_databrowser(cfg,data);
  % bug when calling rejectartifact right after databrowser, pause first
  pause(1);
  
  keepRepairingChannels = true;
  while keepRepairingChannels
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\n\nDo you want to see whether there are channels to repair? (0 or 1, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      cfgChannelRepair.viewmode = 'butterfly';
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\n\nDo you want to plot a particular subset of channels (1) or all channels (0)? (0 or 1, then press ''return''):\n\n');
      end
      if repair_chanNum
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
      else
        cfgChannelRepair.channel = 'all';
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('\tUse the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\n\nWere there channels to repair? (0 or 1, then press ''return''):\n\n');
      end
      
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = 'all';
        cfgChannelRepair.badchannel = badchannel;
        cfgChannelRepair.method = 'spline';
        cfgChannelRepair.elecfile = elecfile;
        fprintf('Repairing channels%s using method=''%s''...\n',sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}),cfgChannelRepair.method);
        data = ft_channelrepair(cfgChannelRepair, data);
        
        %     channelsToRepair = [];
        %     while ~iscell(channelsToRepair)
        %       channelsToRepair = input('\nType channel labels to repair (on a single line) and press ''return'' (cell array of strings, e.g., {''E1'',''E4'',''E11''}). If no, type {}.\n\n');
        %     end
        %
        %     if ~isempty(channelsToRepair)
        %       data.elec = elec;
        %
        %       cfg_repair = [];
        %       cfg_repair.badchannel = channelsToRepair;
        %
        %       data = ft_channelrepair(cfg_repair,data);
        %     end
      end
    else
      keepRepairingChannels = false;
    end
  end
end

% do the actual rejection of artifact trials (complete or parial rejection)
if (rejArt_nsAuto || rejArt_zeroVar) && foundArt
  cfg.artfctdef.reject = 'complete';
  data = ft_rejectartifact(cfg,data);
end

%% visual artifact inspection (manual)

if rejArt_ftManual
  keepRepairingChannels = true;
  while keepRepairingChannels
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\n\nDo you want to see whether there are any really bad channels to repair? (0 or 1, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      cfgChannelRepair.viewmode = 'butterfly';
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\n\nDo you want to plot a particular subset of channels (1) or all channels (0)? (0 or 1, then press ''return''):\n\n');
      end
      if repair_chanNum
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
      else
        cfgChannelRepair.channel = 'all';
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('Use the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('Use the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\n\nWere there channels to repair? (0 or 1, then press ''return''):\n\n');
      end
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = 'all';
        cfgChannelRepair.badchannel = badchannel;
        cfgChannelRepair.method = 'spline';
        cfgChannelRepair.elecfile = elecfile;
        fprintf('Repairing channels%s using method=''%s''...\n',sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}),cfgChannelRepair.method);
        data = ft_channelrepair(cfgChannelRepair, data);
      end
    else
      keepRepairingChannels = false;
    end
  end
  
  if rejArt_ftICA
    fprintf('Before running ICA, you must manually reject artifacts that are not consistent across trials. Do not reject blinks! Please do this now.\n');
  end
  
  % use cursor drag and click to mark artifacts;
  % use arrows to advance to next trial;
  % use the q key to quit the data browser
  
  cfg = [];
  cfg.continuous = 'no';
  %cfg.padding = 0;
  % get the trial definition for automated FT artifact rejection
  cfg.trl = ft_findcfg(data.cfg,'trl');
  
  cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.cutoff = 20;
  cfg.artfctdef.zvalue.trlpadding = 0;
  cfg.artfctdef.zvalue.artpadding = 0.1;
  cfg.artfctdef.zvalue.fltpadding = 0;
  
  % algorithmic parameters
  cfg.artfctdef.zvalue.absdiff = 'yes';
  
  %cfg.artfctdef.zvalue.interactive = 'yes';
  
  fprintf('Checking for zvalue artifacts at z=%d...\n',cfg.artfctdef.zvalue.cutoff);
  
  % auto mark some artifacts
  cfg = ft_artifact_zvalue(cfg, data);
  
  % search for muscle artifacts
  cfg.continuous = 'no';
  %cfg.padding = 0;
  % get the trial definition for automated FT artifact rejection
  cfg.trl = ft_findcfg(data.cfg,'trl');
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection
  cfg.artfctdef.muscle.channel = 'all';
  cfg.artfctdef.muscle.cutoff      = 40;
  cfg.artfctdef.muscle.trlpadding  = 0;
  if strcmp(cfg.continuous,'yes')
    cfg.artfctdef.muscle.artpadding  = 0.1;
  elseif strcmp(cfg.continuous,'no')
    cfg.artfctdef.muscle.artpadding  = 0.1;
  end
  cfg.artfctdef.muscle.fltpadding  = 0;
  
  % algorithmic parameters
  cfg.artfctdef.muscle.bpfilter    = 'yes';
  if data.fsample/2 < 140
    cfg.artfctdef.muscle.bpfreq      = [110 (data.fsample/2 - 1)];
  else
    cfg.artfctdef.muscle.bpfreq      = [110 140];
  end
  cfg.artfctdef.muscle.bpfiltord   = 6;
  cfg.artfctdef.muscle.bpfilttype  = 'but';
  cfg.artfctdef.muscle.hilbert     = 'yes';
  cfg.artfctdef.muscle.boxcar      = 0.2;
  
  %cfg.artfctdef.zvalue.interactive = 'yes';
  
  fprintf('Checking for muscle artifacts at z=%d...\n',cfg.artfctdef.muscle.cutoff);
  
  % auto mark some artifacts
  [cfg] = ft_artifact_muscle(cfg,data);
  
  cfg.continuous = 'no';
  %cfg.viewmode = 'butterfly';
  cfg.viewmode = 'vertical';
  cfg.elecfile = elecfile;
  cfg.plotlabels = 'some';
  cfg.ylim = [-10 10];
  
  % manual rejection
  fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
  fprintf('\n\nManual artifact rejection:\n');
  fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
  fprintf('\tUse arrows to move to next trial.\n');
  if strcmp(cfg.viewmode,'butterfly')
    fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
  end
  fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
  fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n\n');
  
  cfg = ft_databrowser(cfg,data);
  % bug when calling rejectartifact right after databrowser, pause first
  pause(1);
  
  % save a list of trials with artifact status
  %badEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
  badEv = zeros(size(data.sampleinfo,1), 1);
  for i = 1:length(ft_artTypes)
    if isfield(cfg.artfctdef,ft_artTypes{i})
      for j = 1:size(cfg.artfctdef.(ft_artTypes{i}).artifact,1)
        for k = 1:size(badEv,1)
          if cfg.artfctdef.(ft_artTypes{i}).artifact(j,1) >= data.sampleinfo(k,1) && cfg.artfctdef.(ft_artTypes{i}).artifact(j,2) <= data.sampleinfo(k,2)
            %badEv(k,2) = 1;
            badEv(k,1) = 1;
          end
        end
      end
    end
  end
  
  % reject the artifacts (complete or parial rejection)
  cfg.artfctdef.reject = 'complete';
  data = ft_rejectartifact(cfg,data);
  
  keepRepairingChannels = true;
  while keepRepairingChannels
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\n\nDo you want to see whether there are channels to repair? (0 or 1, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      cfgChannelRepair.viewmode = 'butterfly';
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\n\nDo you want to plot a particular subset of channels (1) or all channels (0)? (0 or 1, then press ''return''):\n\n');
      end
      if repair_chanNum
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
      else
        cfgChannelRepair.channel = 'all';
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('Use the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('Use the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\n\nWere there channels to repair? (0 or 1, then press ''return''):\n\n');
      end
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = 'all';
        cfgChannelRepair.badchannel = badchannel;
        cfgChannelRepair.method = 'spline';
        cfgChannelRepair.elecfile = elecfile;
        fprintf('Repairing channels%s using method=''%s''...\n',sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}),cfgChannelRepair.method);
        data = ft_channelrepair(cfgChannelRepair, data);
      end
    else
      keepRepairingChannels = false;
    end
  end
end

%% ICA artifact detection

if rejArt_ftICA
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ICA artifact detection
  %
  % look for eye blink, heart beat, and other artifacts
  %
  % Questionable claim: you should not run ICA on data that has already
  % rejected ICA components
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.channel = 'all';
  cfg.method = 'runica';
  
  comp = ft_componentanalysis(cfg,data);
  
  % % OLD METHOD - view the first 30 components
  % cfg = [];
  % cfg.component = 1:30; % specify the component(s) that should be plotted
  % cfg.layout = elecfile; % specify the layout file that should be used for plotting
  % cfg.comment = 'no';
  % ft_topoplotIC(cfg,comp);
  
  % view some components (aka channels) with time course
  %nComponents = 10;
  cfg = [];
  cfg.viewmode = 'component';
  cfg.continuous = 'yes';
  % number of seconds to display
  cfg.blocksize = 30;
  %cfg.blocksize = 10;
  %cfg.channels = 1:nComponents;
  cfg.plotlabels = 'yes';
  cfg.elecfile = elecfile;
  cfg.ylim = [-10 10];
  ft_databrowser(cfg,comp);
  % bug when calling rejectartifact right after databrowser, pause first
  pause(1);
  
  fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
  %fprintf('\n\nViewing the first %d components.\n',nComponents);
  fprintf('ICA component browsing:\n');
  fprintf('\t1. Look for patterns that are indicative of artifacts.\n');
  fprintf('\t\tPress the ''channel >'' button to see the next set of components.\n');
  fprintf('\t\tComponents may not be numbered, so keep track of where you are (top component has the lowest number). Note component numbers for rejection.\n');
  fprintf('\t2. Manually close the components window when finished browsing.\n');
  
  rej_comp = [];
  while isempty(rej_comp) || (rej_comp ~= 0 && rej_comp ~= 1)
    rej_comp = input('Were there components to reject? (0 or 1, then press ''return''):\n\n');
  end
  if rej_comp
    % prompt the user for the component numbers to reject
    componentsToReject = input(sprintf('\t3. Type component numbers to reject (on a single line) and press ''return'', even if these instructions move up due to output while browsing components (e.g., ''1, 4, 11'' without quotes):\n\n'),'s');
    
    % reject the bad components
    if ~isempty(componentsToReject)
      cfg = [];
      cfg.component = str2double(regexp(componentsToReject,'\d*','match')');
      data = ft_rejectcomponent(cfg, comp, data);
    end
  end
  
  % another auto search for artifacts
  
  % use cursor drag and click to mark artifacts;
  % use arrows to advance to next trial;
  % use the q key to quit the data browser
  
  cfg = [];
  cfg.continuous = 'no';
  %cfg.padding = 0;
  % get the trial definition for automated FT artifact rejection
  cfg.trl = ft_findcfg(data.cfg,'trl');
  
  cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.cutoff = 28;
  cfg.artfctdef.zvalue.trlpadding = 0;
  cfg.artfctdef.zvalue.artpadding = 0.1;
  cfg.artfctdef.zvalue.fltpadding = 0;
  
  %cfg.artfctdef.zvalue.interactive = 'yes';
  
  fprintf('Checking for zvalue artifacts at z=%d...\n',cfg.artfctdef.zvalue.cutoff);
  
  % auto mark some artifacts
  cfg = ft_artifact_zvalue(cfg, data);
  
  % another manual search of the data for artifacts
  
  %cfg.viewmode = 'butterfly';
  cfg.viewmode = 'vertical';
  cfg.continuous = 'no';
  cfg.elecfile = elecfile;
  cfg.plotlabels = 'some';
  cfg.ylim = [-10 10];
  
  fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
  fprintf('\n\nFinal round of manual artifact rejection:\n');
  fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
  fprintf('\tUse arrows to move to next trial.\n');
  if strcmp(cfg.viewmode,'butterfly')
    fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
  end
  fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
  fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n\n');
  
  cfg = ft_databrowser(cfg,data);
  % bug when calling rejectartifact right after databrowser, pause first
  pause(1);
  
  % save a list of trials with artifact status
  if ~exist('badEv','var') || isempty(badEv)
    combineArtLists = false;
    %badEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
    badEv = zeros(size(data.sampleinfo,1), 1);
  else
    combineArtLists = true;
  end
  %postIcaEv = badEv(badEv(:,2) == 0,:);
  %postIcaEv = badEv(badEv == 0);
  %postIcaEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
  postIcaEv = zeros(size(data.sampleinfo,1), 1);
  
  % is size(data.sampleinfo,1) == size(postIcaEv,1)
  % is size(badEv(badEv == 0),1) == size(postIcaEv,1)
  % is size(badEv(badEv == 0),1) == size(data.sampleinfo,1)
  %
  % these should be the same size
  if size(badEv(badEv == 0),1) ~= size(data.sampleinfo,1)
    keyboard
  end
  
  for i = 1:length(ft_artTypes)
    if isfield(cfg.artfctdef,ft_artTypes{i})
      for j = 1:size(cfg.artfctdef.(ft_artTypes{i}).artifact,1)
        for k = 1:size(postIcaEv,1)
          if cfg.artfctdef.(ft_artTypes{i}).artifact(j,1) >= data.sampleinfo(k,1) && cfg.artfctdef.(ft_artTypes{i}).artifact(j,2) <= data.sampleinfo(k,2)
            %postIcaEv(k,2) = 1;
            postIcaEv(k,1) = 1;
          end
        end
      end
    end
  end
  if combineArtLists
    % put the new artifacts into the old list
    rCount = 0;
    for i = 1:size(badEv,1)
      %if badEv(i,2) == 0
      if badEv(i) == 0
        rCount = rCount + 1;
        %if postIcaEv(rCount,2) == 1
        if postIcaEv(rCount) == 1
          %badEv(i,2) = 1;
          badEv(i) = 1;
        end
      end
    end
  else
    badEv = postIcaEv;
  end
  
  % reject the artifacts (complete or parial rejection)
  cfg.artfctdef.remove = 'complete';
  % and reject
  data = ft_rejectartifact(cfg,data);
  
  keepRepairingChannels = true;
  while keepRepairingChannels
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\n\nDo you want to see whether there are channels to repair? (0 or 1, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      cfgChannelRepair.viewmode = 'butterfly';
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\n\nDo you want to plot a particular subset of channels (1) or all channels (0)? (0 or 1, then press ''return''):\n\n');
      end
      if repair_chanNum
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
      else
        cfgChannelRepair.channel = 'all';
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('Use the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('Use the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\n\nWere there channels to repair? (0 or 1, then press ''return''):\n\n');
      end
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = 'all';
        cfgChannelRepair.badchannel = badchannel;
        cfgChannelRepair.method = 'spline';
        cfgChannelRepair.elecfile = elecfile;
        fprintf('Repairing channels%s using method=''%s''...\n',sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}),cfgChannelRepair.method);
        data = ft_channelrepair(cfgChannelRepair, data);
      end
    else
      keepRepairingChannels = false;
    end
  end
end

%% run FieldTrip's automatic artifact detection on the data

if rejArt_ftAuto
  % get the trial definition for automated FT artifact rejection
  trl = ft_findcfg(data.cfg,'trl');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for jump artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.padding = 0;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection
  cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.cutoff = 20;
  %cfg.artfctdef.zvalue.trlpadding = 0.5*cfg.padding;
  %cfg.artfctdef.zvalue.artpadding = 0.5*cfg.padding;
  cfg.artfctdef.zvalue.trlpadding = 0;
  cfg.artfctdef.zvalue.artpadding = 0.1;
  cfg.artfctdef.zvalue.fltpadding = 0;
  
  % algorithmic parameters
  cfg.artfctdef.zvalue.cumulative = 'yes';
  cfg.artfctdef.zvalue.medianfilter = 'yes';
  cfg.artfctdef.zvalue.medianfiltord = 9;
  cfg.artfctdef.zvalue.absdiff = 'yes';
  
  % feedback (artifact viewer)
  %cfg.artfctdef.zvalue.feedback = 'yes';
  cfg.artfctdef.zvalue.interactive = 'yes';
  
  [cfg,artifact_jump] = ft_artifact_zvalue(cfg,data);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for muscle artifacts - mess with bpfiltord
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.padding = 0;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection
  cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.cutoff      = 40;
  %cfg.artfctdef.zvalue.trlpadding  = 0.1*cfg.padding;
  %cfg.artfctdef.zvalue.fltpadding  = 0.1*cfg.padding;
  %cfg.artfctdef.zvalue.artpadding  = 0.1*cfg.padding;
  cfg.artfctdef.zvalue.trlpadding  = 0;
  if strcmp(cfg.continuous,'yes')
    cfg.artfctdef.zvalue.artpadding  = 0.1;
  elseif strcmp(cfg.continuous,'no')
    cfg.artfctdef.zvalue.artpadding  = 0.1;
  end
  cfg.artfctdef.zvalue.fltpadding  = 0;
  
  % algorithmic parameters
  cfg.artfctdef.zvalue.bpfilter    = 'yes';
  if data.fsample/2 < 140
    cfg.artfctdef.zvalue.bpfreq      = [110 (data.fsample/2 - 1)];
  else
    cfg.artfctdef.zvalue.bpfreq      = [110 140];
  end
  %cfg.artfctdef.zvalue.bpfiltord   = 9;
  cfg.artfctdef.zvalue.bpfiltord   = 6;
  cfg.artfctdef.zvalue.bpfilttype  = 'but';
  cfg.artfctdef.zvalue.hilbert     = 'yes';
  cfg.artfctdef.zvalue.boxcar      = 0.2;
  
  % feedback
  %cfg.artfctdef.zvalue.feedback = 'yes';
  cfg.artfctdef.zvalue.interactive = 'yes';
  
  [cfg,artifact_muscle] = ft_artifact_zvalue(cfg,data);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for EOG artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.padding = 0;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
  %cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.channel = {'E127','E126','E128','E125'};
  cfg.artfctdef.zvalue.cutoff      = 4;
  %cfg.artfctdef.zvalue.trlpadding  = 0.5*cfg.padding;
  %cfg.artfctdef.zvalue.artpadding  = 0.1*cfg.padding;
  %cfg.artfctdef.zvalue.fltpadding  = 0.1*cfg.padding;
  cfg.artfctdef.zvalue.trlpadding  = 0;
  if strcmp(cfg.continuous,'yes')
    cfg.artfctdef.zvalue.artpadding  = 0.1;
  elseif strcmp(cfg.continuous,'no')
    cfg.artfctdef.zvalue.artpadding  = 0.1;
  end
  cfg.artfctdef.zvalue.fltpadding  = 0;
  
  % algorithmic parameters
  cfg.artfctdef.zvalue.bpfilter   = 'yes';
  cfg.artfctdef.zvalue.bpfilttype = 'but';
  cfg.artfctdef.zvalue.bpfreq     = [1 15];
  cfg.artfctdef.zvalue.bpfiltord  = 4;
  cfg.artfctdef.zvalue.hilbert    = 'yes';
  
  % feedback
  %cfg.artfctdef.zvalue.feedback = 'yes';
  cfg.artfctdef.zvalue.interactive = 'yes';
  
  [cfg,artifact_EOG] = ft_artifact_zvalue(cfg,data);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % reject the automatically defined artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
  cfg.artfctdef.jump.artifact = artifact_jump;
  cfg.artfctdef.muscle.artifact = artifact_muscle;
  cfg.artfctdef.eog.artifact = artifact_EOG;
  data = ft_rejectartifact(cfg,data);
end

end
