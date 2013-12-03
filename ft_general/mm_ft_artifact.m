function [data,badChan_str,badEv] = mm_ft_artifact(dataroot,subject,sesName,eventValue,ana,exper,elecfile,data,dirs)
%MM_FT_ARTIFACT reject artifacts
% [data,badChan] = mm_ft_artifact(dataroot,subject,sesName,eventValue,ana,exper,elecfile,data,dirs)
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
badChan = [];
badChan_str = {};

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

% if rejArt_nsAuto && rejArt_zeroVar
%   error('Cannot reject both NS artifacts (''nsAuto'') and trials with zero variance (''zeroVar''). Choose one or the other.')
% end

if rejArt_preRejManual && ~rejArt_nsAuto && ~rejArt_zeroVar
  error('To manually inspect prerejected artifacts (''preRejManual''), you must also use either ''nsAuto'' or ''zeroVar''. Otherwise, just use ''ftManual'')');
end

vert_ylim = [-10 10];

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
  badChan_str = cell(length(badChan,1));
  for i = 1:length(badChan)
    badChan_str{i} = sprintf('E%d',badChan(i));
  end
% else
%   badChan = [];
%   badChan_str = {};
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
  cfg.ylim = vert_ylim;
  
  fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
  fprintf('\n\nManual artifact rejection (NS artifacts are marked):\n');
  fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
  fprintf('\tUse arrows to move to next trial.\n');
  fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
  fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
  fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n');
  
  cfg = ft_databrowser(cfg,data);
  % bug when calling rejectartifact right after databrowser, pause first
  pause(1);
  
  keepRepairingChannels = true;
  while keepRepairingChannels
    if ~exist('badChan_str','var')
      badChan_str = {};
    end
    
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\nDo you want to see whether there are channels to repair? (1 or 0, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      cfgChannelRepair.viewmode = 'butterfly';
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\nDo you want to plot all channels (1) or a particular subset of channels (0)? (1 or 0, then press ''return''):\n\n');
      end
      if repair_chanNum
        cfgChannelRepair.channel = 'all';
      else
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
        fprintf('(Manually close any empty figure windows.)\n');
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('\nUse the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('Use the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\nWere there channels to repair? (1 or 0, then press ''return''):\n\n');
      end
      
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        fprintf('(Manually close any empty figure windows.)\n');
        badChan_str = cat(1,badChan_str,badchannel);
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
  % set the file to save after checking artifacts, in case MATLAB crashes
  resumeManArtFT_file = fullfile(dirs.saveDirRaw,subject,sesName,sprintf('%s_%s_manArtFT.mat',subject,sesName));
  %resumeManArtFT_file = fullfile(dataroot,sprintf('%s_%s_manArtFT%s.mat',subject,sesName,sprintf(repmat('_%s',1,length(eventValue)),eventValue{:})));
  if ~isfield(ana.artifact,'resumeManArtFT')
    ana.artifact.resumeManArtFT = false;
  end
  
  % set some defaults
  if isfield(ana.artifact,'basic_art_z')
    basic_art_z_default = ana.artifact.basic_art_z;
  elseif ~isfield(ana.artifact,'basic_art_z')
    basic_art_z_default = 30;
    ana.artifact.basic_art_z = basic_art_z_default;
  end
  if isfield(ana.artifact,'muscle_art_z')
    muscle_art_z_default = ana.artifact.muscle_art_z;
  elseif ~isfield(ana.artifact,'muscle_art_z')
    muscle_art_z_default = 40;
    ana.artifact.muscle_art_z = muscle_art_z_default;
  end
  if isfield(ana.artifact,'jump_art_z')
    jump_art_z_default = ana.artifact.jump_art_z;
  elseif ~isfield(ana.artifact,'jump_art_z')
    jump_art_z_default = 50;
    ana.artifact.jump_art_z = jump_art_z_default;
  end
  
  if isfield(ana.artifact,'threshmin')
    threshmin_default = ana.artifact.threshmin;
  elseif ~isfield(ana.artifact,'threshmin')
    threshmin_default = -150;
    ana.artifact.threshmin = threshmin_default;
  end
  if isfield(ana.artifact,'threshmax')
    threshmax_default = ana.artifact.threshmax;
  elseif ~isfield(ana.artifact,'threshmax')
    threshmax_default = 150;
    ana.artifact.threshmax = threshmax_default;
  end
  
  if ana.artifact.resumeManArtFT
    if exist(resumeManArtFT_file,'file')
      % load the manually processed artifacts
      fprintf('Loading resumable artifacts file: %s...\n',resumeManArtFT_file);
      fprintf('\nIMPORTANT: You must repair the same channels as last time or artifacts may be different!\n');
      load(resumeManArtFT_file,'cfg_manArt');
      
      % make sure number of trials in data.cfg and cfg_manArt match
      if any(size(cfg_manArt.trl) ~= size(data.cfg.trl)) || size(cfg_manArt.trl,1) ~= length(data.trial)
        warning('Resumable artifacts file does not match the number of trials in current data!\nStarting artifact checking over.');
        ana.artifact.resumeManArtFT = false;
      end
    else
      warning('Resumable artifacts file does not exist! %s\nStarting artifact checking over.',resumeManArtFT_file);
      ana.artifact.resumeManArtFT = false;
    end
  end
  
  if ~isfield(ana.artifact,'trlpadding')
    ana.artifact.trlpadding = 0;
  end
  if ~isfield(ana.artifact,'artpadding')
    ana.artifact.artpadding = 0.1;
  end
  if ~isfield(ana.artifact,'fltpadding')
    ana.artifact.fltpadding = 0;
  end
  
  keepRepairingChannels = true;
  while keepRepairingChannels
    if ~exist('badChan_str','var')
      badChan_str = {};
    end
    
    if rejArt_ftICA
      fprintf('\nIMPORTANT!! ICA will not run properly if you repair channels and then include those in ICA.\n');
      fprintf('\tSo, it is ok to repair bad channels, but DO NOT include them in ICA.\n');
      fprintf('\tSee more information about this here: http://sccn.ucsd.edu/pipermail/eeglablist/2003/000182.html\n');
    end
    
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\nDo you want to see whether there are any REALLY BAD channels to repair? (1 or 0, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      cfgChannelRepair.viewmode = 'butterfly';
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\nDo you want to plot all channels (1) or a particular subset of channels (0)? (1 or 0, then press ''return''):\n\n');
      end
      if repair_chanNum
        cfgChannelRepair.channel = 'all';
      else
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
        fprintf('(Manually close any empty figure windows.)\n');
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('\nUse the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('Use the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\nWere there channels to repair? (1 or 0, then press ''return''):\n\n');
      end
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        fprintf('(Manually close any empty figure windows.)\n');
        badChan_str = cat(1,badChan_str,badchannel);
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
  
  if ~ana.artifact.resumeManArtFT
    ft_autoCheckArt = true;
    ft_autoCheckArtNum = 0;
    while ft_autoCheckArt
      ft_autoCheckArtNum = ft_autoCheckArtNum + 1;
      if ft_autoCheckArtNum > 1
        ft_autoCheckArt_prompt = [];
        while isempty(ft_autoCheckArt_prompt) || (ft_autoCheckArt_prompt ~= 0 && ft_autoCheckArt_prompt ~= 1)
          ft_autoCheckArt_prompt = input(sprintf('\nDo you want to run FieldTrip artifact auto-check again? This will reset previously marked artifacts. (1 or 0, then press ''return''):\n\n'));
        end
        if ft_autoCheckArt_prompt
          fprintf('\tThis is run number %d through the FieldTrip artifact auto-check.\n',ft_autoCheckArtNum);
        else
          %ft_autoCheckArt = false;
          break
        end
      end
      
      ft_customZvals_prompt = [];
      while isempty(ft_customZvals_prompt) || (ft_customZvals_prompt ~= 0 && ft_customZvals_prompt ~= 1)
        ft_customZvals_prompt = input('\nDo you want to set your own artifact thresholds (1) or use the defaults (0)? (1 or 0, then press ''return''):\n\n');
      end
      
      if ft_customZvals_prompt
        if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
          threshmin = 1;
          while threshmin >= 0
            threshmin = input(sprintf('\nAt what minimum (negative) voltage do you want to threshold (default=%d)?\n\n',threshmin_default));
          end
          if isempty(threshmin)
            threshmin = threshmin_default;
          end
          ana.artifact.threshmin = threshmin;
          
          threshmax = -1;
          while threshmax <= 0
            threshmax = input(sprintf('\nAt what maximum (positive) voltage do you want to threshold (default=%d)?\n\n',threshmax_default));
          end
          if isempty(threshmax)
            threshmax = threshmax_default;
          end
          ana.artifact.threshmax = threshmax;
        end
        
        basic_art_z = -1;
        while basic_art_z <= 0
          basic_art_z = input(sprintf('\nAt what z-value threshold do you want to check BASIC artifacts (default=%d)?\n\n',basic_art_z_default));
        end
        if isempty(basic_art_z)
          basic_art_z = basic_art_z_default;
        end
        ana.artifact.basic_art_z = basic_art_z;
        
        muscle_art_z = -1;
        while muscle_art_z <= 0
          muscle_art_z = input(sprintf('\nAt what z-value threshold do you want to check MUSCLE artifacts (default=%d)?\n\n',muscle_art_z_default));
        end
        if isempty(muscle_art_z)
          muscle_art_z = muscle_art_z_default;
        end
        ana.artifact.muscle_art_z = muscle_art_z;
        
        jump_art_z = -1;
        while jump_art_z <= 0
          jump_art_z = input(sprintf('\nAt what z-value threshold do you want to check JUMP artifacts (default=%d)?\n\n',jump_art_z_default));
        end
        if isempty(jump_art_z)
          jump_art_z = jump_art_z_default;
        end
        ana.artifact.jump_art_z = jump_art_z;
      end
      
      ft_autoCheckArt_interactive_default = 'no';
      ft_autoCheckArt_interactive = -1;
      while ft_autoCheckArt_interactive < 0 || (ft_autoCheckArt_interactive ~= 0 && ft_autoCheckArt_interactive ~= 1)
        ft_autoCheckArt_interactive = input('\nDo you want to run FieldTrip artifact auto-check in interactive mode (default=0)? (1 or 0, then press ''return''):\n\n');
        if isempty(ft_autoCheckArt_interactive)
          break
        end
      end
      if isempty(ft_autoCheckArt_interactive) || ft_autoCheckArt_interactive == 0
        ft_autoCheckArt_interactive = ft_autoCheckArt_interactive_default;
      elseif ft_autoCheckArt_interactive == 1
        ft_autoCheckArt_interactive = 'yes';
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for threshold artifacts - only if using EGI's HydroCel GSN
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
        cfg = [];
        cfg.continuous = 'no';
        % get the trial definition for automated FT artifact rejection
        cfg.trl = ft_findcfg(data.cfg,'trl');
        
        % exclude eye channels - assumes we're using EGI's HCGSN
        cfg.artfctdef.threshold.channel = {'all', ...
          '-E48', '-E128', '-E127', '-E126', '-E125', '-E119', ...
          '-E43', '-E32', '-E25', '-E21', '-E17', '-E14', '-E8', '-E1', '-E120', ...
          '-E26', '-E22', '-E15', '-E9', '-E2', ...
          '-E23', '-E18', '-E16', '-E10', '-E3', ...
          '-E19', '-E11', '-E4'};
        cfg.artfctdef.threshold.bpfilter = 'yes';
        cfg.artfctdef.threshold.bpfreq = [0.3 30];
        cfg.artfctdef.threshold.bpfiltord = 4;
        
        cfg.artfctdef.threshold.min = ana.artifact.threshmin;
        cfg.artfctdef.threshold.max = ana.artifact.threshmax;
        
        fprintf('\nUsing EGI HydroCel GSN...\nChecking for voltages above %d uV and below %d uV (excludes eye channels and neighbors)...\n',cfg.artfctdef.threshold.max,cfg.artfctdef.threshold.min);
        
        % auto mark zvalue artifacts
        [cfg, artifact_thresh] = ft_artifact_threshold(cfg, data);
      else
        warning('Not using EGI HydroCel GSN 128/129 electrode file! Threshold artifacts are not being assessed!!');
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for zvalue artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg = [];
      cfg.continuous = 'no';
      % get the trial definition for automated FT artifact rejection
      cfg.trl = ft_findcfg(data.cfg,'trl');
      
      if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
        % cfg.artfctdef.zvalue.channel = {'all', '-E25', '-E8', '-E127', '-E126', '-E128', '-E125'};
        cfg.artfctdef.zvalue.channel = {'all', ...
          '-E48', '-E128', '-E127', '-E126', '-E125', '-E119', ...
          '-E43', '-E32', '-E25', '-E21', '-E17', '-E14', '-E8', '-E1', '-E120', ...
          '-E26', '-E22', '-E15', '-E9', '-E2', ...
          '-E23', '-E18', '-E16', '-E10', '-E3', ...
          '-E19', '-E11', '-E4'};
        exclStr = ' (excludes eye channels and neighbors)';
      else
        cfg.artfctdef.zvalue.channel = 'all';
        exclStr = '';
      end
      cfg.artfctdef.zvalue.cutoff = ana.artifact.basic_art_z;
      cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
      cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
      %cfg.artfctdef.zvalue.fltpadding = 0;
      
      % algorithmic parameters
      cfg.artfctdef.zvalue.absdiff = 'yes';
      
      % interactive artifact viewer
      cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
      
      fprintf('\nChecking for (basic) zvalue artifacts at z=%d%s...\n',cfg.artfctdef.zvalue.cutoff,exclStr);
      
      % auto mark zvalue artifacts
      [cfg, artifact_zvalue] = ft_artifact_zvalue(cfg, data);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for muscle artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg = [];
      cfg.continuous = 'no';
      % get the trial definition for automated FT artifact rejection
      cfg.trl = ft_findcfg(data.cfg,'trl');
      
      % cutoff and padding
      % select a set of channels on which to run the artifact detection
      if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
        cfg.artfctdef.zvalue.channel = {'all', '-E25', '-E8', '-E127', '-E126', '-E128', '-E125'};
        exclStr = ' (excludes eye channels)';
      else
        cfg.artfctdef.zvalue.channel = 'all';
        exclStr = '';
      end
      cfg.artfctdef.zvalue.cutoff      = ana.artifact.muscle_art_z;
      cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
      if strcmp(cfg.continuous,'yes')
        cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      elseif strcmp(cfg.continuous,'no')
        cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      end
      cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
      
      % algorithmic parameters
      cfg.artfctdef.zvalue.bpfilter    = 'yes';
      if data.fsample/2 < 140
        cfg.artfctdef.zvalue.bpfreq      = [110 (data.fsample/2 - 1)];
      else
        cfg.artfctdef.zvalue.bpfreq      = [110 140];
      end
      cfg.artfctdef.zvalue.bpfiltord   = 6;
      cfg.artfctdef.zvalue.bpfilttype  = 'but';
      cfg.artfctdef.zvalue.hilbert     = 'yes';
      cfg.artfctdef.zvalue.boxcar      = 0.2;
      
      % interactive artifact viewer
      cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
      
      fprintf('\nChecking for muscle artifacts at z=%d%s...\n',cfg.artfctdef.zvalue.cutoff,exclStr);
      
      % auto mark muscle artifacts
      [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,data);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for jump artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg = [];
      cfg.continuous = 'no';
      % get the trial definition for automated FT artifact rejection
      cfg.trl = ft_findcfg(data.cfg,'trl');
      
      % cutoff and padding
      % select a set of channels on which to run the artifact detection
      cfg.artfctdef.zvalue.channel = 'all';
      cfg.artfctdef.zvalue.cutoff = ana.artifact.jump_art_z;
      cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
      cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
      %cfg.artfctdef.zvalue.fltpadding = 0;
      
      % algorithmic parameters
      cfg.artfctdef.zvalue.cumulative = 'yes';
      cfg.artfctdef.zvalue.medianfilter = 'yes';
      cfg.artfctdef.zvalue.medianfiltord = 9;
      cfg.artfctdef.zvalue.absdiff = 'yes';
      
      % interactive artifact viewer
      cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
      
      fprintf('\nChecking for jump artifacts at z=%d...\n',cfg.artfctdef.zvalue.cutoff);
      
      % auto mark jump artifacts
      [cfg, artifact_jump] = ft_artifact_zvalue(cfg,data);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % manual inspection of artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg_manArt = [];
      cfg_manArt.artfctdef.zvalue.artifact = artifact_zvalue;
      cfg_manArt.artfctdef.muscle.artifact = artifact_muscle;
      cfg_manArt.artfctdef.jump.artifact = artifact_jump;
      if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
        cfg_manArt.artfctdef.threshold.artifact = artifact_thresh;
      end
      
      % get the trial definition for automated FT artifact rejection
      cfg_manArt.trl = ft_findcfg(data.cfg,'trl');
      
      cfg_manArt.continuous = 'no';
      %cfg_manArt.viewmode = 'butterfly';
      cfg_manArt.viewmode = 'vertical';
      cfg_manArt.elecfile = elecfile;
      cfg_manArt.plotlabels = 'some';
      cfg_manArt.ylim = vert_ylim;
      
      % use cursor drag and click to mark artifacts;
      % use arrows to advance to next trial;
      % use the q key to quit the data browser
      
      % manual rejection
      fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
      fprintf('\n\nManual artifact rejection:\n');
      fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
      fprintf('\tUse arrows to move to next trial.\n');
      if strcmp(cfg_manArt.viewmode,'butterfly')
        fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
      end
      fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
      fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n');
      
      if rejArt_ftICA
        fprintf('\nNB: Before running ICA, you must manually reject artifacts that are not consistent across trials.\n');
        fprintf('\tDO NOT reject blinks if you want to remove them with ICA!\n\tPlease reject inconsistent artifacts now.\n\n');
      end
      
      cfg_manArt = ft_databrowser(cfg_manArt,data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      fprintf('\nBacking up artifact data to %s.\n',resumeManArtFT_file);
      fprintf('\tIf MATLAB crashes before finishing, you can resume without re-checking artifacts by setting ana.artifact.resumeManArtFT=true in your main file.\n');
      save(resumeManArtFT_file,'cfg_manArt');
    end
  else
    % use cursor drag and click to mark artifacts;
    % use arrows to advance to next trial;
    % use the q key to quit the data browser
    
    % manual rejection
    fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
    fprintf('\n\nManual artifact rejection:\n');
    fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
    fprintf('\tUse arrows to move to next trial.\n');
    if strcmp(cfg_manArt.viewmode,'butterfly')
      fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
    end
    fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
    fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n');
    
    if rejArt_ftICA
      fprintf('\nNB: Before running ICA, you must manually reject artifacts that are not consistent across trials.\n');
      fprintf('\tDO NOT reject blinks if you want to remove them with ICA!\n\tPlease reject inconsistent artifacts now.\n\n');
    end
    
    fprintf('\nIMPORTANT:\n\tYou loaded saved artifact data. You can exit the artifact browser and the artifacts you marked last time will be rejected.\n');
    cfg_manArt = ft_databrowser(cfg_manArt,data);
    % bug when calling rejectartifact right after databrowser, pause first
    pause(1);
    fprintf('\nBacking up artifact data to %s.\n',resumeManArtFT_file);
    fprintf('\tIf MATLAB crashes before finishing, you can resume without re-checking artifacts by setting ana.artifact.resumeManArtFT=true in your main file.\n');
    save(resumeManArtFT_file,'cfg_manArt');
  end
  
  % initialize to store whether there was an artifact for each trial
  badEv = zeros(size(data.sampleinfo,1), 1);
  % find out what kind of artifacts we're dealing with
  fn = fieldnames(cfg_manArt.artfctdef);
  theseArt = {};
  for i = 1:length(fn)
    if isstruct(cfg_manArt.artfctdef.(fn{i})) && isfield(cfg_manArt.artfctdef.(fn{i}),'artifact') && ~isempty(cfg_manArt.artfctdef.(fn{i}).artifact)
      theseArt = cat(2,theseArt,fn{i});
    end
  end
  % find out which samples were marked as artifacts
  if ~isempty(theseArt)
    artSamp = zeros(data.sampleinfo(end,2),1);
    for i = 1:length(theseArt)
      for j = 1:size(cfg_manArt.artfctdef.(theseArt{i}).artifact,1)
        % mark that it was a particular type of artifact
        artSamp(cfg_manArt.artfctdef.(theseArt{i}).artifact(j,1):cfg_manArt.artfctdef.(theseArt{i}).artifact(j,2)) = find(ismember(theseArt,theseArt{i}));
      end
    end
    % save a list of trials with artifact status
    for k = 1:size(data.sampleinfo,1)
      if any(artSamp(data.sampleinfo(k,1):data.sampleinfo(k,2)) > 0)
        badEv(k,1) = 1;
      end
    end
  end
  
  % reject the artifacts (complete or parial rejection)
  cfg_manArt.artfctdef.reject = 'complete';
  data = ft_rejectartifact(cfg_manArt, data);
  
  keepRepairingChannels = true;
  while keepRepairingChannels
    if ~exist('badChan_str','var')
      badChan_str = {};
    end
    
    if rejArt_ftICA
      fprintf('\nIMPORTANT!! ICA will not run properly if you repair channels and then include those in ICA.\n');
      fprintf('\tSo, it is ok to repair bad channels, but DO NOT include them in ICA.\n');
      fprintf('\tSee more information about this here: http://sccn.ucsd.edu/pipermail/eeglablist/2003/000182.html\n');
    end
    
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\nDo you want to see whether there are channels to repair? (1 or 0, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      %cfgChannelRepair.viewmode = 'butterfly';
      
      % viewmode?
      repair_viewmode = [];
      while isempty(repair_viewmode) || (repair_viewmode ~= 0 && repair_viewmode ~= 1)
        repair_viewmode = input('\nDo you want to plot in butterfly (1) or vertical (0) mode? (1 or 0, then press ''return''):\n\n');
      end
      if repair_viewmode
        cfgChannelRepair.viewmode = 'butterfly';
      else
        cfgChannelRepair.viewmode = 'vertical';
        cfgChannelRepair.ylim = vert_ylim;
      end
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\nDo you want to plot all channels (1) or a particular subset of channels (0)? (1 or 0, then press ''return''):\n\n');
      end
      if repair_chanNum
        cfgChannelRepair.channel = 'all';
      else
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
        fprintf('(Manually close any empty figure windows.)\n');
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('\nUse the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('\nUse the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\nWere there channels to repair? (1 or 0, then press ''return''):\n\n');
      end
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        fprintf('(Manually close any empty figure windows.)\n');
        badChan_str = cat(1,badChan_str,badchannel);
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
  
  if ~exist('badChan_str','var')
    badChan_str = {};
  end
  
  if ~isfield(ana.artifact,'trlpadding')
    ana.artifact.trlpadding = 0;
  end
  if ~isfield(ana.artifact,'artpadding')
    ana.artifact.artpadding = 0.1;
  end
  if ~isfield(ana.artifact,'fltpadding')
    ana.artifact.fltpadding = 0;
  end
  
  if isfield(ana.artifact,'basic_art_z_postICA')
    basic_art_z_postICA_default = ana.artifact.basic_art_z_postICA;
  elseif ~isfield(ana.artifact,'basic_art_z_postICA')
    basic_art_z_postICA_default = 30;
    ana.artifact.basic_art_z_postICA = basic_art_z_postICA_default;
  end
  if isfield(ana.artifact,'muscle_art_z_postICA')
    muscle_art_z_postICA_default = ana.artifact.muscle_art_z_postICA;
  elseif ~isfield(ana.artifact,'muscle_art_z_postICA')
    muscle_art_z_postICA_default = 40;
    ana.artifact.muscle_art_z_postICA = muscle_art_z_postICA_default;
  end
  if isfield(ana.artifact,'jump_art_z_postICA')
    jump_art_z_postICA_default = ana.artifact.jump_art_z_postICA;
  elseif ~isfield(ana.artifact,'jump_art_z_postICA')
    jump_art_z_postICA_default = 50;
    ana.artifact.jump_art_z_postICA = jump_art_z_postICA_default;
  end
  
  if isfield(ana.artifact,'threshmin_postICA')
    threshmin_postICA_default = ana.artifact.threshmin_postICA;
  elseif ~isfield(ana.artifact,'threshmin_postICA')
    threshmin_postICA_default = -150;
    ana.artifact.threshmin_postICA = threshmin_postICA_default;
  end
  if isfield(ana.artifact,'threshmax_postICA')
    threshmax_postICA_default = ana.artifact.threshmax_postICA;
  elseif ~isfield(ana.artifact,'threshmax_postICA')
    threshmax_postICA_default = 150;
    ana.artifact.threshmax_postICA = threshmax_postICA_default;
  end
  
  % set the file to save after running ICA, in case MATLAB crashes
  resumeICACompFT_file = fullfile(dirs.saveDirRaw,subject,sesName,sprintf('%s_%s_ICACompFT.mat',subject,sesName));
  if ~isfield(ana.artifact,'resumeICACompFT')
    ana.artifact.resumeICACompFT = false;
  end
  
  if ana.artifact.resumeICACompFT
    if exist(resumeICACompFT_file,'file')
      % load the manually processed artifacts
      fprintf('Loading resumable ICA components file: %s...\n',resumeICACompFT_file);
      fprintf('\nIMPORTANT: You must have rejected the same channels and artifacts as last time or components may be different!\n');
      load(resumeICACompFT_file,'comp');
    else
      warning('Resumable ICA components file does not exist! %s\nStarting ICA over.',resumeICACompFT_file);
      ana.artifact.resumeICACompFT = false;
    end
  end

  if ~ana.artifact.resumeICACompFT
    cfg_ica = [];
    cfg_ica.method = 'runica';
    %cfg_ica.demean = 'no';
    
    fprintf('\nIf you still have really bad channels, or have repaired channels, you must exclude them from ICA.\n');
    if ~isempty(badChan_str)
      fprintf('\n\tIMPORTANT! You have repaired channels:%s\n\n',sprintf(repmat(' %s',1,length(badChan_str)),badChan_str{:}));
      ica_chanNum = 0;
      fprintf('\tTherefore, you must run ICA on a subset of channels.\n');
    else
      ica_chanNum = [];
      fprintf('\nWe believe that you have NOT repaired any channels. Thus, you can run ICA on all channels (option ''1'').\n');
      fprintf('\tBut if that somehow is not the case, you must run ICA on a subset of channels (option ''0'').\n');
    end
    
    %ica_chanNum = [];
    while isempty(ica_chanNum) || (ica_chanNum ~= 0 && ica_chanNum ~= 1)
      ica_chanNum = input('\nDo you want to run ICA on all channels (1) or only a subset of channels (0)? (1 or 0, then press ''return''):\n\n');
    end
    if ica_chanNum
      cfg_ica.channel = 'all';
    else
      fprintf('\tOnce you see the channel selector:\n');
      fprintf('\t\t1. Add all channels to the right-side list.\n');
      fprintf('\t\t2. Remove the bad channel from the right-side list, into the left-side list.\n');
      fprintf('\t\t3. Manually close any empty figure windows.\n');
      channel = ft_channelselection('gui', data.label);
      cfg_ica.channel = channel;
    end
    
    comp = ft_componentanalysis(cfg_ica,data);
    
    fprintf('\nBacking up ICA component data to %s.\n',resumeICACompFT_file);
    fprintf('\tIf MATLAB crashes before finishing, you can resume without re-running ICA by setting ana.artifact.resumeICACompFT=true in your main file.\n');
    save(resumeICACompFT_file,'comp');
  end
  
  keepChoosingICAcomps = true;
  while keepChoosingICAcomps
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
    cfg.layout = elecfile;
    cfg.elecfile = elecfile;
    cfg.ylim = vert_ylim;
    ft_databrowser(cfg,comp);
    % bug when calling rejectartifact right after databrowser, pause first
    pause(1);
    
    fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
    %fprintf('\n\nViewing the first %d components.\n',nComponents);
    fprintf('ICA component browsing:\n');
    fprintf('\t1. Look for patterns that are indicative of artifacts.\n');
    fprintf('\t\tPress the ''channel >'' button to see the next set of components.\n');
    fprintf('\t\tComponents may not be numbered, so KEEP TRACK of where you are (top component has the lowest number). Write down component numbers for rejection.\n');
    fprintf('\t2. Manually close the components window when finished browsing.\n');
    
    rej_comp = [];
    while isempty(rej_comp) || (rej_comp ~= 0 && rej_comp ~= 1)
      rej_comp = input('Were there components to reject? (1 or 0, then press ''return''):\n\n');
    end
    if rej_comp
      % prompt the user for the component numbers to reject
      componentsToReject = input(sprintf('\t3. Type component numbers to reject (on a single line) and press ''return'',\n\teven if these instructions move up due to output while browsing components (e.g., ''1, 4, 11'' without quotes):\n\n'),'s');
      
      % reject the bad components
      if ~isempty(componentsToReject)
        cfg = [];
        cfg.component = str2double(regexp(componentsToReject,'\d*','match')');
        data_ica_rej = ft_rejectcomponent(cfg, comp, data);
      end
    else
      data_ica_rej = data;
    end
    
    % another auto search for artifacts
    
    % use cursor drag and click to mark artifacts;
    % use arrows to advance to next trial;
    % use the q key to quit the data browser
    
    ft_autoCheckArt = true;
    ft_autoCheckArtNum = 0;
    while ft_autoCheckArt
      ft_autoCheckArtNum = ft_autoCheckArtNum + 1;
      if ft_autoCheckArtNum > 1
        ft_autoCheckArt_prompt = [];
        while isempty(ft_autoCheckArt_prompt) || (ft_autoCheckArt_prompt ~= 0 && ft_autoCheckArt_prompt ~= 1)
          ft_autoCheckArt_prompt = input(sprintf('\nDo you want to run FieldTrip artifact auto-check again? This will reset previously marked artifacts. (1 or 0, then press ''return''):\n\n'));
        end
        if ft_autoCheckArt_prompt
          fprintf('\tThis is run number %d through the FieldTrip artifact auto-check.\n',ft_autoCheckArtNum);
        else
          %ft_autoCheckArt = false;
          break
        end
      end
      
      fprintf('\n\nFinal round of manual artifact rejection:\n');
      
      ft_customZvals_prompt = [];
      while isempty(ft_customZvals_prompt) || (ft_customZvals_prompt ~= 0 && ft_customZvals_prompt ~= 1)
        ft_customZvals_prompt = input('\nDo you want to set your own artifact thresholds (1) or use the defaults (0)? (1 or 0, then press ''return''):\n\n');
      end
      
      if ft_customZvals_prompt
        if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
          threshmin_postICA = 1;
          while threshmin_postICA >= 0
            threshmin_postICA = input(sprintf('\nAt what minimum (negative) voltage do you want to threshold (default=%d)?\n\n',threshmin_postICA_default));
          end
          if isempty(threshmin_postICA)
            threshmin_postICA = threshmin_postICA_default;
          end
          ana.artifact.threshmin_postICA = threshmin_postICA;
          
          threshmax_postICA = -1;
          while threshmax_postICA <= 0
            threshmax_postICA = input(sprintf('\nAt what maximum (positive) voltage do you want to threshold (default=%d)?\n\n',threshmax_postICA_default));
          end
          if isempty(threshmax_postICA)
            threshmax_postICA = threshmax_postICA_default;
          end
          ana.artifact.threshmax_postICA = threshmax_postICA;
        end
        
        basic_art_z = -1;
        while basic_art_z <= 0
          basic_art_z = input(sprintf('\nAt what z-value threshold do you want to check BASIC artifacts (default=%d)?\n\n',basic_art_z_postICA_default));
        end
        if isempty(basic_art_z)
          basic_art_z = basic_art_z_postICA_default;
        end
        ana.artifact.basic_art_z_postICA = basic_art_z;
        
        muscle_art_z = -1;
        while muscle_art_z <= 0
          muscle_art_z = input(sprintf('\nAt what z-value threshold do you want to check MUSCLE artifacts (default=%d)?\n\n',muscle_art_z_postICA_default));
        end
        if isempty(muscle_art_z)
          muscle_art_z = muscle_art_z_postICA_default;
        end
        ana.artifact.muscle_art_z_postICA = muscle_art_z;
        
        jump_art_z = -1;
        while jump_art_z <= 0
          jump_art_z = input(sprintf('\nAt what z-value threshold do you want to check JUMP artifacts (default=%d)?\n\n',jump_art_z_postICA_default));
        end
        if isempty(jump_art_z)
          jump_art_z = jump_art_z_postICA_default;
        end
        ana.artifact.jump_art_z_postICA = jump_art_z;
      end
      
      ft_autoCheckArt_interactive_default = 'no';
      ft_autoCheckArt_interactive = -1;
      while ft_autoCheckArt_interactive < 0 || (ft_autoCheckArt_interactive ~= 0 && ft_autoCheckArt_interactive ~= 1)
        ft_autoCheckArt_interactive = input('\nDo you want to run FieldTrip artifact auto-check in interactive mode (default=0)? (1 or 0, then press ''return''):\n\n');
        if isempty(ft_autoCheckArt_interactive)
          break
        end
      end
      if isempty(ft_autoCheckArt_interactive) || ft_autoCheckArt_interactive == 0
        ft_autoCheckArt_interactive = ft_autoCheckArt_interactive_default;
      elseif ft_autoCheckArt_interactive == 1
        ft_autoCheckArt_interactive = 'yes';
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for threshold artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
        cfg = [];
        cfg.continuous = 'no';
        % get the trial definition for automated FT artifact rejection
        cfg.trl = ft_findcfg(data.cfg,'trl');
        
        % don't exclude eye channels because we want to reject any blinks
        % that ICA didn't catch
        cfg.artfctdef.threshold.channel = {'all'};
        exclStr = '';
        
        % % exclude eye channels - assumes we're using EGI's HCGSN
        % cfg.artfctdef.zvalue.channel = {'all', '-E25', '-E8', '-E127', '-E126', '-E128', '-E125'};
        % exclStr = ' (excludes eye channels)';
        
        % % exclude eye channels - assumes we're using EGI's HCGSN
        % cfg.artfctdef.threshold.channel = {'all', ...
        %   '-E48', '-E128', '-E127', '-E126', '-E125', '-E119', ...
        %   '-E43', '-E32', '-E25', '-E21', '-E17', '-E14', '-E8', '-E1', '-E120', ...
        %   '-E26', '-E22', '-E15', '-E9', '-E2', ...
        %   '-E23', '-E18', '-E16', '-E10', '-E3', ...
        %   '-E19', '-E11', '-E4'};
        % exclStr = ' (excludes eye channels and neighbors)';
        
        cfg.artfctdef.threshold.bpfilter = 'yes';
        cfg.artfctdef.threshold.bpfreq = [0.3 30];
        cfg.artfctdef.threshold.bpfiltord = 4;
        
        cfg.artfctdef.threshold.min = ana.artifact.threshmin_postICA;
        cfg.artfctdef.threshold.max = ana.artifact.threshmax_postICA;
        
        fprintf('\nUsing EGI HydroCel GSN...\nChecking for voltages above %d uV and below %d uV%s...\n',cfg.artfctdef.threshold.max,cfg.artfctdef.threshold.min,exclStr);
        
        % auto mark zvalue artifacts
        [cfg, artifact_thresh] = ft_artifact_threshold(cfg, data_ica_rej);
      else
        warning('Not using EGI HydroCel GSN 128/129 electrode file! Threshold artifacts are not being assessed!!');
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for zvalue artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg = [];
      cfg.continuous = 'no';
      % get the trial definition for automated FT artifact rejection
      cfg.trl = ft_findcfg(data_ica_rej.cfg,'trl');
      
      cfg.artfctdef.zvalue.channel = 'all';
      cfg.artfctdef.zvalue.cutoff = ana.artifact.basic_art_z_postICA;
      cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
      cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
      %cfg.artfctdef.zvalue.fltpadding = 0;
      
      % interactive artifact viewer
      cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
      
      fprintf('Checking for (basic) zvalue artifacts at z=%d...\n',cfg.artfctdef.zvalue.cutoff);
      
      % auto mark some artifacts
      [cfg, artifact_zvalue] = ft_artifact_zvalue(cfg, data_ica_rej);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for muscle artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg = [];
      cfg.continuous = 'no';
      % get the trial definition for automated FT artifact rejection
      cfg.trl = ft_findcfg(data_ica_rej.cfg,'trl');
      
      % cutoff and padding
      % select a set of channels on which to run the artifact detection
      cfg.artfctdef.zvalue.channel = 'all';
      cfg.artfctdef.zvalue.cutoff      = ana.artifact.muscle_art_z_postICA;
      cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
      if strcmp(cfg.continuous,'yes')
        cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      elseif strcmp(cfg.continuous,'no')
        cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      end
      cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
      
      % algorithmic parameters
      cfg.artfctdef.zvalue.bpfilter    = 'yes';
      if data_ica_rej.fsample/2 < 140
        cfg.artfctdef.zvalue.bpfreq      = [110 (data_ica_rej.fsample/2 - 1)];
      else
        cfg.artfctdef.zvalue.bpfreq      = [110 140];
      end
      cfg.artfctdef.zvalue.bpfiltord   = 6;
      cfg.artfctdef.zvalue.bpfilttype  = 'but';
      cfg.artfctdef.zvalue.hilbert     = 'yes';
      cfg.artfctdef.zvalue.boxcar      = 0.2;
      
      % interactive artifact viewer
      cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
      
      fprintf('\nChecking for muscle artifacts at z=%d...\n',cfg.artfctdef.zvalue.cutoff);
      
      % auto mark muscle artifacts
      [cfg, artifact_muscle] = ft_artifact_zvalue(cfg, data_ica_rej);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % look for jump artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg = [];
      cfg.continuous = 'no';
      % get the trial definition for automated FT artifact rejection
      cfg.trl = ft_findcfg(data_ica_rej.cfg,'trl');
      
      % cutoff and padding
      % select a set of channels on which to run the artifact detection
      cfg.artfctdef.zvalue.channel = 'all';
      cfg.artfctdef.zvalue.cutoff = ana.artifact.jump_art_z_postICA;
      cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
      cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
      %cfg.artfctdef.zvalue.fltpadding = 0;
      
      % algorithmic parameters
      cfg.artfctdef.zvalue.cumulative = 'yes';
      cfg.artfctdef.zvalue.medianfilter = 'yes';
      cfg.artfctdef.zvalue.medianfiltord = 9;
      cfg.artfctdef.zvalue.absdiff = 'yes';
      
      % interactive artifact viewer
      cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
      
      fprintf('\nChecking for jump artifacts at z=%d...\n',cfg.artfctdef.zvalue.cutoff);
      
      % auto mark jump artifacts
      [cfg, artifact_jump] = ft_artifact_zvalue(cfg, data_ica_rej);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % manual inspection of artifacts
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg_manArt = [];
      cfg_manArt.artfctdef.zvalue.artifact = artifact_zvalue;
      cfg_manArt.artfctdef.muscle.artifact = artifact_muscle;
      cfg_manArt.artfctdef.jump.artifact = artifact_jump;
      if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
        cfg_manArt.artfctdef.threshold.artifact = artifact_thresh;
      end
      
      % another manual search of the data for artifacts
      
      %cfg_manArt.viewmode = 'butterfly';
      cfg_manArt.viewmode = 'vertical';
      cfg_manArt.continuous = 'no';
      cfg_manArt.elecfile = elecfile;
      cfg_manArt.plotlabels = 'some';
      cfg_manArt.ylim = vert_ylim;
      
      fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue)),eventValue{:}));
      fprintf('\n\nFinal round of manual artifact rejection:\n');
      fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
      fprintf('\tUse arrows to move to next trial.\n');
      if strcmp(cfg_manArt.viewmode,'butterfly')
        fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
      end
      fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
      fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n');
      
      cfg_manArt = ft_databrowser(cfg_manArt, data_ica_rej);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
    end
    done_with_ica = [];
    while isempty(done_with_ica) || (done_with_ica ~= 0 && done_with_ica ~= 1)
      done_with_ica = input('\nAre you happy with your post-ICA results? 1 to move on to next step, 0 to redo ICA component rejection and rerun artifact detection. (1 or 0, then press ''return''):\n\n');
    end
    
    if done_with_ica
      data = data_ica_rej;
      keepChoosingICAcomps = false;
    end
  end
  
  % initialize to store whether there was an artifact for each trial
  if ~exist('badEv','var') || isempty(badEv)
    combineArtLists = false;
    %badEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
    badEv = zeros(size(data.sampleinfo,1), 1);
  else
    combineArtLists = true;
  end
  postIcaEv = zeros(size(data.sampleinfo,1), 1);
  
  %   % debug
  %   % is size(data.sampleinfo,1) == size(postIcaEv,1)
  %   % is size(badEv(badEv == 0),1) == size(postIcaEv,1)
  %   % is size(badEv(badEv == 0),1) == size(data.sampleinfo,1)
  %   %
  %   % these should be the same size
  %   if size(badEv(badEv == 0),1) ~= size(data.sampleinfo,1)
  %     keyboard
  %   end
  
  % find out what kind of artifacts we're dealing with
  fn = fieldnames(cfg_manArt.artfctdef);
  theseArt = {};
  for i = 1:length(fn)
    if isstruct(cfg_manArt.artfctdef.(fn{i})) && isfield(cfg_manArt.artfctdef.(fn{i}),'artifact') && ~isempty(cfg_manArt.artfctdef.(fn{i}).artifact)
      theseArt = cat(2,theseArt,fn{i});
    end
  end
  % find out which samples were marked as artifacts
  if ~isempty(theseArt)
    artSamp = zeros(data.sampleinfo(end,2),1);
    for i = 1:length(theseArt)
      for j = 1:size(cfg_manArt.artfctdef.(theseArt{i}).artifact,1)
        % mark that it was a particular type of artifact
        artSamp(cfg_manArt.artfctdef.(theseArt{i}).artifact(j,1):cfg_manArt.artfctdef.(theseArt{i}).artifact(j,2)) = find(ismember(theseArt,theseArt{i}));
      end
    end
    % save a list of trials with artifact status
    for k = 1:size(data.sampleinfo,1)
      if any(artSamp(data.sampleinfo(k,1):data.sampleinfo(k,2)) > 0)
        postIcaEv(k,1) = 1;
      end
    end
  end
  
  if combineArtLists
    % put the new artifacts into the old list
    rCount = 0;
    for i = 1:size(badEv,1)
      %if badEv(i,2) == 0
      if badEv(i,1) == 0
        rCount = rCount + 1;
        %if postIcaEv(rCount,2) == 1
        if postIcaEv(rCount) == 1
          %badEv(i,2) = 1;
          badEv(i,1) = 1;
        end
      end
    end
  else
    badEv = postIcaEv;
  end
  
  % reject the artifacts (complete or parial rejection)
  cfg_manArt.artfctdef.remove = 'complete';
  % and reject
  data = ft_rejectartifact(cfg_manArt,data);
  
  keepRepairingChannels = true;
  while keepRepairingChannels
    if ~exist('badChan_str','var')
      badChan_str = {};
    end
    
    % see if there were any channels to repair first
    rejArt_repair = [];
    while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
      rejArt_repair = input('\nDo you want to see whether there are channels to repair? (1 or 0, then press ''return''):\n\n');
    end
    
    if rejArt_repair
      cfgChannelRepair = [];
      cfgChannelRepair.continuous = 'no';
      cfgChannelRepair.elecfile = elecfile;
      cfgChannelRepair.viewmode = 'butterfly';
      
      % subset?
      repair_chanNum = [];
      while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
        repair_chanNum = input('\nDo you want to plot all channels (1) or a particular subset of channels (0)? (1 or 0, then press ''return''):\n\n');
      end
      if repair_chanNum
        cfgChannelRepair.channel = 'all';
      else
        channel = ft_channelselection('gui', data.label);
        cfgChannelRepair.channel = channel;
        fprintf('(Manually close any empty figure windows.)\n');
      end
      
      if strcmp(cfgChannelRepair.viewmode,'butterfly')
        fprintf('\nUse the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
      end
      fprintf('\nUse the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
      cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      
      rejArt_repair_really = [];
      while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
        rejArt_repair_really = input('\nWere there channels to repair? (1 or 0, then press ''return''):\n\n');
      end
      if rejArt_repair_really
        badchannel = ft_channelselection('gui', data.label);
        fprintf('(Manually close any empty figure windows.)\n');
        badChan_str = cat(1,badChan_str,badchannel);
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
  if ~isfield(ana.artifact,'trlpadding')
    ana.artifact.trlpadding = 0;
  end
  if ~isfield(ana.artifact,'artpadding')
    ana.artifact.artpadding = 0.1;
  end
  if ~isfield(ana.artifact,'fltpadding')
    ana.artifact.fltpadding = 0;
  end
  
  % get the trial definition for automated FT artifact rejection
  trl = ft_findcfg(data.cfg,'trl');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for jump artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection
  cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.cutoff = 20;
  cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
  cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
  cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
  %cfg.artfctdef.zvalue.fltpadding = 0;
  
  % algorithmic parameters
  cfg.artfctdef.zvalue.cumulative = 'yes';
  cfg.artfctdef.zvalue.medianfilter = 'yes';
  cfg.artfctdef.zvalue.medianfiltord = 9;
  cfg.artfctdef.zvalue.absdiff = 'yes';
  
  % interactive artifact viewer
  cfg.artfctdef.zvalue.interactive = 'yes';
  
  [cfg, artifact_jump] = ft_artifact_zvalue(cfg,data);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for muscle artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection
  cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.cutoff      = 40;
  cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
  if strcmp(cfg.continuous,'yes')
    cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
  elseif strcmp(cfg.continuous,'no')
    cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
  end
  cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
  
  % algorithmic parameters
  cfg.artfctdef.zvalue.bpfilter    = 'yes';
  if data.fsample/2 < 140
    cfg.artfctdef.zvalue.bpfreq      = [110 (data.fsample/2 - 1)];
  else
    cfg.artfctdef.zvalue.bpfreq      = [110 140];
  end
  cfg.artfctdef.zvalue.bpfiltord   = 6;
  cfg.artfctdef.zvalue.bpfilttype  = 'but';
  cfg.artfctdef.zvalue.hilbert     = 'yes';
  cfg.artfctdef.zvalue.boxcar      = 0.2;
  
  % interactive artifact viewer
  %cfg.artfctdef.zvalue.interactive = 'yes';
  
  [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,data);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for EOG artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
  %cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.channel = {'E127','E126','E128','E125'};
  cfg.artfctdef.zvalue.cutoff      = 4;
  cfg.artfctdef.zvalue.trlpadding = ana.artifact.artpadding;
  if strcmp(cfg.continuous,'yes')
    cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
  elseif strcmp(cfg.continuous,'no')
    cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
  end
  cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
  
  % algorithmic parameters
  cfg.artfctdef.zvalue.bpfilter   = 'yes';
  cfg.artfctdef.zvalue.bpfilttype = 'but';
  cfg.artfctdef.zvalue.bpfreq     = [1 15];
  cfg.artfctdef.zvalue.bpfiltord  = 4;
  cfg.artfctdef.zvalue.hilbert    = 'yes';
  
  % interactive artifact viewer
  %cfg.artfctdef.zvalue.interactive = 'yes';
  
  [cfg, artifact_EOG] = ft_artifact_zvalue(cfg,data);
  
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
