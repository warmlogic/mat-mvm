function [ft_raw] = seg2ft(dataroot,nsFileExt,subject,session,eventValue,prepost,elecfile,artifactType)
%SEG2FT: take segmented EEG data and put it in FieldTrip format
%
% [ft_raw] = seg2ft(dataroot,nsFileExt,subject,session,eventValue,prepost,elecfile,artifactType)
%
% Output:
%   ft_raw = struct with one field for each event value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP:
%
% Export Net Station data as either EGIS (129 channels, where the final
% channel [Cz] is the reference; extension: 'egis') or NS Simple Binary
% (extension: 'raw' or 'sbin') calibrated, including the reference channel.
% These options are in the File Export tool. create_ft_struct, the function
% that calls seg2ft, expects EGIS files to be stored in a 'ns_egis'
% directory at the level of dirs.dataDir. If using raw files, they should
% be in 'ns_raw' instead.
%
% This function can deal with event values that have zero events. It is
% probably better to check on your event count and exclude those subjects
% with zero events for any of the event values you're trying to keep before
% trying to run this script. Nonetheless, it will insert an empty event
% entry for an empty eventValue, and the subjects will be excluded when
% using mm_threshSubs.m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFACT INFORMATION:
%
% artifactType can be 'none', 'ns', 'ft_man', or 'ft_ica'; it can be one of
% those strings or a cell array of multiple strings (e.g., {'ns','ft_man'}
% to do both Net Station artifact rejection and FieldTrip manual
% rejection). 'ft_ica' also includes manual rejection after assessing
% components.
%
% 'ns' is processed first, then 'ft_man', then 'ft_ica'.  Subquent
% processing will not include earlier rejected artifacts.  Note that any
% FT artifact processing requires manual intervention, while NS artifact
% processing does not.
%
% If 'ns', this function expects to find a Net Station segment information
% file with a .bci extension; this contains artifact information. It is
% exported from Net Station using the File Export tool. To set up the tool,
% export format is metadata and check the segment information option. Run
% the tool on the file that was exported to egis/raw (i.e., the baseline
% correction file or the average rereference file). The bci files should be
% stored in a 'ns_bci' directory at the same level as 'ns_egis' or
% 'ns_raw'.
%
% If 'ft_man', a visualization of all channels for each event will appear,
% where each trial is shown one-by-one.
%
% If 'ft_ica', ICA will run on all trials across all event values.
% Individual components can be rejected after this.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB DATA
%
% If using eeglab data, no artifact detection is done and no bci file is
% expected to exist. Also, the directory structure is different and can be
% gleaned by examining the code here, but right now it is only set up to
% process Erika Nyhus's KAHN2 data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also: create_ft_struct
%

%% set the artifact processing parameters

if ischar(artifactType)
  artifactType = {artifactType};
end

artifactOpts = {'none','ns','ft_man','ft_ica'};

if any(~ismember(artifactType,artifactOpts))
  error('an artifact option was not set correctly (it was set to ''%s'')',cell2mat(artifactType(~ismember(artifactType,artifactOpts))))
end

% set artifact defaults
rejArt = 0;
rejArt_ns = 0;
rejArt_ft_man = 0;
rejArt_ft_ica = 0;

% figure out which artifact options we're using
if ~ismember('none',artifactType)
  rejArt = 1;
end
if ismember('ns',artifactType)
  rejArt_ns = 1;
end
if ismember('ft_man',artifactType)
  rejArt_ft_man = 1;
end
if ismember('ft_ica',artifactType)
  rejArt_ft_ica = 1;
end

%% set up other parameters

% make sure eventValue is set up correctly
if ~iscell(eventValue)
  eventValue = {eventValue};
end
% if length(eventValue) > 1
%   error('Expecting only one eventValue.');
% end

if ~iscell(session)
  session = {session};
end
if length(session) > 1
  append_data = struct;
end

if strcmpi(nsFileExt,'raw') || strcmpi(nsFileExt,'sbin')
  ftype = 'egi_sbin';
  nsDir = 'ns_raw';
elseif strcmpi(nsFileExt,'egis')
  ftype = 'egi_egis';
  nsDir = 'ns_egis';
elseif strcmpi(nsFileExt,'set')
  ftype = 'eeglab_set';
  nsDir = subject;
else
  error('ftype not set because extension was not properly set.');
end

% make sure the chan locs file exists
%
% if this is an EGI electrode location file included with FieldTrip, the 3
% Fid (fiduciary) points are included, meaning there are 3 non-electrodes
% in elec.label
if ~exist(elecfile,'file')
  error('Cannot find channel locs file at %s',elecfile);
else
  [cpath,cname,cext] = fileparts(elecfile);
  if strcmpi(cext,'.sfp')
    locsFormat = 'besa_sfp';
  else
    locsFormat = [];
  end
  elec = ft_read_sens(elecfile,'fileformat',locsFormat);
end
nChan_elecfile = size(elec.label,1);

%% for each session, read in the EEG file

for ses = 1:length(session)
  % set ses_str to make sure it starts with a character, not a #, etc.
  ses_str = sprintf('ses_%s',session{ses});
  
  if strcmpi(nsFileExt,'sbin') || strcmpi(nsFileExt,'raw') || strcmpi(nsFileExt,'egis')
    % make sure the EEG file exists
    nsfile = dir(fullfile(dataroot,session{ses},nsDir,[subject,'*.',nsFileExt]));
    if isempty(nsfile)
      error('Cannot find %s*.%s file in %s',subject,nsFileExt,fullfile(dataroot,session{ses},nsDir));
    elseif length(nsfile) > 1
      error('More than one %s*.%s file found in %s',subject,nsFileExt,fullfile(dataroot,session{ses},nsDir));
    elseif length(nsfile) == 1
      infile_ns = fullfile(dataroot,session{ses},nsDir,nsfile.name);
    end
    
  elseif strcmpi(nsFileExt,'set')
    % this is really just set up to analyze Erika Nyhus's KAHN2 data
    
    isclean = 1;
    
    if isclean
      clean_str = 'clean';
    else
      clean_str = '';
    end
    
    nsfile = dir(fullfile(dataroot,nsDir,[subject,sprintf('%s%s%s.',session{ses},cell2mat(eventValue),clean_str),nsFileExt]));
    if isempty(nsfile)
      error('Cannot find %s file in %s',[subject,sprintf('%s%s%s.',session{ses},cell2mat(eventValue),clean_str),nsFileExt],fullfile(dataroot,nsDir));
    elseif length(nsfile) > 1
      error('More than one %s file found in %s',[subject,sprintf('%s%s%s.',session{ses},cell2mat(eventValue),clean_str),nsFileExt],fullfile(dataroot,nsDir));
    elseif length(nsfile) == 1
      infile_ns = fullfile(dataroot,nsDir,nsfile.name);
    end
  end
  
  % % debug
  % hdr = ft_read_header(infile_ns,'headerformat',ftype);
  % data = ft_read_data(infile_ns,'dataformat',ftype);
  % event = ft_read_event(infile_ns,'eventformat',ftype);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initial parameters for reading the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.dataset = infile_ns;
  cfg.headerfile = infile_ns;
  if ~isempty(ftype)
    cfg.dataformat = ftype;
    cfg.headerformat = ftype;
  end
  cfg.continuous = 'no';
  
  % % debug
  % data = preprocessing(cfg);
  
  %% Select events
  
  % % debug
  % cfg = [];
  % cfg.dataset = infile_ns;
  % cfg.trialdef.eventtype = '?';
  % allEv = ft_definetrial(cfg);
  
  % find out which events are in infile_ns and throw an error if eventValue
  % is not one of these
  cfg_noEv = [];
  cfg_noEv.dataset = infile_ns;
  cfg_noEv.trialdef.eventtype = '?';
  allEv = ft_definetrial(cfg_noEv);
  evVals = cell(size(allEv.event));
  for i = 1:length(allEv.event)
    evVals{i} = allEv.event(i).value;
  end
  evVals = unique(evVals);
  if ~ismember(eventValue,evVals)
    fprintf('The available event values in %s are: %s\n',infile_ns,sprintf(repmat('''%s'' ',1,length(evVals)),evVals{:}));
    error('%s is not in the EEG file. You should redefine exper.eventValues.',cell2mat(eventValue));
  elseif ismember(eventValue,evVals)
    fprintf('You can safely ignore the warning about ''no trialfun was specified''.\n')
  end
  
  % set up for defining the trials based on file type
  cfg.trialdef.eventvalue = eventValue;
  cfg.trialdef.prestim = abs(prepost(1)); % in seconds; must be positive
  cfg.trialdef.poststim = prepost(2); % in seconds; must be positive
  if strcmpi(nsFileExt,'sbin') || strcmpi(nsFileExt,'raw') || strcmpi(nsFileExt,'egis')
    cfg.trialfun = 'seg_trialfun';
    cfg.trialdef.eventtype = 'trial';
  elseif strcmpi(nsFileExt,'set')
    cfg.trialfun = 'trialfun_general';
    cfg.trialdef.eventtype = 'trigger';
  end
  % define the trials
  try
    cfg = ft_definetrial(cfg);
  catch ME
    % if there were zero trials for this event type
    if strfind(ME.message,'no trials were defined')
      fprintf('No %s events found!\n',cell2mat(eventValue));
    end
    fprintf('Returning an empty dataset for %s. This will save an error file when running the ft_*analysis function.\n',cell2mat(eventValue));
    
    % set an empty cell and return to the calling function
    data.trial = {};
    return
  end
  
  %% check on NS artifacts
  
  if rejArt && rejArt_ns
    % make sure the file with NS artifact info exists
    summaryFile = dir(fullfile(dataroot,session{ses},'ns_bci',[subject,'*.bci']));
    if ~isempty(summaryFile)
      summaryFile = fullfile(dataroot,session{ses},'ns_bci',summaryFile.name);
      
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
      format_str = ['%s%d8%d8%s',repmat('%d8',[1,nChan_summary*2]),'%s'];
      fid = fopen(summaryFile,'r');
      sesSummary = textscan(fid,format_str,'Headerlines',1,'delimiter','\t');
      fclose(fid);
    else
      error('Cannot find %s*.%s file in %s. Use the File Export tool to export Metadata > Segment Information.',subject,'bci',fullfile(dataroot,'ns_bci'));
    end
    
    % select only the good trials for this event
    thisEv = find(ismember(sesSummary{1},eventValue));
    badEv = strcmp(sesSummary{4},'bad');
    %goodEv = strcmp(sesSummary{4},'good');
    %cfg.trials = logical(goodEv(min(thisEv):max(thisEv)));
    
    % remove the trials that have artifacts from the trl matrix
    cfg.trl(logical(badEv(min(thisEv):max(thisEv))),:) = [];
  %elseif rejArt && rejArt_ft_auto
  %  % get the trial definition for automated FT artifact rejection
  %  trl = cfg.trl;
  elseif ~rejArt
    fprintf('Not performing any artifact rejection.\n');
  end
  
  %% Get the data and process it if necessary
  
  % get the actual data
  data = ft_preprocessing(cfg);
  
  % find out how many channels are in the data
  nChan_data = length(data.label);
  
  %% Check on channel information
  
  % check on whether we have the reference channel (we want to have it);
  % assuming that the last channel in the series is the reference channel
  %
  % The channel files included with FieldTrip have 3 "extra" (fiduciary)
  % channels defined, so we need to also check using an extra 3 chans
  % subtracted off
  if (nChan_data == nChan_elecfile - 1) || (nChan_data == nChan_elecfile - 4)
    % one less channel because we're checking to see if the reference
    % channel is missing
    error('This dataset is either not rereferenced or the reference channel was not exported. Go back and rereference or export the reference channel in Net Station before running this script!');
  elseif (nChan_data == nChan_elecfile || nChan_data == nChan_elecfile - 3) && var(data.trial{1}(nChan_data,:)) == 0
    % var=0 means that the final (reference) electrode is flat and this
    % data set has not been (average) rereferenced
    error('This dataset is not rereferenced. Go back and rereference in Net Station before running this script!');
  elseif (nChan_data == nChan_elecfile || nChan_data == nChan_elecfile - 3) && var(data.trial{1}(nChan_data,:)) ~= 0
    % has full number of channels and is already rereferenced (final channel is not flat)
    fprintf('Channels are already rereferenced.\n');
    
    % depending on whether the channel string was capitalized or lowercase
    % in the electrode template, make the data elec label match. This is
    % actually important for how FieldTrip deals with electrode numbers.
    %
    % TODO: We now always want to use capital letters, so this should
    % probably be changed.
    if strcmp(elec.label{ceil(nChan_data/2)}(1),'E')
      isCapital = 1;
    elseif strcmp(elec.label{ceil(nChan_data/2)}(1),'e')
      isCapital = 0;
    else
      %error('There is no ''E'' or ''e'' at the start of the electrode number!')
      warning([mfilename,':electrodeCapitalization'],'There is no ''E'' or ''e'' at the start of the electrode number! Going with uppercase.')
      isCapital = 1;
    end
    
    if isCapital
      % capitalize the E for each electrode, or add it in if it's not there
      for c = 1:nChan_data
        if strcmp(data.label{c}(1),'e')
          data.label{c} = upper(data.label{c});
        elseif ~strcmp(data.label{c}(1),'e') && ~strcmp(data.label{c}(1),'E')
          data.label{c} = ['E' data.label{c}];
        end
      end
    elseif ~isCapital
      % make sure the e for each electrode is lowercase, or add it in if
      % it's not there
      for c = 1:nChan_data
        if strcmp(data.label{c}(1),'E')
          data.label{c} = lower(data.label{c});
        elseif ~strcmp(data.label{c}(1),'e') && ~strcmp(data.label{c}(1),'E')
          data.label{c} = ['e' data.label{c}];
        end
      end
    end
    
    % set the last channel name to 'Cz' if that's what was set in
    % elec.label (e.g., instead of 'E129')
    if strcmp(elec.label{end},'Cz')
      if isCapital
        lastChanStr = sprintf('E%d',nChan_data);
      elseif ~isCapital
        lastChanStr = sprintf('e%d',nChan_data);
      end
      %lastChanStr = 'Cz';
      chanindx = find(strcmpi(data.label,lastChanStr));
      if ~isempty(chanindx)
        % set the label for the reference channel
        %data.label{chanindx} = elec.label{chanindx};
        data.label{chanindx} = elec.label{end};
      end
    end
    
  else
    error('Not sure what to do about rereferencing!');
  end
  
%   %% run FieldTrip's automatic artifact detection on the data
%   
%   if rejArt && rejArt_ft_auto
%     %error('FieldTrip artifact detection has not yet been implemented. Set artifactType to ''ns'' or ''none''.');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % look for jump artifacts
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     cfg = [];
%     cfg.trl = trl;
%     cfg.padding = 0;
%     cfg.continuous = 'no';
%     
%     % cutoff and padding
%     % select a set of channels on which to run the artifact detection
%     cfg.artfctdef.zvalue.channel= 'all';
%     cfg.artfctdef.zvalue.cutoff= 30;
%     cfg.artfctdef.zvalue.trlpadding= 0.5*cfg.padding;
%     cfg.artfctdef.zvalue.artpadding= 0.5*cfg.padding;
%     cfg.artfctdef.zvalue.fltpadding= 0;
%     
%     % algorithmic parameters
%     cfg.artfctdef.zvalue.cumulative= 'yes';
%     cfg.artfctdef.zvalue.medianfilter= 'yes';
%     cfg.artfctdef.zvalue.medianfiltord= 9;
%     cfg.artfctdef.zvalue.absdiff= 'yes';
%     
%     % feedback (artifact viewer)
%     cfg.artfctdef.zvalue.feedback= 'yes';
%     
%     [cfg,artifact_jump] = ft_artifact_zvalue(cfg,data);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % look for muscle artifacts
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     cfg = [];
%     cfg.trl = trl;
%     cfg.padding = 0;
%     cfg.continuous = 'no';
%     
%     % cutoff and padding
%     % select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
%     cfg.artfctdef.zvalue.channel = 'all';
%     cfg.artfctdef.zvalue.cutoff      = 4;
%     cfg.artfctdef.zvalue.trlpadding  = 0.1*cfg.padding;
%     cfg.artfctdef.zvalue.fltpadding  = 0.1*cfg.padding;
%     cfg.artfctdef.zvalue.artpadding  = 0.1*cfg.padding;
%     
%     % algorithmic parameters
%     cfg.artfctdef.zvalue.bpfilter    = 'yes';
%     cfg.artfctdef.zvalue.bpfreq      = [110 124];
%     cfg.artfctdef.zvalue.bpfiltord   = 9;
%     cfg.artfctdef.zvalue.bpfilttype  = 'but';
%     cfg.artfctdef.zvalue.hilbert     = 'yes';
%     cfg.artfctdef.zvalue.boxcar      = 0.2;
%     
%     % feedback
%     cfg.artfctdef.zvalue.feedback = 'yes';
%     
%     [cfg,artifact_muscle] = ft_artifact_zvalue(cfg,data);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % look for EOG artifacts
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     cfg = [];
%     cfg.trl = trl;
%     cfg.padding = 0;
%     cfg.continuous = 'no';
%     
%     % cutoff and padding
%     % select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
%     cfg.artfctdef.zvalue.channel = 'all';
%     cfg.artfctdef.zvalue.channel = {'E127','E126','E128','E125'};
%     cfg.artfctdef.zvalue.cutoff      = 6;
%     cfg.artfctdef.zvalue.trlpadding  = 0.5*cfg.padding;
%     cfg.artfctdef.zvalue.artpadding  = 0.1*cfg.padding;
%     cfg.artfctdef.zvalue.fltpadding  = 0.1*cfg.padding;
%     
%     % algorithmic parameters
%     cfg.artfctdef.zvalue.bpfilter   = 'yes';
%     cfg.artfctdef.zvalue.bpfilttype = 'but';
%     cfg.artfctdef.zvalue.bpfreq     = [1 15];
%     cfg.artfctdef.zvalue.bpfiltord  = 4;
%     cfg.artfctdef.zvalue.hilbert    = 'yes';
%     
%     % feedback
%     cfg.artfctdef.zvalue.feedback = 'yes';
%     
%     [cfg,artifact_EOG] = ft_artifact_zvalue(cfg,data);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % reject the automatically defined artifacts
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     cfg=[];
%     cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
%     cfg.artfctdef.jump.artifact = artifact_jump;
%     %cfg.artfctdef.muscle.artifact = artifact_muscle;
%     cfg.artfctdef.eog.artifact =artifact_EOG;
%     data = ft_rejectartifact(cfg,data);
%   end
  
  %% visual artifact inspection (manual)
  
  if rejArt && rejArt_ft_man
    % use cursor drag and click to mark artifacts;
    % use arrows to advance to next trial;
    % use the q key to quit the data browser
    
    fprintf('Processing %s...\n',sprintf(repmat('''%s'' ',1,length(eventValue)),eventValue{:}));
    fprintf('\n\nManual artifact rejection:\n');
    fprintf('Drag mouse to select artifact area; click area to mark an artifact.\n');
    fprintf('Use arrows to move to next trial.\n');
    fprintf('Use the ''i'' key and mouse to identify channels in the data browser.\n');
    fprintf('Use the ''q'' key to quit the data browser when finished.\n\n\n');
    
    cfg = [];
    cfg.continuous = 'no';
    cfg.viewmode = 'butterfly';
    %cfg.viewmode = 'vertical';
    cfg = ft_databrowser(cfg,data);
    
    % reject the artifacts (complete or parial rejection)
    cfg.artfctdef.reject = 'complete';
    data = ft_rejectartifact(cfg,data);
  end
  
  % % HCGSN 129 Eye channels
  % EOGV_upper = [25 8]; % left, right
  % EOGV_lower = [127 126]; % left, right
  % eog = {[25 127], [8 126]}; % left, right
  % EOGH_left = 128;
  % EOGH_right = 125;
  
  %% ICA artifact detection
  
  if rejArt && rejArt_ft_ica
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ICA artifact detection
    %
    % look for eye blink, heart beat, and other artifacts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfg = [];
    cfg.channel = 'all';
    
    data_ic = ft_componentanalysis(cfg,data);
    
    % % OLD METHOD - view the first 20 components
    % cfg = [];
    % cfg.component = 1:20;       % specify the component(s) that should be plotted
    % cfg.layout = elecfile; % specify the layout file that should be used for plotting
    % cfg.comment = 'no';
    % ft_topoplotIC(cfg,data_ic);

    % view the first 20 components (aka channels)
    nComponents = 20;
    cfg = [];
    cfg.viewmode = 'component';
    %cfg.continuous = 'no';
    cfg.continuous = 'yes';
    cfg.blocksize = 10;
    cfg.channels = 1:nComponents;
    cfg.layout = elecfile;
    ft_databrowser(cfg,data_ic);
    
    fprintf('Processing %s...\n',sprintf(repmat('''%s'' ',1,length(eventValue)),eventValue{:}));
    fprintf('\n\nViewing the first %d components.\n',nComponents);
    fprintf('\nLook for patterns that are indicative of artifacts.\n');
    % prompt the user for the component numbers to reject
    componentsToReject = input('\n\nType component numbers to reject and press ''enter'', even if this line moves up due to browsing components (e.g., ''1, 4, 11''):\n\n','s');
    
    % reject the bad components
    cfg = [];
    cfg.component = str2num(componentsToReject);
    data_ic_cleaned = ft_rejectcomponent(cfg,data_ic);
    
    % another manual search of the data for artifacts
    fprintf('Processing %s...\n',sprintf(repmat('''%s'' ',1,length(eventValue)),eventValue{:}));
    fprintf('\n\nManual artifact rejection:\n');
    fprintf('Drag mouse to select artifact area; click area to mark an artifact.\n');
    fprintf('Use arrows to move to next trial.\n');
    fprintf('Use the ''i'' key and mouse to identify channels in the data browser.\n');
    fprintf('Use the ''q'' key to quit the data browser when finished.\n\n\n');
    cfg = [];
    cfg.viewmode = 'butterfly';
    %cfg.viewmode = 'vertical';
    cfg.continuous = 'no';
    cfg = ft_databrowser(cfg,data_ic_cleaned);
    
    % and reject
    cfg.artfctdef.remove = 'complete';
    data = ft_rejectartifact(cfg,data_ic_cleaned);
  end
  
  %% if we're combining multiple sessions, add the data to the append struct
  if length(session) > 1
    append_data.(ses_str) = data;
  end
end

%% Append sessions, if necessary
  
% run ft_appenddata if we're combining multiple sessions
if length(session) > 1
  ses_str = sprintf('ses_%s',session{1});
  append_str = sprintf('append_data.%s',ses_str);
  
  for ses = 2:length(session)
    ses_str = sprintf('ses_%s',session{ses});
    append_str = cat(2,append_str,sprintf(',append_data.%s',ses_str));
  end
  
  data = eval(sprintf('ft_appenddata([],%s);',append_str));
end

%% Separate the event values

% initialize the struct to return
ft_raw = struct;

if length(eventValue) > 1
  for evVal = 1:length(eventValue)
    
    cfg = [];
    trl = ft_findcfg(data.cfg,'trl');
    cfg.trl = trl(trl(:,4) == evVal,1:3);
    
    if ~isempty(cfg.trl)
      fprintf('Selecting %d trials for %s...\n',size(cfg.trl,1),eventValue{evVal});
      ft_raw.(eventValue{evVal}) = ft_redefinetrial(cfg,data);
      fprintf('Done.\n');
    else
      fprintf('No trials found for %s!\n',eventValue{evVal});
      ft_raw.(eventValue{evVal}).trial = {};
      %keyboard
    end
  end
elseif length(eventValue) == 1
  ft_raw.(eventValue) = data;
end


end
