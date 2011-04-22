function [data] = seg2ft(dataroot,nsFileExt,subject,session,eventValue,prepost,elecfile)
%SEG2FT: take segmented Net Station data and put it in FieldTrip format
%
% [data] = seg2ft(dataroot,nsFileExt,subject,session,eventValue,prepost,elecfile)
%
% The Net Station file should be exported as either EGIS (which has 129
% channels, where the final channel [Cz] is the reference; extension:
% 'egis') or NS Simple Binary (extension: 'raw' or 'sbin') calibrated,
% including the reference channel. These options are in the file export
% tool. create_ft_struct, the function that calls seg2ft, expects EGIS
% files to be stored in a 'ns_egis' directory at the level of dirs.dataDir.
% If using raw files, they should be in 'ns_raw' instead.
%
% This function can deal with event values that have zero events. It's
% probably better to check on your event count and exclude those subjects
% with zero events for any of the event values you're trying to keep before
% trying to run this script, but we will insert an empty event entry for
% that eventValue, and the subjects should get excluded when using
% mm_threshSubs.m.
%
% Currently this function searches for an eye artifact information file.
% This is exported from Net Station using the File Export tool. To set up
% the tool, export format is metadata and check the segment information
% option (segment information is output with a .bci extension). Run it on
% the file that was exported to egis/raw (i.e., the baseline correction
% file or the average rereference file). The bci files should be stored in
% a 'ns_bci' directory at the same level as 'ns_egis' or 'ns_raw'.
%
% If using eeglab data, no artifact detection is done and no bci file is
% expected to exist. The directory structure is different and can be
% gleaned by examining the code here, but right now it is only set up to
% process Erika Nyhus's KAHN2 data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things to add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TODO: add a case for when we want to use NS bci artifact information
% versus no artifact information
%
% TODO: add support for session names in the EEG file path; currently it
% ignores the session string
%


% make sure eventValue is set up correctly
if ~iscell(eventValue)
  eventValue = {eventValue};
end
if length(eventValue) > 1
  error('Expecting eventValue to have length of one.');
end
% for i = 1:length(eventValue)
%   if length(eventValue{i}) > 4
%     error('Each entry in eventValue must be <= 4 characters!');
%   end
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
nChan = size(elec.label,1);

for ses = 1:length(session)
  ses_str = sprintf('ses%s',session{ses});
  
  if strcmpi(nsFileExt,'sbin') || strcmpi(nsFileExt,'raw') || strcmpi(nsFileExt,'egis')
    % make sure the EEG file exists
    nsfile = dir(fullfile(dataroot,nsDir,[subject,'*.',nsFileExt]));
    if isempty(nsfile)
      error('Cannot find %s*.%s file in %s',subject,nsFileExt,fullfile(dataroot,nsDir));
    elseif length(nsfile) > 1
      error('More than one %s*.%s file found in %s',subject,nsFileExt,fullfile(dataroot,nsDir));
    elseif length(nsfile) == 1
      infile_ns = fullfile(dataroot,nsDir,nsfile.name);
    end
    
    % make sure the file with NS artifact info exists
    summaryFile = dir(fullfile(dataroot,'ns_bci',[subject,'*.bci']));
    if ~isempty(summaryFile)
      summaryFile = fullfile(dataroot,'ns_bci',summaryFile.name);
      
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
        nChan_summary = nChan;
      end
      
      % read in the artifact file
      format_str = ['%s%d8%d8%s',repmat('%d8',[1,nChan_summary*2]),'%s'];
      fid = fopen(summaryFile,'r');
      sesSummary = textscan(fid,format_str,'Headerlines',1,'delimiter','\t');
      fclose(fid);
    else
      error('Cannot find %s*.%s file in %s. Use the File Export tool to export Metadata > Segment Information.',subject,'bci',fullfile(dataroot,'ns_bci'));
    end
    
  elseif strcmpi(nsFileExt,'set')
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
  % Step 0: Initial parameters
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
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Step 1: Select events
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % % debug
  % cfg = [];
  % cfg.dataset = infile_ns;
  % cfg.trialdef.eventtype = '?';
  % allev = ft_definetrial(cfg);
  
  % define the trials
  if strcmpi(nsFileExt,'sbin') || strcmpi(nsFileExt,'raw') || strcmpi(nsFileExt,'egis')
    cfg.trialdef.prestim = abs(prepost(1)); % in seconds; must be positive
    cfg.trialdef.poststim = prepost(2); % in seconds; must be positive
    cfg.trialfun = 'seg_trialfun';
    cfg.trialdef.eventtype = 'trial';
    cfg.trialdef.eventvalue = eventValue;
    try
      cfg = ft_definetrial(cfg);
    catch ME
      % if there were zero trials for this event type
      if strfind(ME.message,'no trials were defined')
        fprintf('No %s events found!\n',eventValue{:});
      end
      data.trial = {};
      return
    end
    
    % select only the good trials for this event
    %cfg = [];
    thisEv = find(strcmp(eventValue,sesSummary{1}));
    badEv = strcmp('bad',sesSummary{4});
    %goodEv = strcmp('good',sesSummary{4});
    %cfg.trials = logical(goodEv(min(thisEv):max(thisEv)));
    
    % remove the trials that have artifacts from the trl
    cfg.trl(logical(badEv(min(thisEv):max(thisEv))),:) = [];
  elseif strcmpi(nsFileExt,'set')
    cfg.trialdef.prestim = abs(prepost(1));
    cfg.trialdef.poststim = prepost(2);
    cfg.trialfun = 'trialfun_general';
    cfg.trialdef.eventtype = 'trigger';
    cfg.trialdef.eventvalue = eventValue;
    try
      cfg = ft_definetrial(cfg);
    catch ME
      % if there were zero trials for this event type
      if strfind(ME.message,'no trials were defined')
        fprintf('No %s events found!\n',eventValue{:});
      end
      data.trial = {};
      return
    end
  end
  
  % get the actual data
  data = ft_preprocessing(cfg);
  
  if strcmpi(nsFileExt,'set')
    nChan_summary = length(data.label);
  end
  
  % check on whether we have the reference channel (we want to have it);
  % assuming that the last channel in the series is the reference channel
  %
  % The channel files included with FieldTrip have 3 "extra" channels
  % defined, so an extra 3 needs to be subtracted off when checking
  if (length(data.label) == nChan - 1) || (length(data.label) == nChan - 4)
    error('This dataset is either not rereferenced or the reference channel was not exported. Go back and rereference or export the reference channel in Net Station before running this script!');
  elseif (length(data.label) == nChan || length(data.label) == nChan - 3) && var(data.trial{1}(nChan_summary,:)) == 0
    error('This dataset is not rereferenced. Go back and rereference in Net Station before running this script!');
  elseif (length(data.label) == nChan || length(data.label) == nChan - 3) && var(data.trial{1}(nChan_summary,:)) ~= 0
    fprintf('Channels are already rereferenced.\n');
    % has full number of channels; already rereferenced (E129 is not flat)
    
    % depending on whether the channel string was capitalized or lowercase
    % in the electrode template, make the data elec label match. This is
    % actually important for how fieldtrip deals with electrode numbers. We
    % want to use capital letters, so maybe this should be changed. FIXME.
    if strcmp(elec.label{ceil(nChan_summary/2)}(1),'E')
      isCapital = 1;
    elseif strcmp(elec.label{ceil(nChan_summary/2)}(1),'e')
      isCapital = 0;
    else
      %error('There is no ''E'' or ''e'' at the start of the electrode number!')
      warning([mfilename,':electrodeCapitalization'],'There is no ''E'' or ''e'' at the start of the electrode number! Going with uppercase.')
      isCapital = 1;
    end
    
    if isCapital
      % capitalize the E for each electrode, or add it in if it's not there
      for c = 1:length(data.label)
        if strcmp(data.label{c}(1),'e')
          data.label{c}(1) = 'E';
        elseif ~strcmp(data.label{c}(1),'e') && ~strcmp(data.label{c}(1),'E')
          data.label{c} = ['E' data.label{c}];
        end
      end
    elseif ~isCapital
      % make sure the e for each electrode is lowercase, or add it in if
      % it's not there
      for c = 1:length(data.label)
        if strcmp(data.label{c}(1),'E')
          data.label{c}(1) = 'e';
        elseif ~strcmp(data.label{c}(1),'e') && ~strcmp(data.label{c}(1),'E')
          data.label{c} = ['e' data.label{c}];
        end
      end
    end
    
    % set the last channel name to Cz instead of E129
    if strcmp(elec.label{end},'Cz')
      if isCapital
        lastChanStr = sprintf('E%d',nChan_summary);
      elseif ~isCapital
        lastChanStr = sprintf('e%d',nChan_summary);
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
  
  % add the data to the append struct if we're combining multiple sessions
  if length(session) > 1
    append_data.(ses_str) = data;
  end
end

% run ft_appenddata if we're combining multiple sessions
if length(session) > 1
  ses_str = sprintf('ses%s',session{1});
  append_str = sprintf('append_data.%s',ses_str);
  
  for ses = 2:length(session)
    ses_str = sprintf('ses%s',session{ses});
    append_str = cat(2,append_str,sprintf(',append_data.%s',ses_str));
  end
  
  data = eval(sprintf('ft_appenddata([],%s);',append_str));
end
