function [data,elec] = seg2ft_kahn2(dataroot,nsFileExt,subject,session,eventValue,prepost,elecfile)

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
  % make sure the EEG file exists
  isclean = 1;
  if isclean
    clean_str = 'clean';
  else
    clean_str = '';
  end
  
  nsfile = dir(fullfile(dataroot,nsDir,[subject,sprintf('%s%s%s.',session{ses},cell2mat(eventValue),clean_str),nsFileExt]));
  if isempty(nsfile)
    error('Cannot find % file in %s',[subject,sprintf('%s%s%s.',session{ses},cell2mat(eventValue),clean_str),nsFileExt],fullfile(dataroot,nsDir));
  elseif length(nsfile) > 1
    error('More than one %s file found in %s',[subject,sprintf('%s%s%s.',session{ses},cell2mat(eventValue),clean_str),nsFileExt],fullfile(dataroot,nsDir));
  elseif length(nsfile) == 1
    infile_ns = fullfile(dataroot,nsDir,nsfile.name);
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
  % ft_definetrial(cfg);
  
  % define the trials
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
  
  % get the actual data
  data = ft_preprocessing(cfg);
  
  % check on whether we have the reference channel (we want to have it);
  % assuming that the last channel in the series is the reference channel
  if (length(data.label) == nChan - 1) || (length(data.label) == nChan - 4)
    error('This dataset is either not rereferenced or the reference channel was not exported. Go back and rereference or export the reference channel in Net Station before running this script!');
  elseif (length(data.label) == nChan || length(data.label) == nChan - 3) && var(data.trial{1}(nChan_summary,:)) == 0
    error('This dataset is not rereferenced. Go back and rereference in Net Station before running this script!');
  elseif (length(data.label) == nChan || length(data.label) == nChan - 3) && var(data.trial{1}(nChan_summary,:)) ~= 0
    fprintf('Channels are already rereferenced.\n');
    % has full number of channels; already rereferenced (E129 is not flat)
    
    % last channel name
    lastChanStr = sprintf('E%d',nChan);
    %lastChanStr = 'Cz';
    chanindx = find(strcmpi(data.label,lastChanStr));
    if ~isempty(chanindx)
      % set the label for the reference channel
      data.label{chanindx} = elec.label{chanindx};
    end
  else
    error('Not sure what to do about rereferencing!');
  end
  
  % add the data to the append struct if we're combining multiple sessions
  if length(session) > 1
    append_data.(session{ses}) = data;
  end
end

% run ft_appenddata if we're combining multiple sessions
if length(session) > 1
  append_str = [];
  for ses = 1:length(session)
    append_str = cat(2,append_str,sprintf(',append_data.%s',session{ses}));
  end
  
  data = eval(sprintf('ft_appenddata([]%s);',append_str));
end
