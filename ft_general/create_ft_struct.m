function [exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files)
%CREATE_FT_STRUCT: Create a FieldTrip data structure from NS files
%
% [exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files)
%
% The default function (set in ana.segFxn) that this function calls is
% SEG2FT. See 'HELP SEG2FT' for more information.
%
% For multiple sessions, SEG2FT will collapse across them if they're input
% as {{'ses1','ses2'}}. If they're input as {{'ses1'},{'ses2'}} then they
% will be kept separate.
%
% ana.artifact.type = the type of artifact info to use: 'none', 'nsAuto',
%                     'zeroVar', 'preRejManual','ftManual', or 'ftICA'.
%
%                     Some types can be combined.
%                     
%                     If nsAuto, SEG2FT expects to find a NS bci metadata
%                     file. This file denotes the reject artifact trials.
%                     Use the File Export tool to export Metadata > Segment
%                     Information; put it in a 'ns_bci' directory in
%                     dataroot/session.
%                     
%                     See SEG2FT for more information.
%                     (default = 'none')
%
% ana.overwrite.raw = prevent overwriting of raw data, and load in existing
%                     raw data (any otherFxns will be run). Binary.
%                     Default: 1
%                     DO NOT USE YET, NOT FULLY IMPLEMENTED. (TODO)
%
% TODO: when running SEG2FT, choice for processing all event values at once
% or one at a time (mostly useful for manual artifact rejection)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on equating trial counts when also creating new extra event values:
%
% exper.eventValuesExtra.equateExtrasSeparately = equate extras based on
% their own minimum event number so you can still have bigger numbers with
% the normal events
%
% if exper.equateTrials==1 and exper.eventValuesExtra.onlyKeepExtras==1,
% the trial counts of the unkept event values will not necessarily add up
% to those of the kept event values; this is due to the way unkept and kept
% event values are processed.
%
% if exper.equateTrials==1 and event values are being combined to form
% exper.eventValuesExtra.newValue, there will not necessarily be equal
% contibution from each exper.eventValue to the new values because this
% could cut down on the number of trials if one event value had a low trial
% count. Instead, the event values are combined and them randomly sampled
% to create equal bin sizes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also: SEG2FT, PROCESS_FT_DATA, MM_FT_ARTIFACT
%

%% make sure some required fields are set

% overwrite by default
if ~isfield(ana,'overwrite')
  ana.overwrite.raw = 1;
elseif isfield(ana,'overwrite') && ~isfield(ana.overwrite,'raw')
  ana.overwrite.raw = 1;
elseif isfield(ana,'overwrite') && isfield(ana.overwrite,'raw') && ana.overwrite.raw == 0
  % TODO
  error('Prevention of overwriting is not fully implemented. You must set overwrite to 1.');
end

% make sure event values are sorted
if ~issorted(exper.eventValues)
  exper.eventValues = sort(exper.eventValues);
end

% make sure exper.sessions is a cell
if ~iscell(exper.sessions)
  exper.sessions = {exper.sessions};
end

% need a segmentation function
if ~isfield(ana,'segFxn')
  ana.segFxn = 'seg2ft';
end

% need an artifact detection type ('none', or: 'nsAuto', 'preRejManual', 'ftManual', 'ftICA')
if ~isfield(ana,'artifact') || (isfield(ana,'artifact') && ~isfield(ana.artifact,'type'))
  ana.artifact.type = {'none'};
elseif isfield(ana,'artifact') && isfield(ana.artifact,'type') && ischar(ana.artifact.type)
  ana.artifact.type = {ana.artifact.type};
end

% check on any extra events
if ~isfield(exper,'eventValuesExtra')
  exper.eventValuesExtra = {};
end
if ~isfield(exper.eventValuesExtra,'newValue')
  exper.eventValuesExtra.newValue = {};
elseif isfield(exper.eventValuesExtra,'newValue') && ~isempty(exper.eventValuesExtra.newValue) && ~iscell(exper.eventValuesExtra.newValue{1})
  exper.eventValuesExtra.newValue = {exper.eventValuesExtra.newValue};
end
if ~isfield(exper.eventValuesExtra,'toCombine')
  exper.eventValuesExtra.toCombine = {};
elseif isfield(exper.eventValuesExtra,'toCombine') && ~isempty(exper.eventValuesExtra.toCombine) && ~iscell(exper.eventValuesExtra.toCombine{1})
  exper.eventValuesExtra.toCombine = {exper.eventValuesExtra.toCombine};
end
if ~isfield(exper.eventValuesExtra,'onlyKeepExtras')
  exper.eventValuesExtra.onlyKeepExtras = 0;
end

% check on event equating parameters
if ~isfield(exper,'equateTrials')
  exper.equateTrials = 0;
end

% equate the extra events based on their own minimum number?
if ~isfield(exper.eventValuesExtra,'equateExtrasSeparately')
  exper.eventValuesExtra.equateExtrasSeparately = 0;
end
% can't equate the extra event values separately if trials aren't being
% equated in the first place; this could be done if the code was modified,
% but it doesn't make sense to.
if exper.eventValuesExtra.equateExtrasSeparately == 1 && exper.equateTrials == 0
  warning([mfilename,':equateExtrasSeparately'],'Cannot equate the extra event values separately if trials are not being equated in the first place. Setting equateExtrasSeparately=0\n');
  exper.eventValuesExtra.equateExtrasSeparately = 0;
end

% check on the other functions
if isfield(ana,'otherFxn')
  if ~iscell(ana.otherFxn)
    ana.otherFxn = {ana.otherFxn};
  end
  if ~isfield(ana,'cfg_other')
    error('If ana.otherFxn is set, ana.cfg_other must be defined.');
  else
    for i = 1:length(ana.cfg_other)
      if ~isfield(ana.cfg_other{i},'ftype')
        error('Need to set a file type in ana.cfg_other{%d}.ftype to name the outputfile (e.g., ''tla'', ''pow'', ''powandcsd'', ''fourier'')',i);
      end
    end
  end
end

% make sure we have the right electrode information
if ~isfield(files,'elecfile')
  error('Must set files.elecfile');
end
if ~isempty(strfind(files.elecfile,'GSN'))
  ana.isNS = true;
else
  ana.isNS = false;
end
if ~isfield(exper,'refChan')
  if ana.isNS
    warning([mfilename,':refChanDefault'],'Setting exper.refChan to average rereference default: ''Cz''.');
    exper.refChan = 'Cz';
  else
    error('Must set exper.refChan (e.g., ''Cz'' or channel number)');
  end
end

%% get all the exper.eventValues and exper.eventValuesExtra.newValue together; make sure the extra event values aren't already in the list

if ~isempty(exper.eventValuesExtra.newValue)
  eventValuesWithExtra = exper.eventValues;
  for nVal = 1:length(exper.eventValuesExtra.newValue)
    if ~ismember(exper.eventValuesExtra.newValue{nVal},exper.eventValues)
      eventValuesWithExtra = cat(2,eventValuesWithExtra,exper.eventValuesExtra.newValue{nVal});
    else
      error('%s is already in the event value list!',exper.eventValuesExtra.newValue{nVal}{1});
    end
  end
end

%% initialize for storing some information

% store the number of events for each subject and session all of the event
% values regardless of whether we're keeping them
if ~isempty(exper.eventValuesExtra.newValue)
  for evVal = 1:length(eventValuesWithExtra)
    exper.nTrials.(eventValuesWithExtra{evVal}) = zeros(length(exper.subjects),length(exper.sessions));
  end
else
  for evVal = 1:length(exper.eventValues)
    exper.nTrials.(exper.eventValues{evVal}) = zeros(length(exper.subjects),length(exper.sessions));
  end
end

% store the bad channel information
exper.badChan = cell(length(exper.subjects),length(exper.sessions));

%% Store the raw data in FieldTrip format

fprintf('Converting NetStation to FieldTrip...\n');

if exper.equateTrials && length(exper.eventValues) > 1
  % if we're going to equate trials across conditions, need to use a random
  % selection
  rng('shuffle','twister');
elseif exper.equateTrials && length(exper.eventValues) == 1
  % don't equate
  fprintf('Not equating trials because there is only 1 event.\n');
  exper.equateTrials = 0;
  exper.eventValuesExtra.equateExtrasSeparately = 0;
end

% initialize for the sesStr
exper.sesStr = cell(size(exper.sessions));

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    
    % turn the session name into a string
    if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
      sesStr = exper.sessions{ses}{1};
      for i = 2:length(exper.sessions{ses})
        sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
      end
    elseif ~iscell(exper.sessions{ses}) || (iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1)
      sesStr = exper.sessions{ses};
    end
    % store the sesStr
    exper.sesStr{ses} = sesStr;
    
    % set the location to save the data and make sure it exists
    saveDirRawFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
    if ~exist(saveDirRawFile,'dir')
      mkdir(saveDirRawFile);
    end
    
    % if we're going to create new values out of the events...
    if ~isempty(exper.eventValuesExtra.newValue)
      % initialize the append_data structure
      append_data = struct;
    end
    
    % if we don't want to overwrite any files...
    if ~ana.overwrite.raw
      % see if any raw files already exist; run only the missing ones
      rawFiles = dir(fullfile(saveDirRawFile,'data_raw_*.mat'));
      if ~isempty(rawFiles)
        rawEvents = cell(size(rawFiles))';
        for rf = 1:length(rawFiles)
          rawEvents{rf} = strrep(strrep(rawFiles(rf).name,'data_raw_',''),'.mat','');
        end
        
        if ~exper.eventValuesExtra.onlyKeepExtras
          if isempty(exper.eventValuesExtra.newValue)
            eventValuesToProcess = exper.eventValues(~ismember(exper.eventValues,rawEvents));
          elseif ~isempty(exper.eventValuesExtra.newValue)
            eventValuesToProcess = eventValuesWithExtra(~ismember(eventValuesWithExtra,rawEvents));
          end
        elseif exper.eventValuesExtra.onlyKeepExtras && ~isempty(exper.eventValuesExtra.newValue)
          eventValuesToProcess = exper.eventValuesExtra.toCombine(~ismember(cat(2,exper.eventValuesExtra.newValue{:}),rawEvents));
          eventValuesToProcess = cat(2,eventValuesToProcess{:});
        else
          error('eventValuesToProcess was not set correctly.');
        end
      else
        eventValuesToProcess = exper.eventValues;
      end
    else
      eventValuesToProcess = exper.eventValues;
    end
    
    if isempty(eventValuesToProcess)
      fprintf('All event values already processed. Loading raw files.\n');
      
      ft_raw = struct;
      
      if ~ana.overwrite.raw
        % load in the ones we didn't process
        if ~isempty(rawFiles)
          for rf = 1:length(rawFiles)
            ft_raw.(rawEvents{rf}) = load(fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',rawEvents{rf})));
          end
        end
      end
    elseif ~strcmpi(exper.eegFileExt,'set') && ~isempty(eventValuesToProcess)
      fprintf('Creating FT struct of raw EEG data: %s, %s%s.\n',exper.subjects{sub},sesStr,sprintf(repmat(', ''%s''',1,length(eventValuesToProcess)),eventValuesToProcess{:}));
      
      % collect all the raw data
      [ft_raw,badChan] = feval(str2func(ana.segFxn),fullfile(dirs.dataroot,dirs.dataDir),exper.subjects{sub},exper.sessions{ses},eventValuesToProcess,files.elecfile,ana,exper);
      
      if ~ana.overwrite.raw
        % load in the ones we didn't process
        if ~isempty(rawFiles)
          for rf = 1:length(rawFiles)
            ft_raw.(rawEvents{rf}) = load(fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',rawEvents{rf})));
          end
        end
      end
      
      % if an extra value was defined store each of its sub-values in the
      % toCombine field for appending them together at a later point
      if ~isempty(exper.eventValuesExtra.newValue)
        % get the name of this event type
        for evVal = 1:length(exper.eventValues)
          eventVal = exper.eventValues{evVal};
          for cVal = 1:length(exper.eventValuesExtra.toCombine)
            if ismember(eventVal,cat(2,exper.eventValuesExtra.toCombine{cVal})) && ~isempty(ft_raw.(eventVal).trial)
              append_data.(cell2mat(exper.eventValuesExtra.newValue{cVal})).(eventVal) = ft_raw.(eventVal);
            end
          end
        end
      end
      
    elseif strcmp(exper.eegFileExt,'set') && ~isempty(eventValuesToProcess)
      % if we're processing EEGLAB data; I don't remember why this is
      % separate
      
      % initialize the data structure
      ft_raw = struct;
      
      if ~ana.overwrite.raw
        % load in the ones we didn't process
        if ~isempty(rawFiles)
          for rf = 1:length(rawFiles)
            ft_raw.(rawEvents{rf}) = load(fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',rawEvents{rf})));
          end
        end
      end
      
      % collect all the raw data
      for evVal = 1:length(eventValuesToProcess)
        % get the name of this event type
        eventVal = eventValuesToProcess{evVal};
        
        fprintf('Creating FT struct of raw EEG data: %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
        
        % turn the NS segments into raw FT data; uses the NS bci metadata
        % file to reject artifact trials before returning good trials in
        % ft_raw
        [ft_raw.(eventVal),badChan] = feval(str2func(ana.segFxn),fullfile(dirs.dataroot,dirs.dataDir),exper.subjects{sub},exper.sessions{ses},eventValuesToProcess(evVal),files.elecfile,ana,exper);
        
        % if an extra value was defined and the current event is one of its
        % sub-values, store the current event in the toCombine field for
        % appending the sub-values together at a later point
        if ~isempty(exper.eventValuesExtra.newValue)
          for cVal = 1:length(exper.eventValuesExtra.toCombine)
            if ismember(eventVal,cat(2,exper.eventValuesExtra.toCombine{cVal})) && ~isempty(ft_raw.(eventVal).trial)
              append_data.(cell2mat(exper.eventValuesExtra.newValue{cVal})).(eventVal) = ft_raw.(eventVal);
            end
          end
        end
      end % for evVal
    end % if
    
    % store the bad channel information
    exper.badChan{sub,ses} = badChan;
    
    % check on any empty events
    for evVal = 1:length(eventValuesToProcess)
      % get the name of this event type
      eventVal = eventValuesToProcess{evVal};
      if ~isempty(ft_raw.(eventVal).trial)
        % store the number of trials for this event value
        exper.nTrials.(eventVal)(sub,ses) = size(ft_raw.(eventVal).trial,2);
        
        fprintf('Created FT struct of raw EEG data: %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
      else
        warning([mfilename,':noTrialsFound'],'Did not create FT struct of raw EEG data: %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
        if ~exist('emptyEvents','var')
          emptyEvents = {eventVal};
        elseif exist('emptyEvents','var')
          emptyEvents = cat(2,emptyEvents,eventVal);
        end
      end
    end
    
    % if there were no trials for an eventValue and it is a sub-value of an
    % extra value, remove that evVal from the toCombine field
    if ~isempty(exper.eventValuesExtra.newValue) && exist('emptyEvents','var')
      orig = struct;
      % keep original eventValuesExtra to restore it for the next subject
      orig.eventValuesExtra = exper.eventValuesExtra;
      for evVal = 1:length(emptyEvents)
        eventVal = emptyEvents{evVal};
        for cVal = 1:length(exper.eventValuesExtra.toCombine)
          if ismember(eventVal,cat(2,exper.eventValuesExtra.toCombine{cVal}))
            exper.eventValuesExtra.toCombine{cVal} = exper.eventValuesExtra.toCombine{cVal}(~ismember(exper.eventValuesExtra.toCombine{cVal},emptyEvents));
          end
        end
      end
    end
    
    % if we want to equate the numbers of events within each subject, find
    % the minimum number of trials across these events. If we're not only
    % keeping the extras, the numbers of events in extraEventValues
    % won't matter because each extra is, at minimum, a combination of
    % these events.
    if exper.equateTrials && ~exper.eventValuesExtra.onlyKeepExtras
      % set an Inf event size
      minNumEv = Inf;
      % check the event sizes
      for evVal = 1:length(exper.eventValues)
        if size(ft_raw.(exper.eventValues{evVal}).trial,2) < minNumEv && ~isempty(ft_raw.(exper.eventValues{evVal}).trial)
          minNumEv = size(ft_raw.(exper.eventValues{evVal}).trial,2);
        end
      end
    end
    
    % if we want to keep the events in exper.eventValues (regardless of
    % whether there were any extras)
    if ~exper.eventValuesExtra.onlyKeepExtras
      for evVal = 1:length(exper.eventValues)
        % get the name of this event type
        eventVal = exper.eventValues{evVal};
        
        if size(ft_raw.(eventVal).trial,2) > 0
          % if we want to equate the numbers of events within each subject,
          % select a random subset
          if exper.equateTrials
            cfg_eq = [];
            fprintf('Randomly selecting %d events (out of %d) for %s, %s, %s.\n',minNumEv,size(ft_raw.(eventVal).trial,2),exper.subjects{sub},sesStr,eventVal);
            % do a random permutation of all of this eventVal's trial numbers
            cfg_eq.evNums = randperm(size(ft_raw.(eventVal).trial,2));
            % select only the number of events that we need
            cfg_eq.evNums = sort(cfg_eq.evNums(1:minNumEv));
            
            % get the raw data's trl for our new configuration
            cfg_eq.trl = ft_raw.(eventVal).cfg.trl;
            % remove the trials were not randomly chosen
            cfg_eq.trl(logical(~ismember(1:size(ft_raw.(eventVal).trial,2),cfg_eq.evNums)),:) = [];
            
            ft_raw.(eventVal) = ft_redefinetrial(cfg_eq,ft_raw.(eventVal));
          end
          
          % store the number of trials for this event value
          exper.nTrials.(eventVal)(sub,ses) = size(ft_raw.(eventVal).trial,2);
          
          % set outputfile so the raw data is saved; especially useful for
          % implementing analyses with the peer toolbox
          cfg_pp.outputfile = fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',eventVal));
          
          % do any preprocessing
          fprintf('Running ft_preprocessing on %s, %s, %s...\n',exper.subjects{sub},sesStr,eventVal);
          ft_preprocessing(cfg_pp,ft_raw.(eventVal));
          fprintf('Done with %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
          
          % any other functions to run? (e.g., ft_resampledata)
          data_suffix = [];
          if isfield(ana,'otherFxn')
            for i = 1:length(ana.otherFxn)
              fprintf('Running %s...\n',ana.otherFxn{i});
              cfg = ana.cfg_other{i};
              % set the file to be read (data_suffix should be [] on 1st
              % time through)
              cfg.inputfile = fullfile(saveDirRawFile,sprintf('data_raw%s_%s.mat',data_suffix,eventVal));
              % run the function
              ft_raw.(eventVal) = feval(str2func(ana.otherFxn{i}),cfg);
              
              % add the new data_suffix for the outputfile
              data_suffix = cat(2,data_suffix,sprintf('_%s',cfg.ftype));
              
              cfg = [];
              % set the file to save
              cfg.outputfile = fullfile(saveDirRawFile,sprintf('data_raw%s_%s.mat',data_suffix,eventVal));
              % preprocess it again
              ft_raw.(eventVal) = ft_preprocessing(cfg,ft_raw.(eventVal));
              fprintf('Done with %s.\n',ana.otherFxn{i});
            end
          end
        else
          fprintf('No good trials found for %s, %s, %s!\n',exper.subjects{sub},sesStr,eventVal);
        end % if there were trials
      end % for exper.eventValues
    end % if ~exper.eventValuesExtra.onlyKeepExtras
    
    % if we have an extra value, append the sub-values together with
    % ft_appenddata
    if ~isempty(exper.eventValuesExtra.newValue)
      for cVal = 1:length(exper.eventValuesExtra.toCombine)
        % create a cell of the events to combine
        evValToCombine = exper.eventValues(ismember(exper.eventValues,cat(2,exper.eventValuesExtra.toCombine{cVal})));
        
        % get the name of the new event
        eventVal = cell2mat(exper.eventValuesExtra.newValue{cVal});
        
        if length(evValToCombine) > 1
          % create the string to pass into ft_appenddata
          append_str = sprintf('append_data.%s.%s',cell2mat(exper.eventValuesExtra.newValue{cVal}),evValToCombine{1});
          append_vals = evValToCombine{1};
          for cInd = 2:length(evValToCombine)
            append_str = cat(2,append_str,sprintf(',append_data.%s.%s',cell2mat(exper.eventValuesExtra.newValue{cVal}),evValToCombine{cInd}));
            append_vals = cat(2,append_vals,'+',evValToCombine{cInd});
          end
          % append the data
          ft_raw.(eventVal) = eval(sprintf('ft_appenddata([],%s);',append_str));
          fprintf('Created appended FT struct: %s, %s, %s (from %s).\n',exper.subjects{sub},sesStr,eventVal,append_vals);
        elseif length(evValToCombine) == 1
          fprintf('Only one event value (%s) was specified to create %s. Not running ft_append_data, setting %s=%s.\n',evValToCombine{1},eventVal,eventVal,evValToCombine{1});
          ft_raw.(eventVal) = ft_raw.(evValToCombine{1});
        elseif isempty(evValToCombine)
          % make sure it's not empty
          fprintf('No event values found to create %s.\n',eventVal);
          ft_raw.(eventVal).trial = {};
          continue
        end
      end % for cVal
      
      % if we want to equate the numbers of events within each subject,
      % find the minimum number of trials across these events. This time we
      % need to take the extras into account because either they alone are
      % being kept, or they're being equated separately.
      if exper.equateTrials && (exper.eventValuesExtra.onlyKeepExtras || exper.eventValuesExtra.equateExtrasSeparately)
        % set an Inf event size
        minNumEv = Inf;
        % check the event sizes
        for evVal = 1:length(exper.eventValuesExtra.toCombine)
          eventVal = cell2mat(exper.eventValuesExtra.newValue{evVal});
          if size(ft_raw.(eventVal).trial,2) < minNumEv && ~isempty(ft_raw.(eventVal).trial)
            minNumEv = size(ft_raw.(eventVal).trial,2);
          end
        end
      end
      
      for cVal = 1:length(exper.eventValuesExtra.toCombine)
        % get the name of the new event
        eventVal = cell2mat(exper.eventValuesExtra.newValue{cVal});
        
        if size(ft_raw.(eventVal).trial,2) > 0
          % if we want to equate the numbers of events within each subject,
          % select a random subset
          if exper.equateTrials
            cfg_eq = [];
            fprintf('Randomly selecting %d events (out of %d) for %s, %s, %s.\n',minNumEv,size(ft_raw.(eventVal).trial,2),exper.subjects{sub},sesStr,eventVal);
            % do a random permutation of all of this eventVal's trial numbers
            cfg_eq.evNums = randperm(size(ft_raw.(eventVal).trial,2));
            % select only the number of events that we need
            cfg_eq.evNums = sort(cfg_eq.evNums(1:minNumEv));
            
            % get the raw data's trl for our new configuration
            cfg_eq.trl = ft_raw.(eventVal).cfg.trl;
            % remove the trials were not randomly chosen
            cfg_eq.trl(logical(~ismember(1:size(ft_raw.(eventVal).trial,2),cfg_eq.evNums)),:) = [];
            
            ft_raw.(eventVal) = ft_redefinetrial(cfg_eq,ft_raw.(eventVal));
          end
          
          % store the number of trials for this event value
          exper.nTrials.(eventVal)(sub,ses) = size(ft_raw.(eventVal).trial,2);
          
          % set outputfile so the raw data is saved; especially useful for
          % implementing analyses with the peer toolbox
          cfg_pp.outputfile = fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',eventVal));
          
          % do any preprocessing
          fprintf('Running ft_preprocessing on %s, %s, %s...\n',exper.subjects{sub},sesStr,eventVal);
          ft_preprocessing(cfg_pp,ft_raw.(eventVal));
          fprintf('Done with %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
          
          % any other functions to run? (e.g., ft_resampledata)
          data_suffix = [];
          if isfield(ana,'otherFxn')
            for i = 1:length(ana.otherFxn)
              fprintf('Running %s...\n',ana.otherFxn{i});
              cfg = ana.cfg_other{i};
              % set the file to be read (data_suffix should be [] on 1st
              % time through)
              cfg.inputfile = fullfile(saveDirRawFile,sprintf('data_raw%s_%s.mat',data_suffix,eventVal));
              % run the function
              ft_raw.(eventVal) = feval(str2func(ana.otherFxn{i}),cfg);
              
              % add the new data_suffix for the outputfile
              data_suffix = cat(2,data_suffix,sprintf('_%s',cfg.ftype));
              
              cfg = [];
              % set the file to save
              cfg.outputfile = fullfile(saveDirRawFile,sprintf('data_raw%s_%s.mat',data_suffix,eventVal));
              % preprocess it again
              ft_raw.(eventVal) = ft_preprocessing(cfg,ft_raw.(eventVal));
              fprintf('Done with %s.\n',ana.otherFxn{i});
            end
          end
        else
          fprintf('No good trials found for %s, %s, %s!\n',exper.subjects{sub},sesStr,eventVal);
        end % if there were trials
      end % for cVal
    end % if ~isempty(exper.eventValuesExtra.newValue)
    
    % if we found an eventValue with 0 trials, put original events back so
    % the next subject/session will be processed correctly
    if exist('emptyEvents','var')
      %exper.eventValues = orig.eventValues;
      if ~isempty(exper.eventValuesExtra.newValue)
        exper.eventValuesExtra = orig.eventValuesExtra;
      end
      clear emptyEvents
    end
    
  end % for exper.sessions
end % for exper.subjects

%% if necessary, overwrite exper.eventValues so it includes the extra values
if ~isempty(exper.eventValuesExtra.newValue) && ~exper.eventValuesExtra.onlyKeepExtras
  exper.eventValues = sort(eventValuesWithExtra);
elseif ~isempty(exper.eventValuesExtra.newValue) && exper.eventValuesExtra.onlyKeepExtras
  exper.eventValues = sort(cat(2,exper.eventValuesExtra.newValue{:}));
end

end
