function [exper] = create_ft_struct_multiSes(ana,cfg_pp,exper,dirs,files)
%CREATE_FT_STRUCT: Create a FieldTrip data structure from NS files
%
% [exper] = create_ft_struct_multiSes(ana,cfg_pp,exper,dirs,files)
%
% The default function (set in ana.segFxn) that this function calls is
% SEG2FT. See 'HELP SEG2FT' for more information.
%
% For multiple sessions, SEG2FT will collapse across them if they're input
% as {{'ses1','ses2'}}. If they're input as {{'ses1'},{'ses2'}} then they
% will be kept separate.
%
% ana.artifact.type = the type of artifact info to use: 'none', 'nsAuto',
% 'zeroVar', 'preRejManual', ftManual', 'ftICA', 'badChanManual',
% 'badChanEP', and/or 'rmBadChan'.
%
%                     Some types can be combined.
%
%                     If 'nsAuto', SEG2FT expects to find a NS bci metadata
%                     file. This file denotes the reject artifact trials.
%                     Use the File Export tool to export Metadata > Segment
%                     Information; put it in a 'ns_bci' directory in
%                     dataroot/session.
%
%                     See SEG2FT for more artifact information.
%                     (default = {'none'})
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

% make sure exper.eventValues is a cell
if ~iscell(exper.eventValues)
  exper.eventValues = {{exper.eventValues}};
elseif iscell(exper.eventValues) && ~iscell(exper.eventValues{1})
  exper.eventValues = {exper.eventValues};
elseif iscell(exper.eventValues) && iscell(exper.eventValues{1})
  fprintf('good eventValues setup!\n');
end

% make sure exper.sessions is a cell
if ~iscell(exper.sessions)
  exper.sessions = {{exper.sessions}};
elseif iscell(exper.sessions) && ~iscell(exper.sessions{1})
  exper.sessions = {exper.sessions};
elseif iscell(exper.sessions) && iscell(exper.sessions{1})
  fprintf('good sessions setup!\n');
end

% need a segmentation function
if ~isfield(ana,'segFxn')
  ana.segFxn = 'seg2ft';
end

if ~isfield(ana,'trialFxn');
  ana.trialFxn = 'seg_trialfun';
end

if ~isfield(ana,'continuous')
  ana.continuous = 'no';
end

% maximum size of FieldTrip cfg struct
if ~isfield(ana,'checksize')
  ana.checksize = 1e5;
end

% for adding metadata; false by default
if ~isfield(ana,'useMetadata')
  ana.useMetadata = false;
  ana.metadata.types = {};
elseif isfield(ana,'useMetadata') && ana.useMetadata
  if ~isfield(ana,'metadata')
    error('No ana.metadata.types field');
  elseif isfield(ana,'metadata') && ~isfield(ana.metadata,'types')
    error('No ana.metadata.types field');
  elseif isfield(ana,'metadata') && isfield(ana.metadata,'types') && isempty(ana.metadata.types)
    error('ana.metadata.types field is empty');
  end
elseif isfield(ana,'useMetadata') && ~ana.useMetadata
  ana.metadata.types = {};
end

% use info like ana.trl_order, ana.trl_order, ana.sessionNames,
% ana.phaseNames, exper.eventValues, and exper.prepost (usually when
% segmenting continuous data)
if ~isfield(ana,'useExpInfo')
  ana.useExpInfo = false;
end

if ~isfield(ana,'usePhotodiodeDIN')
  ana.usePhotodiodeDIN = false;
end
if ana.usePhotodiodeDIN && ~isfield(ana,'photodiodeDIN_thresholdMS')
  ana.photodiodeDIN_thresholdMS = 75;
end
if ana.usePhotodiodeDIN && ~isfield(ana,'photodiodeDIN_str')
  ana.photodiodeDIN_str = 'DIN ';
end

% need an artifact detection type ('none', or: 'nsAuto', 'preRejManual', 'ftManual', 'ftICA')
if ~isfield(ana,'artifact') || (isfield(ana,'artifact') && ~isfield(ana.artifact,'type'))
  ana.artifact.type = {'none'};
elseif isfield(ana,'artifact') && isfield(ana.artifact,'type') && ischar(ana.artifact.type)
  ana.artifact.type = {ana.artifact.type};
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

%% make the field names legal

% some of the many characters that cannot be used in a struct field, to be
% replaced with an underscore
illegalStructFieldChars = {' ','-','+','=','!','@','#','$','%','^','&','*','?','~','(',')','{','}','[',']'};
replaceIllegalCharWith = '_';

% appent this character in front of a field that starts with a number
appendInFrontOfNum = 'a';

%% initialize for storing some information

% store the bad channel information
exper.badChan = cell(length(exper.subjects),length(exper.sessions));
% store the bad event information
exper.badEv = cell(length(exper.subjects),length(exper.sessions));

%% Store the raw data in FieldTrip format

fprintf('Converting NetStation to FieldTrip...\n');

% initialize for the sesStr
exper.sesStr = cell(size(exper.sessions));

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    
    % make sure event values are sorted
    if ~issorted(exper.eventValues{ses})
      exper.eventValues{ses} = sort(exper.eventValues{ses});
    end
    
    % back up original field values because that is how we access the data
    eventValues_orig = exper.eventValues{ses};
    for i = 1:length(exper.eventValues{ses})
      % look for illegal characters
      for ic = 1:length(illegalStructFieldChars)
        if ~isempty(strfind(exper.eventValues{ses}{i},illegalStructFieldChars{ic}))
          exper.eventValues{ses}{i} = strrep(exper.eventValues{ses}{i},illegalStructFieldChars{ic},replaceIllegalCharWith);
        end
      end
      % cannot start a struct field with a number
      if isstrprop(exper.eventValues{ses}{i}(1),'digit')
        exper.eventValues{ses}{i} = [appendInFrontOfNum,exper.eventValues{ses}{i}];
      end
    end
    
    % turn the session name into a string
    if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
      sesStr = exper.sessions{ses}{1};
      for i = 2:length(exper.sessions{ses})
        sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
      end
    elseif ~iscell(exper.sessions{ses})
      sesStr = exper.sessions{ses};
    elseif iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1
      sesStr = cell2mat(exper.sessions{ses});
    end
    % store the sesStr
    exper.sesStr{ses} = sesStr;
    
    % initialize to store trial counts
    for evVal = 1:length(exper.eventValues{ses})
      exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = zeros(length(exper.subjects),length(exper.sessions{ses}));
    end
    
    % set the location to save the data and make sure it exists
    saveDirRawFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
    if ~exist(saveDirRawFile,'dir')
      mkdir(saveDirRawFile);
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
        
        eventValuesToProcess = exper.eventValues{ses}(~ismember(exper.eventValues{ses},rawEvents));
        eventValuesToProcess_orig = eventValues_orig(~ismember(exper.eventValues{ses},rawEvents));
        
      else
        eventValuesToProcess = exper.eventValues{ses};
        eventValuesToProcess_orig = eventValues_orig;
      end
    else
      eventValuesToProcess = exper.eventValues{ses};
      eventValuesToProcess_orig = eventValues_orig;
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
      [ft_raw,badChan,badEv] = feval(str2func(ana.segFxn),fullfile(dirs.dataroot,dirs.dataDir),exper.subjects{sub},exper.sessions{ses},eventValuesToProcess,eventValuesToProcess_orig,files.elecfile,ana,exper,dirs);
      
      if ~ana.overwrite.raw
        % load in the ones we didn't process
        if ~isempty(rawFiles)
          for rf = 1:length(rawFiles)
            ft_raw.(rawEvents{rf}) = load(fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',rawEvents{rf})));
          end
        end
      end
    end % if
    
    % store the bad channel information
    exper.badChan{sub,ses} = badChan;
    % store the bad event information
    exper.badEv{sub,ses} = badEv;
    
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
    
    % if we want to keep the events in exper.eventValues (regardless of
    % whether there were any extras)
    for evVal = 1:length(exper.eventValues{ses})
      % get the name of this event type
      eventVal = exper.eventValues{ses}{evVal};
      
      if size(ft_raw.(eventVal).trial,2) > 0
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
