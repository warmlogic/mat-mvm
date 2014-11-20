function [exper] = create_ft_struct_multiSes(ana,cfg_pp,exper,dirs,files)
%CREATE_FT_STRUCT_MULTISES: Create a FieldTrip data structure from NS files
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
  %error('Prevention of overwriting is not fully implemented. You must set overwrite to 1.');
  warning('Prevention of overwriting raw files is not fully implemented. Be careful.');
end

% make sure exper.eventValues is a cell containing one cell per session
if ~iscell(exper.eventValues)
  exper.eventValues = {{exper.eventValues}};
elseif iscell(exper.eventValues) && ~iscell(exper.eventValues{1})
  exper.eventValues = {exper.eventValues};
  % elseif iscell(exper.eventValues) && iscell(exper.eventValues{1})
  %   fprintf('good eventValues setup!\n');
end

% make sure exper.sessions is a cell
if ~iscell(exper.sessions)
  exper.sessions = {{exper.sessions}};
elseif iscell(exper.sessions) && ~iscell(exper.sessions{1})
  exper.sessions = {exper.sessions};
  % elseif iscell(exper.sessions) && iscell(exper.sessions{1})
  %   fprintf('good sessions setup!\n');
end

if length(exper.eventValues) ~= length(exper.sessions)
  error('Length of exper.eventValues must equal the number of sessions in exper.sessions!');
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

% checks if sample numbers in sampleinfo field overlap (after segmentation
% of continuous data)
if ~isfield(ana,'allowTrialOverlap')
  ana.allowTrialOverlap = false;
end

% if true, will overwrite sampleinfo field  (after segmentation of
% continuous data) so trial sample indices are contiguously numbered
if ~isfield(ana,'renumberSamplesContiguous')
  ana.renumberSamplesContiguous = false;
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

% saving FT events for a MFF file can take a long time, so this is a way to
% have a multi-node cluster exit early if you want other EEG processing.
if strcmpi(exper.eegFileExt,'mff')
  if ~isfield(ana,'returnAfterSavingFtEvents')
    ana.returnAfterSavingFtEvents = false;
  end
  if ana.returnAfterSavingFtEvents
    warning('ana.returnAfterSavingFtEvents=true; Changing settings so EEG is not processed!');
    % set these so we don't do any unnecessary processing
    ana.cfg_cont = [];
    ana.artifact.continuousRepair = false;
    ana.artifact.continuousReject = false;
    ana.artifact.continuousICA = false;
  end
else
  ana.returnAfterSavingFtEvents = false;
end

% use info like ana.trl_order, ana.trl_order, ana.sessionNames,
% ana.phaseNames, exper.eventValues, and exper.prepost (usually when
% segmenting continuous data)
if ~isfield(ana,'useExpInfo')
  ana.useExpInfo = false;
end

% millisecond tolerance for how close together two events can appear in the
% evt file, but fieldtrip will collapse them together into a single event
% when reading the flags from the raw file; default = 2 samples at this
% sampling rate
if ana.useExpInfo
  if ~isfield(ana,'evtToleranceMS')
    ana.evtToleranceMS = ceil((2 * 1000) / exper.sampleRate);
  end
end

if ~isfield(ana,'usePhotodiodeDIN')
  ana.usePhotodiodeDIN = false;
end
% how long after the event of interest can the DIN show up? Seems to occur
% within about 6 to 9 ms, but can be delayed as much as 25 (or more)
if ana.usePhotodiodeDIN && ~isfield(ana,'photodiodeDIN_toleranceMS')
  ana.photodiodeDIN_toleranceMS = 20;
end
if ana.usePhotodiodeDIN && ~isfield(ana,'photodiodeDIN_str')
  ana.photodiodeDIN_str = 'DIN ';
end

if ~isfield(ana,'offsetMS')
  ana.offsetMS = 0;
end

% need an artifact detection type ('none', or: 'nsAuto', 'preRejManual', 'ftManual', 'ftICA')
if ~isfield(ana,'artifact')
  ana.artifact.type = {'none'};
elseif isfield(ana,'artifact')
  if isfield(ana.artifact,'type')
    if isempty(ana.artifact.type)
      ana.artifact.type = {'none'};
    elseif ~isempty(ana.artifact.type) && ischar(ana.artifact.type)
      ana.artifact.type = {ana.artifact.type};
    end
  elseif ~isfield(ana.artifact,'type')
    ana.artifact.type = {'none'};
  end
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

%% Store the raw data in FieldTrip format

fprintf('Converting NetStation to FieldTrip...\n');

% initialize for the sesStr
exper.sesStr = cell(size(exper.sessions));

for ses = 1:length(exper.sessions)
  
  if size(exper.prepost{ses},1) ~= length(exper.eventValues{ses})
    error('Number of entries in exper.prepost{ses} does not match that of exper.eventValues{ses}!\nConstruct exper.prepost as a cell with one Nx2 matrix per session where N is length(exper.eventValues{ses}); e.g. exper.prepost = {[-1.0 2.0; -1.0 2.0], [-1.0 2.0; -1.0 2.0]} for two sessions with two events each.');
  end
  
  % make sure event values are sorted
  if ~issorted(exper.eventValues{ses})
    [exper.eventValues{ses}, evInd] = sort(exper.eventValues{ses});
    exper.prepost{ses} = exper.prepost{ses}(evInd,:);
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
  
  % initialize to store the bad channel information
  exper.badChan.(sesStr) = cell(length(exper.subjects),1);
  
  % initialize to store trial counts and bad events (if using artifacts)
  for evVal = 1:length(exper.eventValues{ses})
    % get the name of this event type
    eventVal = exper.eventValues{ses}{evVal};
    
    exper.nTrials.(sesStr).(eventVal) = zeros(length(exper.subjects),1);
    exper.badEv.(sesStr).(eventVal) = cell(length(exper.subjects),1);
    if ~ismember('none',ana.artifact.type)
      %exper.artifacts.(sesStr).(eventVal) = cell(length(exper.subjects),1);
      exper.artifacts.(sesStr).(eventVal) = [];
    end
    if ana.useExpInfo
      exper.trialinfo_allEv.(sesStr).(eventVal) = cell(length(exper.subjects),1);
    end
  end
  
  for sub = 1:length(exper.subjects)
    % set the location to save the data and make sure it exists
    saveDirRawFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
    if ~exist(saveDirRawFile,'dir')
      mkdir(saveDirRawFile);
    end
    
    % if we don't want to overwrite any files...
    if ~ana.overwrite.raw
      subDetailFile = dir(fullfile(saveDirRawFile,'subjectDetails.mat'));
      if ~isempty(subDetailFile)
        if length(subDetailFile) == 1
          sd = load(fullfile(saveDirRawFile,subDetailFile.name));
        else
          error('More than one subject detail file found. This should not be possible.');
        end
      end
      
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
        if ~isempty(rawEvents)
          for rf = 1:length(rawEvents)
            ft_raw.(rawEvents{rf}) = load(fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',rawEvents{rf})));
            ft_raw.(rawEvents{rf}) = ft_raw.(rawEvents{rf}).data;
          end
        end
        
        % reconstruct other data for this single subject
        badChan = sd.exper.badChan.(sesStr){1};
        for evVal = 1:length(exper.eventValues{ses})
          % get the name of this event type
          eventVal = exper.eventValues{ses}{evVal};
          
          badEv.(eventVal) = sd.exper.badEv.(sesStr).(eventVal){1};
          if isfield(sd.exper,'artifacts')
            artTypes = fieldnames(sd.exper.artifacts.(sesStr).(eventVal));
            for at = 1:length(artTypes)
              artifacts.(eventVal).(artTypes{at}) = sd.exper.artifacts.(sesStr).(eventVal).(artTypes{at}){1};
            end
          end
          if ana.useExpInfo
            trialinfo_allEv.(eventVal) = sd.exper.trialinfo_allEv.(sesStr).(eventVal){1};
          end
        end
      end
      
    else
      if ~ana.overwrite.raw
        error('Support not implemented for adding event values!');
      end
      
      fprintf('Creating FT struct of raw EEG data: %s, %s%s.\n',exper.subjects{sub},sesStr,sprintf(repmat(', ''%s''',1,length(eventValuesToProcess)),eventValuesToProcess{:}));
      
      % collect all the raw data
      [ft_raw,badChan,badEv,artifacts,trialinfo_allEv] = feval(str2func(ana.segFxn),fullfile(dirs.dataroot,dirs.dataDir),exper.subjects{sub},exper.sessions{ses},ses,eventValuesToProcess,eventValuesToProcess_orig,exper.prepost{ses},files.elecfile,ana,exper,dirs);
      if ana.returnAfterSavingFtEvents
        warning('Returning after having saved FieldTrip events. Not running %s to completion!',mfilename);
        continue
      end
      
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
    exper.badChan.(sesStr){sub} = badChan;
    
    %     % check on any empty events
    %     for evVal = 1:length(eventValuesToProcess)
    %       % get the name of this event type
    %       eventVal = eventValuesToProcess{evVal};
    %       if ~isempty(ft_raw.(eventVal).trial)
    %         % store the number of trials for this event value
    %         exper.nTrials.(eventVal)(sub,ses) = size(ft_raw.(eventVal).trial,2);
    %
    %         fprintf('Created FT struct of raw EEG data: %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
    %       else
    %         warning([mfilename,':noTrialsFound'],'Did not create FT struct of raw EEG data: %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
    %         if ~exist('emptyEvents','var')
    %           emptyEvents = {eventVal};
    %         elseif exist('emptyEvents','var')
    %           emptyEvents = cat(2,emptyEvents,eventVal);
    %         end
    %       end
    %     end
    
    % store data about events and do any additional preprocessing
    for evVal = 1:length(exper.eventValues{ses})
      % get the name of this event type
      eventVal = exper.eventValues{ses}{evVal};
      
      if size(ft_raw.(eventVal).trial,2) > 0
        % store the number of trials for this event value
        exper.nTrials.(sesStr).(eventVal)(sub) = size(ft_raw.(eventVal).trial,2);
        exper.badEv.(sesStr).(eventVal){sub} = badEv.(eventVal);
        if ~isempty(artifacts)
          artTypes = fieldnames(artifacts.(eventVal));
          for at = 1:length(artTypes)
            exper.artifacts.(sesStr).(eventVal).(artTypes{at}){sub} = artifacts.(eventVal).(artTypes{at});
          end
        end
        if ana.useExpInfo && ~isempty(trialinfo_allEv)
          exper.trialinfo_allEv.(sesStr).(eventVal){sub} = trialinfo_allEv.(eventVal);
        end
        
        if ana.overwrite.raw
          % set outputfile so the raw data is saved; especially useful for
          % implementing analyses with the peer toolbox
          cfg_pp.outputfile = fullfile(saveDirRawFile,sprintf('data_raw_%s.mat',eventVal));
          
          % do any preprocessing
          fprintf('Running ft_preprocessing on %s, %s, %s...\n',exper.subjects{sub},sesStr,eventVal);
          ft_preprocessing(cfg_pp,ft_raw.(eventVal));
        else
          cfg_pp = sd.cfg_pp;
          warning('Not re-running ft_preprocessing on %s, %s, %s. Overwriting current cfg_pp with old cfg_pp.',exper.subjects{sub},sesStr,eventVal);
        end
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
    
    % save a subjectDetails.mat file for each session
    saveFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr,'subjectDetails.mat');
    % back up original exper struct
    exper_orig = exper;
    sesEvValues = exper.eventValues{ses};
    % include only this subject
    exper.subjects = exper.subjects(sub);
    exper.sessions = exper.sessions(ses);
    % % do not discard all of this because of how mm_loadAD works
    %exper.eventValues = exper.eventValues(ses);
    %exper.prepost = exper.prepost(ses);
    %exper.sesStr = {sesStr};
    exper.badChan.(sesStr) = exper_orig.badChan.(sesStr)(sub);
    for evVal = 1:length(sesEvValues)
      exper.nTrials.(sesStr).(sesEvValues{evVal}) = exper_orig.nTrials.(sesStr).(sesEvValues{evVal})(sub);
      exper.badEv.(sesStr).(sesEvValues{evVal}) = exper_orig.badEv.(sesStr).(sesEvValues{evVal})(sub);
      if ~isempty(artifacts)
        artTypes = fieldnames(exper_orig.artifacts.(sesStr).(sesEvValues{evVal}));
        for at = 1:length(artTypes)
          exper.artifacts.(sesStr).(sesEvValues{evVal}).(artTypes{at}) = exper_orig.artifacts.(sesStr).(sesEvValues{evVal}).(artTypes{at})(sub);
        end
      end
      if ana.useExpInfo && ~isempty(trialinfo_allEv)
        exper.trialinfo_allEv.(sesStr).(sesEvValues{evVal}) = exper_orig.trialinfo_allEv.(sesStr).(sesEvValues{evVal})(sub);
      end
    end
    %fn = fieldnames(exper.nTrials);
    %for f = 1:length(fn)
    %  if f ~= ses
    %    if isfield(exper.badChan,fn{f})
    %      exper.badChan = rmfield(exper.badChan,fn{f});
    %    end
    %    if isfield(exper.nTrials,fn{f})
    %      exper.nTrials = rmfield(exper.nTrials,fn{f});
    %    end
    %    if isfield(exper.badEv,fn{f})
    %      exper.badEv = rmfield(exper.badEv,fn{f});
    %    end
    %  end
    %end
    save(saveFile,'exper','ana','dirs','files','cfg_pp');
    % restore original exper struct
    exper = exper_orig;
    
  end % for exper.subjects
end % for exper.sessions

end
