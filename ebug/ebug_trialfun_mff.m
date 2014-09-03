function trl = ebug_trialfun_mff(cfg)

% operates using Net Station evt files and event structs

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
fprintf('Reading flags from EEG file using FieldTrip...');
ft_hdr = ft_read_header(cfg.dataset);
[pathstr,name] = fileparts(cfg.dataset);
ftEventsFile = fullfile(pathstr,sprintf('%s_ftEvents.mat',name));
if exist(ftEventsFile,'file')
  ft_event = load(ftEventsFile);
  if isfield(ft_event,'date_string')
    warning('Using pre-saved FT events from this date: %s!',ft_event.date_string);
  else
    warning('Using pre-saved FT events from an unknown date!');
  end
  ft_event = ft_event.ft_event;
else
  tic
  ft_event = ft_read_event(cfg.dataset);
  toc
  date_string = datestr(now);
  fprintf('Saving FT events from MFF (current time: %s): %s...',date_string,ftEventsFile);
  save(ftEventsFile,'ft_event','date_string');
  fprintf('Done.\n');
end
fprintf('Done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in external data, if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg.eventinfo.useMetadata
  md = cfg.eventinfo.metadata;
  
  if ismember('eventStruct', md.types)
    % read the events file
    eventsFile = fullfile(md.dirs.dataroot,md.dirs.behDir,md.subject,'events','events.mat');
    if exist(eventsFile,'file')
      fprintf('Loading events file: %s...',eventsFile);
      events_all = load(eventsFile,'events');
      events_all = events_all.events;
      fprintf('Done.\n');
    else
      error('Cannot find events file: %s\n',eventsFile)
    end
  end
  
  if ismember('expParam', md.types)
    % read the experiment parameters file
    expParamFile = fullfile(md.dirs.dataroot,md.dirs.behDir,md.subject,'experimentParams.mat');
    if exist(expParamFile,'file')
      fprintf('Loading experiment parameters file: %s...',expParamFile);
      load(expParamFile,'expParam');
      fprintf('Done.\n');
    else
      error('Cannot find experiment parameters file: %s\n',expParamFile)
    end
  end
end

if cfg.eventinfo.usePhotodiodeDIN
  triggers = {'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND', cfg.eventinfo.photodiodeDIN_str};
else
  triggers = {'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND'};
end

% initialize the trl matrix
trl = [];

% all trls need to have the same length
maxTrlCols = -Inf;
fn_trl_ord = fieldnames(cfg.eventinfo.trl_order);
for fn = 1:length(fn_trl_ord)
  if ismember(fn_trl_ord{fn},cfg.eventinfo.eventValues)
    if length(cfg.eventinfo.trl_order.(fn_trl_ord{fn})) > maxTrlCols
      maxTrlCols = length(cfg.eventinfo.trl_order.(fn_trl_ord{fn}));
    end
  end
end
if maxTrlCols == -Inf
  error('Did not set maximum number of trialinfo columns!\n');
  %keyboard
end
timeCols = 3;
trl_ini = -1 * ones(1, timeCols + maxTrlCols);

if cfg.eventinfo.usePhotodiodeDIN
  photodiodeDIN_toleranceMS = cfg.eventinfo.photodiodeDIN_toleranceMS;
  photodiodeDIN_toleranceSamp = ceil((photodiodeDIN_toleranceMS / 1000) * ft_hdr.Fs);
end

offsetSamp = ceil((cfg.eventinfo.offsetMS / 1000) * ft_hdr.Fs);

% ft_event_ind = true(1,length(ft_event));
% for i = 1:length(ft_event)
%   if isempty(ft_event(i).value)
%     ft_event_ind(i) = false;
%   end
% end
% ft_event = ft_event(ft_event_ind);

ft_event = ft_event(ismember({ft_event.type},cfg.trialdef.eventtype));

% only keep the ft events with triggers
ft_event = ft_event(ismember({ft_event.value},triggers));

%% things to change for MFF

% ns_evt{
% ft_event(i).orig.keys(

% +1}{ec}
% ).key.data.data

% +1}(ec)
% ).key.data.data

% ec
% i

%% go through events and add metadata to trl matrix

ses = cfg.eventinfo.sessionNum;
sesName = cfg.eventinfo.sessionNames{ses};
sesType = find(ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses}));
% sesType = ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses});

fprintf('FT event count of NS flags (out of %d): %s',length(ft_event),repmat(' ',1,length(num2str(length(ft_event)))));

for i = 1:length(ft_event)
  fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
  
  if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
    % found an EEG event that we might want to process
    
    switch ft_event(i).value
      case 'STIM'
        nKeys = length(ft_event(i).orig.keys);
        
        % set column types because Net Station evt files can vary
        ns_evt_cols = {};
        for ns = 1:nKeys
          ns_evt_cols = cat(1,ns_evt_cols,ft_event(i).orig.keys(ns).key.keyCode);
        end
        cols.isExp = find(strcmp(ns_evt_cols,'expt'));
        cols.phase = find(strcmp(ns_evt_cols,'phas'));
        if isempty(cols.isExp) || isempty(cols.phase)
          error('isExp and phase columns not found!');
          %keyboard
        end
        
        % only get data from the real experiment, not the practice
        if strcmpi(ft_event(i).orig.keys(cols.isExp).key.data.data,'true')
          phaseName = ft_event(i).orig.keys(cols.phase).key.data.data;
          phaseType = find(ismember(cfg.eventinfo.phaseNames{sesType},phaseName));
          
          % set column types because Net Station evt files can vary
          cols.phaseCount = find(strcmp(ns_evt_cols,'pcou'));
          cols.trial = find(strcmp(ns_evt_cols,'trln'));
          
          if length(phaseType) > 1
            % choose the phase that corresponds to this count, not that it
            % really matters because they all have the same name...
            phaseType = phaseType(str2double(ft_event(i).orig.keys(cols.phaseCount).key.data.data));
          end
          
          switch phaseName
            case {'match', 'prac_match'}
              match_responses = {'same', 'diff'};
              
              cols.snum = find(strcmp(ns_evt_cols,'snum'));
              
              if isempty(cols.phaseCount) || isempty(cols.trial) || isempty(cols.snum)
                error('match: phaseCount, trial, or snum columns not found!');
                %keyboard
              end
              
            case {'nametrain', 'prac_nametrain', 'name', 'prac_name'}
         
              cols.block = find(strcmp(ns_evt_cols,'bloc'));
              
              if isempty(cols.phaseCount) || isempty(cols.trial) || isempty(cols.block)
                error('nametrain: phaseCount, trial, or block column not found!');
                %keyboard
              end
          end
          
          % set the phase name with phase count for event struct
          %phaseName_pc = sprintf('%s_%s',phaseName,ft_event(i).orig.keys(cols.phaseCount).key.data.data);
          phaseName_pc = sprintf('%s_%d',phaseName,str2double(ft_event(i).orig.keys(cols.phaseCount).key.data.data));
          
          % if we want to process this phase (set up in top-level script)
          if ismember(phaseName,cfg.eventinfo.phaseNames{sesType})
            
            % Critical: set up the stimulus type, as well as the event
            % string to match eventValues
            switch phaseName
              case {'match', 'prac_match'}
                stimNum = str2double(ft_event(i).orig.keys(cols.snum).key.data.data);
                if stimNum == 1
                  stimType = 'MATCH_STIM1';
                elseif stimNum == 2
                  stimType = 'MATCH_STIM2';
                else
                  error('stimNum %d does not correspond to MATCH_STIM!',stimNum);
                  %keyboard
                end
                evVal = 'match_stim';
                
              case {'nametrain', 'prac_nametrain', 'name', 'prac_name'}
                stimType = 'NAME_STIM';
                if strcmp(phaseName,'nametrain')
                  evVal = 'nametrain_stim';
                elseif strcmp(phaseName,'name')
                  evVal = 'name_stim';
                end
            end
            
            % get the order of trl columns for this phase and event type
            trl_order = cfg.eventinfo.trl_order.(evVal);
            
            % find where this event type occurs in the list
            eventNumber = find(ismember(cfg.trialdef.eventvalue,evVal));
            if isempty(eventNumber)
              eventNumber = -1;
            end
            
            if length(eventNumber) == 1 && eventNumber ~= -1
              % set the times we need to segment before and after the
              % trigger
              prestimSec = abs(cfg.eventinfo.prepost(eventNumber,1));
              poststimSec = cfg.eventinfo.prepost(eventNumber,2);
              
              % prestimulus period should be negative
              prestimSamp = -round(prestimSec * ft_hdr.Fs);
              poststimSamp = round(poststimSec * ft_hdr.Fs);
            else
              error('event number not found for %s!\n',evVal);
              %fprintf('event number not found for %s!\n',evVal);
              %keyboard
            end
            
            % find the entry in the event struct
            switch phaseName
              case {'match', 'prac_match'}
                this_event = events_all.(sesName).(phaseName_pc).data(...
                  [events_all.(sesName).(phaseName_pc).data.isExp] == 1 &...
                  ismember({events_all.(sesName).(phaseName_pc).data.phaseName},{phaseName}) &...
                  [events_all.(sesName).(phaseName_pc).data.phaseCount] == str2double(ft_event(i).orig.keys(cols.phaseCount).key.data.data) &...
                  ismember({events_all.(sesName).(phaseName_pc).data.type},{stimType}) &...
                  [events_all.(sesName).(phaseName_pc).data.trial] == str2double(ft_event(i).orig.keys(cols.trial).key.data.data) ...
                  );
                
              case {'nametrain', 'prac_nametrain', 'name', 'prac_name'}
                this_event = events_all.(sesName).(phaseName_pc).data(...
                  [events_all.(sesName).(phaseName_pc).data.isExp] == 1 &...
                  ismember({events_all.(sesName).(phaseName_pc).data.phaseName},{phaseName}) &...
                  [events_all.(sesName).(phaseName_pc).data.phaseCount] == str2double(ft_event(i).orig.keys(cols.phaseCount).key.data.data) &...
                  ismember({events_all.(sesName).(phaseName_pc).data.type},{stimType}) &...
                  [events_all.(sesName).(phaseName_pc).data.trial] == str2double(ft_event(i).orig.keys(cols.trial).key.data.data) &...
                  [events_all.(sesName).(phaseName_pc).data.block] == str2double(ft_event(i).orig.keys(cols.block).key.data.data) ...
                  );
            end
            
            if length(this_event) > 1
              error('More than one event found (i=%d)! Fix this script before continuing analysis.',i)
              %warning('More than one event found (i=%d)! Fix this script before continuing analysis.',i)
              %keyboard
            elseif isempty(this_event)
              error('No event found (i=%d)! Fix this script before continuing analysis.',i)
              %warning('No event found (i=%d)! Fix this script before continuing analysis.',i)
              %keyboard
            elseif length(this_event) == 1
              
              switch phaseName
                case {'match', 'prac_match'}
                  phaseCount = this_event.phaseCount;
                  trial = this_event.trial;
                  familyNum = this_event.familyNum;
                  speciesNum = this_event.speciesNum;
                  exemplarNum = this_event.exemplarNum;
                  % stimNum = str2double(this_event.type(end));
                  isSubord = this_event.isSubord;
                  trained = this_event.trained;
                  sameSpecies = this_event.sameSpecies;
                  response = ismember(match_responses,this_event.resp);
                  if any(response)
                    response = find(response);
                  elseif ~any(response) && strcmp(this_event.resp,'none')
                    response = 0;
                  else
                    error('response (%s) not found in match_responses list',this_event.resp);
                    %keyboard
                  end
                  rt = this_event.rt;
                  acc = this_event.acc;
                  
                case {'nametrain', 'prac_nametrain', 'name', 'prac_name'}
                  phaseCount = this_event.phaseCount;
                  block = this_event.block;
                  trial = this_event.trial;
                  familyNum = this_event.familyNum;
                  speciesNum = this_event.speciesNum;
                  exemplarNum = this_event.exemplarNum;
                  isSubord = this_event.isSubord;
                  response = this_event.resp;
                  rt = this_event.rt;
                  acc = this_event.acc;
              end
              
              % add it to the trial definition
              this_trl = trl_ini;
              
              % get the time of this event
              this_sample = ft_event(i).sample;
              
              % if we're using the photodiode DIN and we find one within
              % the threshold, replace the current sample time with that of
              % the DIN
              if cfg.eventinfo.usePhotodiodeDIN
                if strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str)
                  % if there is a DIN before and after the stim, pick
                  % the closer one
                  preDiff = (ft_event(i-1).sample - this_sample);
                  postDiff = (ft_event(i+1).sample - this_sample);
                  
                  if preDiff < 0 && abs(preDiff) <= photodiodeDIN_toleranceSamp
                    preFlag = true;
                  else
                    preFlag = false;
                  end
                  if postDiff <= photodiodeDIN_toleranceSamp
                    postFlag = true;
                  else
                    postFlag = false;
                  end
                  
                  if preFlag && ~postFlag
                    % only the pre-DIN makes sense
                    this_sample = ft_event(i-1).sample;
                  elseif ~preFlag && postFlag
                    % only the post-DIN makes sense
                    this_sample = ft_event(i+1).sample;
                  elseif preFlag && postFlag
                    % choose the smaller one
                    if abs(preDiff) < abs(postDiff)
                      this_sample = ft_event(i-1).sample;
                    elseif abs(preDiff) > abs(postDiff)
                      this_sample = ft_event(i+1).sample;
                    elseif abs(preDiff) == abs(postDiff)
                      error('same exact difference offset! figure out what to do.');
                      %keyboard
                    end
                  end
                elseif strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && ~strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i+1).sample - this_sample) <= photodiodeDIN_toleranceSamp
                  this_sample = ft_event(i+1).sample;
                elseif strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i-1).sample - this_sample) < 0 && abs(ft_event(i-1).sample - this_sample) <= photodiodeDIN_toleranceSamp
                  % apparently the ethernet tags can be delayed
                  % enough that the DIN shows up first
                  this_sample = ft_event(i-1).sample;
                end
              end
              
              % prestimulus sample
              this_trl(1) = this_sample + prestimSamp + offsetSamp;
              % poststimulus sample
              this_trl(2) = this_sample + poststimSamp + offsetSamp;
              % offset in samples
              this_trl(3) = prestimSamp;
              
              for to = 1:length(trl_order)
                thisInd = find(ismember(trl_order,trl_order{to}));
                if ~isempty(thisInd)
                  if exist(trl_order{to},'var')
                    this_trl(timeCols + thisInd) = eval(trl_order{to});
                  else
                    error('variable %s does not exist (i = %d)!\n',trl_order{to},i);
                    %fprintf('variable %s does not exist (i = %d)!\n',trl_order{to},i);
                    %keyboard
                  end
                end
              end
              
              % put all the trials together
              trl = cat(1,trl,double(this_trl));
              
            end % only one event
          else
            % keyboard
            continue
          end % is a phase in this session
        end % isExp
        
        %       case 'RESP'
        %
        %       case 'FIXT'
        %
        %       case 'PROM'
        %
        %       case 'REST'
        %
        %       case 'REND'
        %
        %       case cfg.eventinfo.photodiodeDIN_str
    end % switch
  end
end % i

fprintf('\n');
