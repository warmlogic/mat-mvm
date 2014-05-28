function trl = ebird_trialfun(cfg)

% operates using Net Station evt files and event structs

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
fprintf('Reading flags from EEG file using FieldTrip...');
ft_hdr = ft_read_header(cfg.dataset);
ft_event = ft_read_event(cfg.dataset);
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
  
  if ismember('nsEvt', md.types)
    % read the file with information about all events
    evtDir = 'ns_evt';
    % find the evt file
    evtfile = dir(fullfile(md.dataroot,md.sesName,evtDir,[md.subject,'*.evt']));
    if isempty(evtfile)
      error('Cannot find %s*.evt file in %s',md.subject,fullfile(md.dataroot,md.sesName,evtDir));
    elseif length(evtfile) > 1
      error('More than one %s*.evt file found in %s',md.subject,fullfile(md.dataroot,md.sesName,evtDir));
    elseif length(evtfile) == 1
      evtfile = fullfile(md.dataroot,md.sesName,evtDir,evtfile.name);
    end
    fprintf('Reading evt file: %s...',evtfile);
    % figure out how many columns there are
    fid = fopen(evtfile,'r');
    maxNumCols = -Inf;
    while 1
      % get each line of the file
      tline = fgetl(fid);
      if ischar(tline)
        if strcmp(tline(end),sprintf('\t'))
          tline = tline(1:end-1);
        end
        % since it's tab delimited, split it and count the length
        numCols = length(regexp(tline,'\t','split'));
        % change the max number
        if numCols > maxNumCols
          maxNumCols = numCols;
        end
      else
        break
      end
    end
    fclose(fid);
    % read the evt file
    fid = fopen(evtfile,'r');
    ns_evt = textscan(fid,repmat('%s',1,maxNumCols),'Headerlines',3,'Delimiter','\t');
    fclose(fid);
    fprintf('Done.\n');
    
    %evtToleranceMS = cfg.eventinfo.evtToleranceMS;
    %evtToleranceSamp = ceil((cfg.eventinfo.evtToleranceMS / 1000) * ft_hdr.Fs);
  end
  
  if ismember('expParam', md.types)
    error('loading experimentParams.mat is not yet supported');
  end
end

% read the types of triggers (flags) in the Net Station evt file
if ismember('nsEvt', md.types)
  triggers = unique(ns_evt{1});
else
  %warning('Need to set up %s to find trigger types!',mfilename);
  %keyboard
  if cfg.eventinfo.usePhotodiodeDIN
    triggers = {'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND', cfg.eventinfo.photodiodeDIN_str};
  else
    triggers = {'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND'};
  end
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

% only keep the ft events with triggers
ft_event = ft_event(ismember({ft_event.value},triggers));

%% Alignment 1/3: read evt and put events with the same sample in alphabetical order

% The header, as read by FieldTrip, is (usually) sorted alphabetically when
% events with different values occurred at the same sample (but not
% necessarily the same millisecond). Therefore, events in the Net Station
% evt file need to be rearranged, especially because the next section of
% this function will remove "duplicate" events (defined as events with the
% same value that occurred at the same sample). There is a third section of
% this function to catch any events that should not be in alphabetical
% order.

% backup so we can access original data
ns_evt_orig = ns_evt;

% % wait to change the values
% ecInd = 1:length(ns_evt{1});

for ec = 1:length(ns_evt{1})
  if ec > 1 && ~strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
    this_time_ms = (str2double(ns_evt{5}{ec}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec}(8:9)) * 1000) + (str2double(ns_evt{5}{ec}(11:13)));
    this_time_samp = fix((this_time_ms / 1000) * ft_hdr.Fs);
    prev_time_ms = (str2double(ns_evt{5}{ec-1}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(8:9)) * 1000) + (str2double(ns_evt{5}{ec-1}(11:13)));
    prev_time_samp = fix((prev_time_ms / 1000) * ft_hdr.Fs);
    
    if this_time_samp == prev_time_samp || abs(this_time_ms - prev_time_ms) < (1000 / ft_hdr.Fs)
      if ~issorted({ns_evt{1}{ec-1},ns_evt{1}{ec}})
        % change the values on the fly
        for i = 1:length(ns_evt_orig)
          ns_evt{i}(ec - 1) = ns_evt_orig{i}(ec);
          ns_evt{i}(ec) = ns_evt_orig{i}(ec - 1);
        end
        
        % % wait to change the values
        % prevInd = ec - 1;
        % ecInd(ec) = prevInd;
        % ecInd(prevInd) = ec;
      end
    end
  end
  
end

% % wait to change the values
% for i = 1:length(ns_evt)
%   ns_evt{i} = ns_evt{i}(ecInd);
% end

%% Alignment 2/3: remove duplicate events (from FieldTrip's sample-based perspective)

% If two (or more) events with the same value string occurred at the same
% sample, FieldTrip will only keep one of them. Therefore, we need to
% remove the duplicate events from the Net Station evt file.

ecInd = true(length(ns_evt{1}),1);

% keep track of how many real evt events we have counted
ec = 0;
  
for i = 1:length(ft_event)
  if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
    % found a trigger in the EEG file events; increment index if
    % value is correct.
    
    if ismember(ft_event(i).value,triggers)
      ec = ec + 1;
    else
      continue
    end
    
    if ec > 1 && strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
      this_time_ms = (str2double(ns_evt{5}{ec}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec}(8:9)) * 1000) + (str2double(ns_evt{5}{ec}(11:13)));
      this_time_samp = fix((this_time_ms / 1000) * ft_hdr.Fs);
      prev_time_ms = (str2double(ns_evt{5}{ec-1}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(8:9)) * 1000) + (str2double(ns_evt{5}{ec-1}(11:13)));
      prev_time_samp = fix((prev_time_ms / 1000) * ft_hdr.Fs);
      
      if (this_time_samp == prev_time_samp || abs(this_time_ms - prev_time_ms) < (1000 / ft_hdr.Fs))
        % need to check whether these are the same events; can't rely on
        % comparing current NS and FT events with ~strcmp because they
        % might have the same value even though they're not the same exact
        % event. So, check the event value for the next NS event vs the
        % current FT event, and make it more robust by checking the sample
        % offset
        next_time_ms = (str2double(ns_evt{5}{ec+1}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec+1}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec+1}(8:9)) * 1000) + (str2double(ns_evt{5}{ec+1}(11:13)));
        next_time_samp = fix((next_time_ms / 1000) * ft_hdr.Fs);
        
        if strcmp(ns_evt{1}(ec+1),ft_event(i).value) && abs(next_time_samp - prev_time_samp) == abs(ft_event(i).sample - ft_event(i-1).sample)
          ecInd(ec) = false;
          ec = ec + 1;
        end
      end
      
    end
  end
end

% only keep the ns_evt indices that align with ft_event
for i = 1:length(ns_evt)
  ns_evt{i} = ns_evt{i}(ecInd);
end

%% Alignment 3/3: final comparison of Net Station and FieldTrip data

% make sure the Net Station evt Code column and ft_event value field are in
% the same order. If they are not (e.g., if the events should not have been
% rearranged alphabetically), then put the Net Station events in the same
% order as the FieldTrip events, but only if the neighboring events are
% aligned.

ns_evt_backup = ns_evt;

% verify that our lists have the same event values

% if length(ns_evt{1}) == length(ft_event(~ismember({ft_event.value},'epoc')))
if length(ns_evt{1}) == length(ft_event)
  for i = 1:length(ft_event)
    if ~strcmp(ns_evt{1}{i},ft_event(i).value)
      if strcmp(ns_evt{1}{i-1},ft_event(i-1).value) && strcmp(ns_evt{1}{i},ft_event(i+1).value) && strcmp(ns_evt{1}{i+1},ft_event(i).value) && strcmp(ns_evt{1}{i+2},ft_event(i+2).value)
        % they're just out of order for some reason (probably due to
        % alphabetical sorting), so put them back
        for j = 1:length(ns_evt_orig)
          ns_evt{j}(i) = ns_evt_backup{j}(i + 1);
          ns_evt{j}(i + 1) = ns_evt_backup{j}(i);
        end
        
      else
        error('Index %d does not have the same value! ns_evt: %s, ft_event: %s',i,ns_evt{1}{i},ft_event(i).value);
        %warning('Index %d does not have the same value! ns_evt: %s, ft_event: %s',i,ns_evt{1}{i},ft_event(i).value);
        %keyboard
      end
    end
  end
else
  error('ns_evt and ft_event are not the same length! Need to fix alignment code.');
  %warning('ns_evt and ft_event are not the same length! Need to fix alignment code.');
  %keyboard
end

%% go through events and add metadata to trl matrix

ses = cfg.eventinfo.sessionNum;
sesName = cfg.eventinfo.sessionNames{ses};
sesType = find(ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses}));
% sesType = ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses});

% keep track of how many real evt events we have counted
ec = 0;

fprintf('FT event count of NS flags (out of %d): %s',length(ft_event),repmat(' ',1,length(num2str(length(ft_event)))));

for i = 1:length(ft_event)
  fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
  
  if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
    % found a trigger in the EEG file events; increment index if value is
    % correct.
    
    if ismember(ft_event(i).value,triggers)
      ec = ec + 1;
    else
      continue
    end
    
    switch ft_event(i).value
      case 'STIM'
        % set column types because Net Station evt files can vary
        ns_evt_cols = {};
        for ns = 1:size(ns_evt,2)
          ns_evt_cols = cat(1,ns_evt_cols,ns_evt{ns}(ec));
        end
        cols.isExp = find(strcmp(ns_evt_cols,'expt'));
        cols.phase = find(strcmp(ns_evt_cols,'phas'));
        if isempty(cols.isExp) || isempty(cols.phase)
          error('isExp and phase columns not found!');
          %keyboard
        end
        
        % only get data from the real experiment, not the practice
        if strcmpi(ns_evt{cols.isExp+1}{ec},'true')
          phaseName = ns_evt{cols.phase+1}{ec};
          phaseType = find(ismember(cfg.eventinfo.phaseNames{sesType},phaseName));
          
          % set column types because Net Station evt files can vary
          cols.phaseCount = find(strcmp(ns_evt_cols,'pcou'));
          cols.trial = find(strcmp(ns_evt_cols,'trln'));
          
          if length(phaseType) > 1
            % choose the phase that corresponds to this count, not that it
            % really matters because they all have the same name...
            phaseType = phaseType(str2double(ns_evt{cols.phaseCount+1}{ec}));
          end
          
          switch phaseName
            case {'match', 'prac_match'}
              image_conditions = sort({'color', 'g', 'g_hi8', 'g_lo8', 'normal'});
              match_responses = {'same', 'diff'};
              
              cols.snum = find(strcmp(ns_evt_cols,'snum'));
              
              if isempty(cols.phaseCount) || isempty(cols.trial) || isempty(cols.snum)
                error('match: phaseCount, trial, or snum columns not found!');
                %keyboard
              end
              
            case {'nametrain', 'prac_nametrain', 'name', 'prac_name'}
              image_conditions = {'normal'};
              
              cols.block = find(strcmp(ns_evt_cols,'bloc'));
              
              if isempty(cols.phaseCount) || isempty(cols.trial) || isempty(cols.block)
                error('nametrain: phaseCount, trial, or block column not found!');
                %keyboard
              end
          end
          
          % set the phase name with phase count for event struct
          %phaseName_pc = sprintf('%s_%s',phaseName,ns_evt{cols.phaseCount+1}{ec});
          phaseName_pc = sprintf('%s_%d',phaseName,str2double(ns_evt{cols.phaseCount+1}{ec}));
          
          % if we want to process this phase (set up in top-level script)
          if ismember(phaseName,cfg.eventinfo.phaseNames{sesType})
            
            % Critical: set up the stimulus type, as well as the event
            % string to match eventValues
            switch phaseName
              case {'match', 'prac_match'}
                stimNum = str2double(ns_evt{cols.snum+1}{ec});
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
                  [events_all.(sesName).(phaseName_pc).data.phaseCount] == str2double(ns_evt{cols.phaseCount+1}(ec)) &...
                  ismember({events_all.(sesName).(phaseName_pc).data.type},{stimType}) &...
                  [events_all.(sesName).(phaseName_pc).data.trial] == str2double(ns_evt{cols.trial+1}(ec)) ...
                  );
                
              case {'nametrain', 'prac_nametrain', 'name', 'prac_name'}
                this_event = events_all.(sesName).(phaseName_pc).data(...
                  [events_all.(sesName).(phaseName_pc).data.isExp] == 1 &...
                  ismember({events_all.(sesName).(phaseName_pc).data.phaseName},{phaseName}) &...
                  [events_all.(sesName).(phaseName_pc).data.phaseCount] == str2double(ns_evt{cols.phaseCount+1}(ec)) &...
                  ismember({events_all.(sesName).(phaseName_pc).data.type},{stimType}) &...
                  [events_all.(sesName).(phaseName_pc).data.trial] == str2double(ns_evt{cols.trial+1}(ec)) &...
                  [events_all.(sesName).(phaseName_pc).data.block] == str2double(ns_evt{cols.block+1}(ec)) ...
                  );
            end
            
            if length(this_event) > 1
              error('More than one event found (ec=%d)! Fix this script before continuing analysis.',ec)
              %warning('More than one event found (ec=%d)! Fix this script before continuing analysis.',ec)
              %keyboard
            elseif isempty(this_event)
              error('No event found (ec=%d)! Fix this script before continuing analysis.',ec)
              %warning('No event found (ec=%d)! Fix this script before continuing analysis.',ec)
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
                  imgCond = find(ismember(image_conditions,this_event.imgCond));
                  isSubord = this_event.isSubord;
                  trained = this_event.trained;
                  sameTrained = this_event.sameTrained;
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
                  imgCond = find(ismember(image_conditions,this_event.imgCond));
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
              this_trl(1) = this_sample + prestimSamp;
              % poststimulus sample
              this_trl(2) = this_sample + poststimSamp;
              % offset in samples
              this_trl(3) = prestimSamp;
              
              for to = 1:length(trl_order)
                thisInd = find(ismember(trl_order,trl_order{to}));
                if ~isempty(thisInd)
                  if exist(trl_order{to},'var')
                    this_trl(timeCols + thisInd) = eval(trl_order{to});
                  else
                    error('variable %s does not exist (ec = %d)!\n',trl_order{to},ec);
                    %fprintf('variable %s does not exist (ec = %d)!\n',trl_order{to},ec);
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
