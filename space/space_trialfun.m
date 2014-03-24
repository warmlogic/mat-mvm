function trl = space_trialfun(cfg)

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
    
    evtToleranceMS = cfg.eventinfo.evtToleranceMS;
    evtToleranceSamp = ceil((cfg.eventinfo.evtToleranceMS / 1000) * ft_hdr.Fs);
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
  fprintf('Did not set maximum number of trialinfo columns!\n');
  keyboard
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
        warning('Index %d does not have the same value! ns_evt: %s, ft_event: %s',i,ns_evt{1}{i},ft_event(i).value);
        keyboard
      end
    end
  end
else
  warning('ns_evt and ft_event are not the same length! Need to fix alignment code.');
  keyboard
end

%% go through events and add metadata to trl matrix

ses = cfg.eventinfo.sessionNum;
sesName = cfg.eventinfo.sessionNames{ses};
% sesType = find(ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses}));
sesType = ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses});

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
          keyboard
        end
        
        % only get data from the real experiment, not the practice
        if strcmpi(ns_evt{cols.isExp+1}{ec},'true')
          phaseName = ns_evt{cols.phase+1}{ec};
          phaseType = find(ismember(cfg.eventinfo.phaseNames{sesType},phaseName));
          
          % set column types because Net Station evt files can vary
          cols.phaseCount = find(strcmp(ns_evt_cols,'pcou'));
          cols.trial = find(strcmp(ns_evt_cols,'trln'));
          switch phaseName
            case {'multistudy', 'prac_multistudy', 'cued_recall', 'prac_cued_recall'}
              cols.type = find(strcmp(ns_evt_cols,'type'));
              if isempty(cols.phaseCount) || isempty(cols.trial) || isempty(cols.type)
                keyboard
              end
            case {'expo', 'prac_expo', 'distract_math', 'prac_distract_math'}
              if isempty(cols.phaseCount) || isempty(cols.trial)
                keyboard
              end
          end
          
          % if we want to process this phase (set up in top-level script)
          if ismember(phaseName,cfg.eventinfo.phaseNames{sesType})
            
            % Critical: set up the stimulus type, as well as the event
            % string to match eventValues
            switch phaseName
              case {'expo', 'prac_expo'}
                
                stimType = 'EXPO_IMAGE';
                evVal = 'expo_stim';
                
              case {'multistudy', 'prac_multistudy'}
                
                if strcmp(ns_evt{cols.type+1}{ec},'word')
                  stimType = 'STUDY_WORD';
                  evVal = 'multistudy_word';
                elseif strcmp(ns_evt{cols.type+1}{ec},'image')
                  stimType = 'STUDY_IMAGE';
                  evVal = 'multistudy_image';
                end
                
              case {'distract_math', 'prac_distract_math'}
                
                stimType = 'MATH_PROB';
                evVal = 'distract_math_stim';
                
              case {'cued_recall', 'prac_cued_recall'}
                
                if strcmp(ns_evt{cols.type+1}{ec},'recognition')
                  stimType = 'RECOGTEST_STIM';
                  evVal = 'cued_recall_stim';
                else
                  keyboard
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
              fprintf('event number not found for %s!\n',evVal);
              keyboard
            end
            
            % find the entry in the event struct
            this_event = events_all.(sesName).(phaseName).data(...
              [events_all.(sesName).(phaseName).data.isExp] == 1 &...
              ismember({events_all.(sesName).(phaseName).data.phaseName},{phaseName}) &...
              [events_all.(sesName).(phaseName).data.phaseCount] == str2double(ns_evt{cols.phaseCount+1}{ec}) &...
              ismember({events_all.(sesName).(phaseName).data.type},{stimType}) &...
              [events_all.(sesName).(phaseName).data.trial] == str2double(ns_evt{cols.trial+1}{ec}) ...
              );
            
            if length(this_event) > 1
              warning('More than one event found! Fix this script before continuing analysis.')
              keyboard
            elseif isempty(this_event)
              warning('No event found! Fix this script before continuing analysis.')
              keyboard
            elseif length(this_event) == 1
              
              switch phaseName
                case {'expo', 'prac_expo'}
                  
                  phaseCount = this_event.phaseCount;
                  trial = this_event.trial;
                  stimNum = this_event.stimNum;
                  i_catNum = this_event.i_catNum;
                  targ = this_event.targ;
                  spaced = this_event.spaced;
                  lag = this_event.lag;
                  expo_response = this_event.resp;
                  cr_recog_acc = this_event.cr_recog_acc;
                  cr_recall_resp = this_event.cr_recall_resp;
                  cr_recall_spellCorr = this_event.cr_recall_spellCorr;
                  rt = this_event.rt;
                  
                case {'multistudy', 'prac_multistudy'}
                  
                  phaseCount = this_event.phaseCount;
                  trial = this_event.trial;
                  stimNum = this_event.stimNum;
                  catNum = this_event.catNum;
                  targ = this_event.targ;
                  spaced = this_event.spaced;
                  lag = this_event.lag;
                  presNum = this_event.presNum;
                  pairOrd = this_event.pairOrd;
                  pairNum = this_event.pairNum;
                  cr_recog_acc = this_event.cr_recog_acc;
                  cr_recall_resp = this_event.cr_recall_resp;
                  cr_recall_spellCorr = this_event.cr_recall_spellCorr;
                  
                case {'distract_math', 'prac_distract_math'}
                  
                  phaseCount = this_event.phaseCount;
                  trial = this_event.trial;
                  response = this_event.resp;
                  acc = this_event.acc;
                  rt = this_event.rt;
                  
                case {'cued_recall', 'prac_cued_recall'}
                  
                  recog_responses = {'old', 'new'};
                  new_responses = {'sure', 'maybe'};
                  
                  phaseCount = this_event.phaseCount;
                  trial = this_event.trial;
                  stimNum = this_event.stimNum;
                  i_catNum = this_event.i_catNum;
                  targ = this_event.targ;
                  spaced = this_event.spaced;
                  lag = this_event.lag;
                  pairNum = this_event.pairNum;
                  
                  if ~isempty(this_event.recog_resp) && ismember(this_event.recog_resp,recog_responses)
                    recog_resp = find(ismember(recog_responses,this_event.recog_resp));
                  else
                    recog_resp = 0;
                  end
                  recog_acc = this_event.recog_acc;
                  recog_rt = this_event.recog_rt;
                  
                  if ~isempty(this_event.new_resp) && ismember(this_event.new_resp,new_responses)
                    new_resp = find(ismember(new_responses,this_event.new_resp));
                  else
                    new_resp = 0;
                  end
                  new_acc = this_event.new_acc;
                  new_rt = this_event.new_rt;
                  
                  if ~isempty(this_event.recall_resp) && ~ismember(this_event.recall_resp,{'NO_RESPONSE'})
                    recall_resp = 1;
                  else
                    recall_resp = 0;
                  end
                  recall_spellCorr = this_event.recall_spellCorr;
                  recall_rt = this_event.recall_rt;
              end
              
              % add it to the trial definition
              this_trl = trl_ini;
              
              % get the time of this event
              this_sample = ft_event(i).sample;
              
              % if we're using the photodiode DIN and we find one
              % within the threshold, replace the current sample time
              % with that of the DIN
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
                      keyboard
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
                    fprintf('variable %s does not exist!\n',trl_order{to});
                    keyboard
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
