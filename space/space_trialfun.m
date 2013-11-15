function trl = space_trialfun(cfg)

% operates using Net Station evt files

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
ft_hdr = ft_read_header(cfg.dataset);
ft_event = ft_read_event(cfg.dataset);

% offset should be negative
offsetSamp = round(-cfg.trialdef.prestim*ft_hdr.Fs);
% duration should be 1 sample less than the whole length of an event
durationSamp = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*ft_hdr.Fs) - 1;
% TODO: should this be ceil instead of round?

ns_evt = cfg.eventinfo.ns_evt;
events_all = cfg.eventinfo.events;
% expParam = cfg.eventinfo.expParam;

% initialize the trl matrix
trl = [];

%% set up the exposure phase

cols.expo = [];
cols.expo.sub = 7;
cols.expo.ses = 9;
cols.expo.phase = 11;
cols.expo.phaseCount = 13;
cols.expo.isExp = 15;
cols.expo.trial = 17;
cols.expo.type = 19;
% cols.expo.inum = 23;
cols.expo.icat_str = 25;

% cols.expo.icat_num = 27;
% cols.expo.targ = 29;
% cols.expo.spaced = 31;
% cols.expo.lag = 33;
% cols.expo.resp_str = 35;
% cols.expo.rt = 39;
% cols.expo.keypress = 41;

%% set up the multistudy phase

cols.multistudy = [];
cols.multistudy.sub = 7;
cols.multistudy.ses = 9;
cols.multistudy.phase = 11;
cols.multistudy.phaseCount = 13;
cols.multistudy.isExp = 15;
cols.multistudy.trial = 17;

cols.multistudy.type_image = 19;
cols.multistudy.type_word = 27;

cols.multistudy.lag = 32;

%% go through the events

for ses = 1:length(cfg.eventinfo.sessionNames)
  sesName = cfg.eventinfo.sessionNames{ses};
  sesType = find(ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses}));
  
  for pha = 1:length(cfg.eventinfo.phaseNames)
    phaseName = cfg.eventinfo.phaseNames{sesType}{pha};
    phaseType = find(ismember(cfg.eventinfo.phaseNames{sesType},cfg.eventinfo.phaseNames{sesType}{pha}));
    
    switch phaseName
      
      case {'expo'}
        %% process the exposure phase
        
        %expo_events = events_all.(sessionNames{ses}).(phaseName).data;
        
        %trl_order_expo = cfg.eventinfo.trl_order.(phaseName);
        
        % keep track of how many real evt events we have counted
        ec = 0;
        
        for i = 1:length(ft_event)
          if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
            % found a trigger in the evt file; increment index if value is correct.
            
            %if ~ismember(event(i).value,{'epoc'})
            if ismember(ft_event(i).value,{'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND'})
              ec = ec + 1;
            end
            
            switch ft_event(i).value
              case 'STIM'
                if strcmp(ns_evt{cols.(phaseName).isExp}(ec),'expt') && strcmpi(ns_evt{cols.(phaseName).isExp+1}(ec),'true') &&...
                    strcmp(ns_evt{cols.(phaseName).phase}(ec),'phas') && strcmp(ns_evt{cols.(phaseName).phase+1}(ec),phaseName) &&...
                    strcmp(ns_evt{cols.(phaseName).phaseCount}(ec),'pcou') &&...
                    strcmp(ns_evt{cols.(phaseName).type}(ec),'type') && strcmp(ns_evt{cols.(phaseName).type+1}(ec),'image') &&...
                    strcmp(ns_evt{cols.(phaseName).icat_str}(ec),'icts')
                  
                  % type is not actually necessary; I don't think icat_str
                  % is either
                  
                  % find the entry in the event struct
                  this_event = events_all.(sesName).(phaseName).data(...
                    [events_all.(sesName).(phaseName).data.isExp] == 1 &...
                    ismember({events_all.(sesName).(phaseName).data.type},{'EXPO_IMAGE'}) &...
                    ismember({events_all.(sesName).(phaseName).data.phaseName},{phaseName}) &...
                    [events_all.(sesName).(phaseName).data.phaseCount] == str2double(ns_evt{cols.(phaseName).phaseCount+1}(ec)) &...
                    [events_all.(sesName).(phaseName).data.trial] == str2double(ns_evt{cols.(phaseName).trial+1}(ec))...
                    );
                  
                  if length(this_event) > 1
                    warning('More than one event found! Fix this script before continuing analysis.')
                    keyboard
                    %elseif isempty(this_event)
                    %  warning('No event found! Fix this script before continuing analysis.')
                    %  keyboard
                  else
                    
                    if strcmp(this_event.i_catStr,'Faces')
                      category = 'Face';
                    elseif strcmp(this_event.i_catStr,'HouseInside')
                      category = 'House';
                    end
                    
                    phaseCount = this_event.phaseCount;
                    trial = this_event.trial;
                    stimNum = this_event.stimNum;
                    i_catNum = this_event.i_catNum;
                    targ = this_event.targ;
                    spaced = this_event.spaced;
                    lag = this_event.lag;
                    expo_response = this_event.resp;
                    
                    expo_keypress = 1;
                    if expo_response == -1
                      expo_keypress = 0;
                      %             rating = '';
                      %           elseif response == 4
                      %             rating = ' VA';
                      %           elseif response == 3
                      %             rating = ' SA';
                      %           elseif response == 2
                      %             rating = ' SU';
                      %           elseif response == 1
                      %             rating = ' VU';
                      %             %else
                      %             %  rating = '';
                      %             %  keypress = 0;
                    end
                    
                    rt = this_event.rt;
                    
                    % Critical: set up the event string to match eventValues
                    %evVal = sprintf('%s%s',category,rating);
                    evVal = category;
                    
                    % find where this event type occurs in the list
                    eventNumber = find(ismember(cfg.trialdef.eventvalue,evVal));
                    if isempty(eventNumber)
                      eventNumber = -1;
                    end
                    
                    % how many time columns are there?
                    timeCols = 3;
                    % add it to the trial definition
                    this_trl = nan(1, timeCols + length(cfg.eventinfo.trl_order.(phaseName)));
                    
                    this_trl(1) = ft_event(i).sample;
                    this_trl(2) = (ft_event(i).sample + durationSamp);
                    this_trl(3) = offsetSamp;
                    
                    for to = 1:length(cfg.eventinfo.trl_order.(phaseName))
                      thisInd = find(ismember(cfg.eventinfo.trl_order.(phaseName),cfg.eventinfo.trl_order.(phaseName){to}));
                      if ~isempty(thisInd)
                        if exist(cfg.eventinfo.trl_order.(phaseName){to},'var')
                          this_trl(timeCols + thisInd) = eval(cfg.eventinfo.trl_order.(phaseName){to});
                        else
                          fprintf('variable %s does not exist!\n',cfg.eventinfo.trl_order.(phaseName){to});
                          keyboard
                        end
                      end
                    end
                    
                    % put all the trials together
                    trl = cat(1,trl,double(this_trl));
                    
                    % hardcoded old method
                    % trl = cat(1,trl,double([event(i).sample, (event(i).sample + durationSamp), offsetSamp,...
                    %   eventNumber, sesType, phaseType, phaseCount, trial, stimNum, i_catNum, targ, spaced, lag, expo_response, expo_keypress, rt]));
                    
                  end % check the event struct
                end % check the evt event
                
              case 'RESP'
                
              case 'FIXT'
                
              case 'PROM'
                
              case 'REST'
                
              case 'REND'
                
            end
            
          end
        end
        
      case {'multistudy'}
        
        %% process the multistudy phase
        
        %multistudy_events = events_all.(sesName).(phaseName).data;
        
        %trl_order_multistudy = cfg.eventinfo.trl_order.multistudy;
        
        % keep track of how many real evt events we have counted
        ec = 0;
        
        for i = 1:length(ft_event)
          if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
            % found a trigger in the evt file; increment index if value is correct.
            
            %if ~ismember(event(i).value,{'epoc'})
            if ismember(ft_event(i).value,{'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND'})
              ec = ec + 1;
            end
            
            switch ft_event(i).value
              case 'STIM'
                if strcmp(ns_evt{cols.(phaseName).isExp}(ec),'expt') && strcmpi(ns_evt{cols.(phaseName).isExp+1}(ec),'true') &&...
                    strcmp(ns_evt{cols.(phaseName).phase}(ec),'phas') && strcmp(ns_evt{cols.(phaseName).phase+1}(ec),phaseName) &&...
                    strcmp(ns_evt{cols.(phaseName).phaseCount}(ec),'pcou')
                  %&&...
                  %strcmp(ns_evt{cols.(phaseName).lag}(ec),'slag') && ~strcmp(ns_evt{cols.(phaseName).lag+1}(ec),'-1') &&...
                  %strcmp(evt{cols.(phaseName).type}(ec),'type') && strcmp(evt{cols.(phaseName).type+1}(ec),'image') &&...
                  %strcmp(evt{cols.(phaseName).icat_str}(ec),'icts')
                  
                  % Critical: set up the stimulus type, as well as the
                  % event string to match eventValues
                  if strcmp(ns_evt{cols.(phaseName).type_word+1}(ec),'word')
                    stimType = 'STUDY_WORD';
                    evVal = 'study_word';
                  elseif strcmp(ns_evt{cols.(phaseName).type_image+1}(ec),'image')
                    stimType = 'STUDY_IMAGE';
                    evVal = 'study_image';
                  end
                  
                  % find the entry in the event struct; not buffers (lag~=-1)
                  
                  this_event = events_all.(sesName).(phaseName).data(...
                    [events_all.(sesName).(phaseName).data.isExp] == 1 &...
                    ismember({events_all.(sesName).(phaseName).data.type},{stimType}) &...
                    ismember({events_all.(sesName).(phaseName).data.phaseName},{phaseName}) &...
                    [events_all.(sesName).(phaseName).data.phaseCount] == str2double(ns_evt{cols.(phaseName).phaseCount+1}(ec)) &...
                    [events_all.(sesName).(phaseName).data.trial] == str2double(ns_evt{cols.(phaseName).trial+1}(ec)) &...
                    [events_all.(sesName).(phaseName).data.lag] ~= 1 ...
                    );
                  
                  if length(this_event) > 1
                    warning('More than one event found! Fix this script before continuing analysis.')
                    keyboard
                    %elseif isempty(this_event)
                    %  warning('No event found! Fix this script before continuing analysis.')
                    %  keyboard
                  else
                    %                   if ~isempty(this_event.catStr) && strcmp(stimType,'STUDY_IMAGE')
                    %                     if strcmp(this_event.catStr,'Faces')
                    %                       category = 'Face';
                    %                     elseif strcmp(this_event.catStr,'HouseInside')
                    %                       category = 'House';
                    %                     end
                    %                   else
                    %                     category = 'Word';
                    %                   end
                    
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
                    
                    % find where this event type occurs in the list
                    eventNumber = find(ismember(cfg.trialdef.eventvalue,evVal));
                    if isempty(eventNumber)
                      eventNumber = -1;
                    end
                    
                    % how many time columns are there?
                    timeCols = 3;
                    % add it to the trial definition
                    this_trl = nan(1, timeCols + length(cfg.eventinfo.trl_order.(phaseName)));
                    
                    this_trl(1) = ft_event(i).sample;
                    this_trl(2) = (ft_event(i).sample + durationSamp);
                    this_trl(3) = offsetSamp;
                    
                    for to = 1:length(cfg.eventinfo.trl_order.(phaseName))
                      thisInd = find(ismember(cfg.eventinfo.trl_order.(phaseName),cfg.eventinfo.trl_order.(phaseName){to}));
                      if ~isempty(thisInd)
                        if exist(cfg.eventinfo.trl_order.(phaseName){to},'var')
                          this_trl(timeCols + thisInd) = eval(cfg.eventinfo.trl_order.(phaseName){to});
                        else
                          fprintf('variable %s does not exist!\n',cfg.eventinfo.trl_order.(phaseName){to});
                          keyboard
                        end
                      end
                    end
                    
                    % put all the trials together
                    trl = cat(1,trl,double(this_trl));
                    
                  end % check the event struct
                end % check the evt event
                
              case 'RESP'
                
              case 'FIXT'
                
              case 'PROM'
                
              case 'REST'
                
              case 'REND'
                
            end
            
          end
        end
        
    end % switch phase
  end % pha
end % ses

%% old

%           if strcmp(evt{cols.icat_str+1}(ec),'Faces')
%             category = 'Face';
%           elseif strcmp(evt{cols.icat_str+1}(ec),'HouseInside')
%             category = 'House';
%           end
%           catNum = str2double(evt{cols.icat_num+1}(ec));
%
%           targ = str2double(evt{cols.targ+1}(ec));
%
%           if strcmp(evt{cols.spaced+1}(ec),'true')
%             spaced = 1;
%           elseif strcmp(evt{cols.spaced+1}(ec),'false')
%             spaced = 0;
%           end
%
%           lag = str2double(evt{cols.lag+1}(ec));
%
%           if strcmp(evt{cols.resp_str+1}(ec),'v_appeal')
%             response = 4;
%             rating = ' VA';
%           elseif strcmp(evt{cols.resp_str+1}(ec),'s_appeal')
%             response = 3;
%             rating = ' SA';
%           elseif strcmp(evt{cols.resp_str+1}(ec),'s_unappeal')
%             response = 2;
%             rating = ' SU';
%           elseif strcmp(evt{cols.resp_str+1}(ec),'v_unappeal')
%             response = 1;
%             rating = ' VU';
%           else
%             response = -1;
%             rating = '';
%           end
%
%           if strcmp(evt{cols.keypress+1}(ec),'true')
%             keypress = 1;
%           elseif strcmp(evt{cols.keypress+1}(ec),'false')
%             keypress = 0;
%           end
%
%           rt = str2double(evt{cols.rt+1}(ec));
