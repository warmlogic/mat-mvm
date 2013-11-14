function trl = space_trialfun(cfg)

% operates using Net Station evt files

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% offset should be negative
offsetSamp = round(-cfg.trialdef.prestim*hdr.Fs);
% duration should be 1 sample less than the whole length of an event
durationSamp = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
% TODO: should this be ceil instead of round?

evt = cfg.eventinfo.evt;
events_all = cfg.eventinfo.events;
sessionNames = cfg.eventinfo.sessionNames;
phaseNames = cfg.eventinfo.phaseNames;
% expParam = cfg.eventinfo.expParam;

% initialize the trl matrix
trl = [];

for ses = 1:length(sessionNames)
  sesType = find(ismember(sessionNames,sessionNames{ses}));
  
  for pha = 1:length(phaseNames)
    phaseName = phaseNames{sesType}{pha};
    phaseType = find(ismember(phaseNames{sesType},phaseNames{sesType}{pha}));
    
    switch phaseName
      
      case {'expo'}
        %% process the exposure phase
        expo_events = events_all.(sessionNames{ses}).(phaseName).data;
        
        cols_expo = [];
        cols_expo.sub = 7;
        cols_expo.ses = 9;
        cols_expo.phase = 11;
        cols_expo.phasenum = 13;
        cols_expo.isExp = 15;
        cols_expo.trial = 17;
        cols_expo.type = 19;
        cols_expo.inum = 23;
        cols_expo.icat_str = 25;
        cols_expo.icat_num = 27;
        cols_expo.targ = 29;
        cols_expo.spaced = 31;
        cols_expo.lag = 33;
        cols_expo.resp_str = 35;
        cols_expo.rt = 39;
        cols_expo.keypress = 41;
        
        trl_order_expo = cfg.eventinfo.trl_order.expo;
        
        % keep track of how many real evt events we have counted
        ec = 0;
        
        for i = 1:length(event)
          if strcmp(event(i).type,cfg.trialdef.eventtype)
            % found a trigger in the evt file; increment index if value is correct.
            
            %if ~ismember(event(i).value,{'epoc'})
            if ismember(event(i).value,{'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND'})
              ec = ec + 1;
            end
            
            switch event(i).value
              case 'STIM'
                if strcmp(evt{cols_expo.isExp}(ec),'expt') && strcmpi(evt{cols_expo.isExp+1}(ec),'true') &&...
                    strcmp(evt{cols_expo.phase}(ec),'phas') && strcmp(evt{cols_expo.phase+1}(ec),phaseName) &&...
                    strcmp(evt{cols_expo.type}(ec),'type') && strcmp(evt{cols_expo.type+1}(ec),'image') &&...
                    strcmp(evt{cols_expo.phasenum}(ec),'pcou') &&...
                    strcmp(evt{cols_expo.icat_str}(ec),'icts')
                  
                  % find the entry in the event struct
                  this_event = expo_events(...
                    [expo_events.isExp] == 1 &...
                    ismember({expo_events.type},{'EXPO_IMAGE'}) &...
                    ismember({expo_events.phaseName},{phaseName}) &...
                    [expo_events.phaseCount] == str2double(evt{cols_expo.phasenum+1}(ec)) &...
                    [expo_events.trial] == str2double(evt{cols_expo.trial+1}(ec))...
                    );
                  
                  if length(this_event) ~= 1
                    warning('More than one event found! Fix this script before continuing analysis.')
                    keyboard
                  end
                  
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
                  
                  % add it to the trial definition
                  this_trl = nan(1, 3 + length(trl_order_expo));
                  
                  this_trl(1) = event(i).sample;
                  this_trl(2) = (event(i).sample + durationSamp);
                  this_trl(3) = offsetSamp;
                  
                  for to = 1:length(trl_order_expo)
                    thisInd = find(ismember(trl_order_expo,trl_order_expo{to}));
                    if ~isempty(thisInd)
                      if exist(trl_order_expo{to},'var')
                        this_trl(thisInd) = eval(trl_order_expo{to});
                      else
                        fprintf('variable %s does not exist!\n',trl_order_expo{to});
                        keyboard
                      end
                    end
                  end
                  
                  % put all the trials together
                  trl = cat(1,trl,double(this_trl));
                  
                  % hardcoded old method
                  % trl = cat(1,trl,double([event(i).sample, (event(i).sample + durationSamp), offsetSamp,...
                  %   eventNumber, sesType, phaseType, phaseCount, trial, stimNum, i_catNum, targ, spaced, lag, expo_response, expo_keypress, rt]));
                  
                end
                
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
        multistudy_events = events_all.(sesName).(phaseName).data;
        
        cols_ms = [];
        cols_ms.sub = 7;
        cols_ms.ses = 9;
        cols_ms.phase = 11;
        cols_ms.phasenum = 13;
        cols_ms.isExp = 15;
        cols_ms.trial = 17;
        cols_ms.type = 19;
        cols_ms.inum = 23;
        cols_ms.icat_str = 25;
        cols_ms.icat_num = 27;
        cols_ms.targ = 29;
        cols_ms.spaced = 31;
        cols_ms.lag = 33;
        cols_ms.resp_str = 35;
        cols_ms.rt = 39;
        cols_ms.keypress = 41;
        
        trl_order_multistudy = cfg.eventinfo.trl_order.multistudy;
        
        % keep track of how many real evt events we have counted
        ec = 0;
        
        for i = 1:length(event)
          if strcmp(event(i).type,cfg.trialdef.eventtype)
            % found a trigger in the evt file; increment index if value is correct.
            
            %if ~ismember(event(i).value,{'epoc'})
            if ismember(event(i).value,{'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND'})
              ec = ec + 1;
            end
            
            switch event(i).value
              case 'STIM'
                if strcmp(evt{cols_ms.isExp}(ec),'expt') && strcmpi(evt{cols_ms.isExp+1}(ec),'true') &&...
                    strcmp(evt{cols_ms.phase}(ec),'phas') && strcmp(evt{cols_ms.phase+1}(ec),phaseName) &&...
                    strcmp(evt{cols_ms.type}(ec),'type') && strcmp(evt{cols_ms.type+1}(ec),'image') &&...
                    strcmp(evt{cols_ms.phasenum}(ec),'pcou') &&...
                    strcmp(evt{cols_ms.icat_str}(ec),'icts')
                  
                  % find the entry in the event struct
                  this_event = multistudy_events(...
                    [multistudy_events.isExp] == 1 &...
                    ismember({multistudy_events.type},{'EXPO_IMAGE'}) &...
                    ismember({multistudy_events.phaseName},{phaseName}) &...
                    [multistudy_events.phaseCount] == str2double(evt{cols_ms.phasenum+1}(ec)) &...
                    [multistudy_events.trial] == str2double(evt{cols_ms.trial+1}(ec))...
                    );
                  
                  if length(this_event) ~= 1
                    warning('More than one event found! Fix this script before continuing analysis.')
                    keyboard
                  end
                  
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
                  ms_response = this_event.resp;
                  
                  ms_keypress = 1;
                  if ms_response == -1
                    ms_keypress = 0;
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
                  
                  % add it to the trial definition
                  this_trl = nan(1, 3 + length(trl_order_multistudy));
                  
                  this_trl(1) = event(i).sample;
                  this_trl(2) = (event(i).sample + durationSamp);
                  this_trl(3) = offsetSamp;
                  
                  for to = 1:length(trl_order_multistudy)
                    thisInd = find(ismember(trl_order_multistudy,trl_order_multistudy{to}));
                    if ~isempty(thisInd)
                      if exist(trl_order_multistudy{to},'var')
                        this_trl(thisInd) = eval(trl_order_multistudy{to});
                      else
                        fprintf('variable %s does not exist!\n',trl_order_multistudy{to});
                        keyboard
                      end
                    end
                  end
                  
                  % put all the trials together
                  trl = cat(1,trl,double(this_trl));
                  
                  % hardcoded old method
                  % trl = cat(1,trl,double([event(i).sample, (event(i).sample + durationSamp), offsetSamp,...
                  %   eventNumber, sesType, phaseType, phaseCount, trial, stimNum, i_catNum, targ, spaced, lag, ms_response, ms_keypress, rt]));
                  
                end
                
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
