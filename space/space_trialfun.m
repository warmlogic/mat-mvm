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

%% process the exposure phase

sesName = 'oneDay';
sesType = find(ismember(sessionNames,sesName));
% sesType = find(ismember(expParam.sesTypes,sesName));
phaseName = 'expo';
phaseType = find(ismember(phaseNames,phaseName));
% phaseType = find(ismember(expParam.session.(sesName).phases,phaseName));

expo_events = events_all.(sesName).(phaseName).data;

cols = [];
cols.sub = 7;
cols.ses = 9;
cols.phase = 11;
cols.phasenum = 13;
cols.isExp = 15;
cols.trial = 17;
cols.type = 19;
cols.inum = 23;
cols.icat_str = 25;
cols.icat_num = 27;
cols.targ = 29;
cols.spaced = 31;
cols.lag = 33;
cols.resp_str = 35;
cols.rt = 39;
cols.keypress = 41;

expo_trl_order = cfg.eventinfo.expo_trl_order;

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
        if strcmp(evt{cols.isExp}(ec),'expt') && strcmpi(evt{cols.isExp+1}(ec),'true') &&...
            strcmp(evt{cols.phase}(ec),'phas') && strcmp(evt{cols.phase+1}(ec),phaseName) &&...
            strcmp(evt{cols.type}(ec),'type') && strcmp(evt{cols.type+1}(ec),'image') &&...
            strcmp(evt{cols.phasenum}(ec),'pcou') &&...
            strcmp(evt{cols.icat_str}(ec),'icts')
          
          % find the entry in the event struct
          this_event = expo_events(...
            [expo_events.isExp] == 1 &...
            ismember({expo_events.type},{'EXPO_IMAGE'}) &...
            ismember({expo_events.phaseName},{phaseName}) &...
            [expo_events.phaseCount] == str2double(evt{cols.phasenum+1}(ec)) &...
            [expo_events.trial] == str2double(evt{cols.trial+1}(ec))...
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
          this_trl = nan(1, 3 + length(expo_trl_order));
          
          this_trl(1) = event(i).sample;
          this_trl(2) = (event(i).sample + durationSamp);
          this_trl(3) = offsetSamp;
          
          for to = 1:length(expo_trl_order)
            thisInd = find(ismember(expo_trl_order,expo_trl_order{to}));
            if ~isempty(thisInd)
              if exist(expo_trl_order{to},'var')
                this_trl(thisInd) = eval(expo_trl_order{to});
              else
                fprintf('variable %s does not exist!\n',expo_trl_order{to});
                keyboard
              end
            end
          end
          
          % put all the trials together
          trl = cat(1,trl,double(this_trl));
          
          % trl = cat(1,trl,double([event(i).sample, (event(i).sample + durationSamp), offsetSamp,...
          %   eventNumber, sesType, phaseType, phaseCount, trial, stimNum, i_catNum, targ, spaced, lag, expo_response, expo_keypress, rt]));
          %
          % %{'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'stimNum', 'i_catNum', 'targ', 'spaced', 'lag', 'expo_response', 'expo_keypress', 'rt'};

        end
        
      case 'RESP'
        
      case 'FIXT'
        
      case 'PROM'
        
      case 'REST'
        
      case 'REND'
        
    end

  end
end

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
