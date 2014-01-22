function trl = nemo_trialfun(cfg)

% operates using Net Station evt files and event structs

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in external data, if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg.eventinfo.useMetadata
  md = cfg.eventinfo.metadata;
  
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
  end
  
  if ismember('expParam', md.types)
    error('loading experimentParams.mat is not yet supported');
  end
end

% get the header and event information
ft_hdr = ft_read_header(cfg.dataset);
ft_event = ft_read_event(cfg.dataset);

ns_evt = cfg.eventinfo.ns_evt;
events_all = cfg.eventinfo.events;
% expParam = cfg.eventinfo.expParam;

%triggers = {'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND', 'DIN '};
triggers = unique(ns_evt{1});

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

%% go through the events

for pha = 1:length(cfg.eventinfo.phaseNames)
%   phaseName = cfg.eventinfo.phaseNames{pha};
%   phaseType = find(ismember(cfg.eventinfo.phaseNames{sesType},phaseName));
%   phaseType = pha;
  
  switch phaseName
    
    case {'TC_NEMO_fN400study'}
      %% process the study phase
      
      fprintf('Processing %s...\n',phaseName);
      
      % keep track of how many real evt events we have counted
      ec = 0;
      
      fprintf('%s NS event flag count: %s',phaseName,repmat(' ',1,length(num2str(length(ft_event)))));
      
      for i = 1:length(ft_event)
        fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
        
        if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
          % found a trigger in the EEG file events; increment index if
          % value is correct.
          
          %if ~ismember(event(i).value,{'epoc'})
          if ismember(ft_event(i).value,triggers)
            ec = ec + 1;
          end
          
          switch ft_event(i).value
            case 'TRSP'
              
              % hack: 2 special cases: evt events occurred at the same
              % sample (but possibly at different MS). Since the evt
              % records MS and FieldTrip events records samples, this can
              % cause some screwy things to happen. Both events always
              % exist in the evt file; however, when FT reads events, it
              % seems to respect events with different codes, but it
              % ignores one of the two with the same code. In the former
              % case, sometimes they are in a slightly different order in
              % the evt file compared to the events that FT reads due to
              % two events having the same sample time but different MS
              % times, so ec needs to get reset to its previous state. In
              % the latter case, since FT completely skips duplicate events
              % at the same sample, we simply need to increment ec by 1.
              ec_add = 0;
              if ~strcmp(ns_evt{1}(ec),ft_event(i).value)% && (strcmp(ns_evt{1}(ec-1),ft_event(i).value) || strcmp(ns_evt{1}(ec+1),ft_event(i).value))
                this_time_ms_str = ns_evt{5}(ec);
                this_time_ms = (str2double(this_time_ms_str{1}(2:3)) * 60 * 60 * 1000) + (str2double(this_time_ms_str{1}(5:6)) * 60 * 1000) + (str2double(this_time_ms_str{1}(8:9)) * 1000) + (str2double(this_time_ms_str{1}(11:13)));
                this_time_samp = fix((this_time_ms / 1000) * ft_hdr.Fs);
                prev_time_ms_str = ns_evt{5}(ec-1);
                prev_time_ms = (str2double(prev_time_ms_str{1}(2:3)) * 60 * 60 * 1000) + (str2double(prev_time_ms_str{1}(5:6)) * 60 * 1000) + (str2double(prev_time_ms_str{1}(8:9)) * 1000) + (str2double(prev_time_ms_str{1}(11:13)));
                prev_time_samp = fix((prev_time_ms / 1000) * ft_hdr.Fs);
                
                if this_time_samp == prev_time_samp || abs(this_time_ms - prev_time_ms) == 1
                  % events occurred at the same sample or one ms apart
                  
                  if strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
                    % don't put ec back in its prior state if the event
                    % codes were the same
                    ec = ec + 1;
                  elseif ~strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
                    % put ec back in its prior state if the event codes
                    % were not the same
                    if strcmp(ns_evt{1}(ec-1),ft_event(i).value)
                      ec = ec - 1;
                      ec_add = 1;
                    elseif strcmp(ns_evt{1}(ec+1),ft_event(i).value)
                      ec = ec + 1;
                      ec_add = -1;
                    end
                  end
                end
              end
              
              if strcmp(ns_evt{1}(ec),ft_event(i).value)
                
                % set column types because Net Station evt files can vary
                ns_evt_cols = {};
                for ns = 1:size(ns_evt,2)
                  ns_evt_cols = cat(1,ns_evt_cols,ns_evt{ns}(ec));
                end
                cols.(phaseName).isExp = find(strcmp(ns_evt_cols,'expt'));
                cols.(phaseName).phase = find(strcmp(ns_evt_cols,'exid'));     %string but the only phase identifier
                if isempty(cols.(phaseName).isExp) || isempty(cols.(phaseName).phase)
                  keyboard
                end
                
                if strcmp(ns_evt{cols.(phaseName).phase+1}(ec),phaseName) &&...
                    strcmpi(ns_evt{cols.(phaseName).isExp+1}(ec),'true')
                              
                 % set column types because Net Station evt files can vary
                  %cols.(phaseName).phaseCount = find(strcmp(ns_evt_cols,'pcou'));
                  cols.(phaseName).trial = find(strcmp(ns_evt_cols,'trl#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cel#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'obs#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pos#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rsp#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'eval'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rtim'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'prim')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'targ')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'fstg'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'bstg'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'relt')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cnct'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'frqt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'lent'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phot'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'sylt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'ornt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phnt'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'ccat')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'aoat'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'bfqt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'famt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'imgt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cmnt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pmnt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'tfrt'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'valt')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'rspt')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'widt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cncp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'frqp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'lenp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phop'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'sylp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'ornp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phnp'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'ccap')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'aoap'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'bfqp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'famp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'imgp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cmnp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pmnp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'tfrp'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'valp')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'rspp')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'widp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pair'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rnum'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'mody')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'task')); %string
                  if isempty(cols.(phaseName).phaseCount) || isempty(cols.(phaseName).trial) || isempty(cols.(phaseName).type)
                    keyboard
                  end         
                 
                  % Critical: set up the stimulus type, as well as the
                  % event string to match eventValues
                  if strcmp(ns_evt{cols.(phaseName).type+1}(ec),'widp')
                    stimType = 'study_prime';
                    evVal = 'widp';
                  elseif strcmp(ns_evt{cols.(phaseName).type+1}(ec),'widt')
                    stimType = 'study_target';
                    evVal = 'widt';
                  end
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
                  
%                   % find the entry in the event struct
%                   this_event = events_all.(sesName).(phaseName).data(...
%                     [events_all.(phaseName).data.isExp] == 1 &...
%                     ismember({events_all.(phaseName).data.phaseName},{phaseName}) &...
%                     [events_all.(phaseName).data.phaseCount] == str2double(ns_evt{cols.(phaseName).phaseCount+1}(ec)) &...
%                     ismember({events_all.(phaseName).data.type},{stimType}) &...
%                     [events_all.(phaseName).data.trial] == str2double(ns_evt{cols.(phaseName).trial+1}(ec)) ...
%                     );
%                   
%                   if length(this_event) > 1
%                     warning('More than one event found! Fix this script before continuing analysis.')
%                     keyboard
%                   elseif isempty(this_event)
%                     warning('No event found! Fix this script before continuing analysis.')
%                     keyboard
%                   elseif length(this_event) == 1
%                     
%                     phaseCount = this_event.phaseCount;
%                     trial = this_event.trial;
%                     stimNum = this_event.stimNum;
%                     i_catNum = this_event.i_catNum;
%                     targ = this_event.targ;
%                     spaced = this_event.spaced;
%                     lag = this_event.lag;
%                     expo_response = this_event.resp;
%                     cr_recog_acc = this_event.cr_recog_acc;
%                     cr_recall_resp = this_event.cr_recall_resp;
%                     cr_recall_spellCorr = this_event.cr_recall_spellCorr;
%                     
%                     rt = this_event.rt;
%                     
%                     % add it to the trial definition
%                     this_trl = trl_ini;
%                     
%                     % get the time of this event
%                     this_sample = ft_event(i).sample;
%                     
% %                     % if we're using the photodiode DIN and we find one
% %                     % within the threshold, replace the current sample time
% %                     % with that of the DIN
% %                     if cfg.eventinfo.usePhotodiodeDIN
% %                       photodiodeDIN_thresholdSamp = round((cfg.eventinfo.photodiodeDIN_thresholdMS / 1000) * ft_hdr.Fs);
% %                       if strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str)
% %                         % if there is a DIN before and after the stim, pick
% %                         % the closer one
% %                         preDiff = (ft_event(i-1).sample - this_sample);
% %                         postDiff = (ft_event(i+1).sample - this_sample);
% %                         
% %                         if preDiff < 0 && abs(preDiff) <= photodiodeDIN_thresholdSamp
% %                           preFlag = true;
% %                         else
% %                           preFlag = false;
% %                         end
% %                         if postDiff <= photodiodeDIN_thresholdSamp
% %                           postFlag = true;
% %                         else
% %                           postFlag = false;
% %                         end
% %                         
% %                         if preFlag && ~postFlag
% %                           % only the pre-DIN makes sense
% %                           this_sample = ft_event(i-1).sample;
% %                         elseif ~preFlag && postFlag
% %                           % only the post-DIN makes sense
% %                           this_sample = ft_event(i+1).sample;
% %                         elseif preFlag && postFlag
% %                           % choose the smaller one
% %                           if abs(preDiff) < abs(postDiff)
% %                             this_sample = ft_event(i-1).sample;
% %                           elseif abs(preDiff) > abs(postDiff)
% %                             this_sample = ft_event(i+1).sample;
% %                           elseif abs(preDiff) == abs(postDiff)
% %                             keyboard
% %                           end
% %                         end
% %                       elseif strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && ~strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i+1).sample - this_sample) <= photodiodeDIN_thresholdSamp
% %                         this_sample = ft_event(i+1).sample;
% %                       elseif strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i-1).sample - this_sample) < 0 && abs(ft_event(i-1).sample - this_sample) <= photodiodeDIN_thresholdSamp
% %                         % apparently the ethernet tags can be delayed
% %                         % enough that the DIN shows up first
% %                         this_sample = ft_event(i-1).sample;
% %                       end
% %                     end
%                     
%                     % prestimulus sample
%                     this_trl(1) = this_sample + prestimSamp;
%                     % poststimulus sample
%                     this_trl(2) = this_sample + poststimSamp;
%                     % offset in samples
%                     this_trl(3) = prestimSamp;
%                     
%                     for to = 1:length(trl_order)
%                       thisInd = find(ismember(trl_order,trl_order{to}));
%                       if ~isempty(thisInd)
%                         if exist(trl_order{to},'var')
%                           this_trl(timeCols + thisInd) = eval(trl_order{to});
%                         else
%                           fprintf('variable %s does not exist!\n',trl_order{to});
%                           keyboard
%                         end
%                       end
%                     end
%                     
%                     % put all the trials together
%                     trl = cat(1,trl);
%                     
%                   end % check the event struct
                end % check the evt event
                
              else
                % the count is off? EEG and evt events don't line up
                keyboard
              end
              
              % put ec back in its prior state
              ec = ec + ec_add;
              
            case 'bgin'
              
            case 'prm+'
              
            case 'resp'
              
            case 'fix+'
              
            case 'trg+'
                            
          end
          
        end
      end % for
      fprintf('\n');
      
    case {'TC_NEMO_AO'}
      
      %% process the auditory oddball phase
      
      fprintf('Processing %s...\n',phaseName);
      
      % keep track of how many real evt events we have counted
      ec = 0;
      
      fprintf('%s NS event flag count: %s',phaseName,repmat(' ',1,length(num2str(length(ft_event)))));
      
      for i = 1:length(ft_event)
        fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
        
        if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
          % found a trigger in the evt file; increment index if value is correct.
          
          %if ~ismember(event(i).value,{'epoc'})
          if ismember(ft_event(i).value,triggers)
            ec = ec + 1;
          end
          
          switch ft_event(i).value
            case 'TRSP'
              
              % hack: 2 special cases: evt events occurred at the same
              % sample (but possibly at different MS). Since the evt
              % records MS and FieldTrip events records samples, this can
              % cause some screwy things to happen. Both events always
              % exist in the evt file; however, when FT reads events, it
              % seems to respect events with different codes, but it
              % ignores one of the two with the same code. In the former
              % case, sometimes they are in a slightly different order in
              % the evt file compared to the events that FT reads due to
              % two events having the same sample time but different MS
              % times, so ec needs to get reset to its previous state. In
              % the latter case, since FT completely skips duplicate events
              % at the same sample, we simply need to increment ec by 1.
              ec_add = 0;
              if ~strcmp(ns_evt{1}(ec),ft_event(i).value)% && (strcmp(ns_evt{1}(ec-1),ft_event(i).value) || strcmp(ns_evt{1}(ec+1),ft_event(i).value))
                this_time_ms_str = ns_evt{5}(ec);
                this_time_ms = (str2double(this_time_ms_str{1}(2:3)) * 60 * 60 * 1000) + (str2double(this_time_ms_str{1}(5:6)) * 60 * 1000) + (str2double(this_time_ms_str{1}(8:9)) * 1000) + (str2double(this_time_ms_str{1}(11:13)));
                this_time_samp = fix((this_time_ms / 1000) * ft_hdr.Fs);
                prev_time_ms_str = ns_evt{5}(ec-1);
                prev_time_ms = (str2double(prev_time_ms_str{1}(2:3)) * 60 * 60 * 1000) + (str2double(prev_time_ms_str{1}(5:6)) * 60 * 1000) + (str2double(prev_time_ms_str{1}(8:9)) * 1000) + (str2double(prev_time_ms_str{1}(11:13)));
                prev_time_samp = fix((prev_time_ms / 1000) * ft_hdr.Fs);
                
                if this_time_samp == prev_time_samp || abs(this_time_ms - prev_time_ms) == 1
                  % events occurred at the same sample or one ms apart
                  
                  if strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
                    % don't put ec back in its prior state if the event
                    % codes were the same
                    ec = ec + 1;
                  elseif ~strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
                    % put ec back in its prior state if the event codes
                    % were not the same
                    if strcmp(ns_evt{1}(ec-1),ft_event(i).value)
                      ec = ec - 1;
                      ec_add = 1;
                    elseif strcmp(ns_evt{1}(ec+1),ft_event(i).value)
                      ec = ec + 1;
                      ec_add = -1;
                    end
                  end
                end
              end
              
              if strcmp(ns_evt{1}(ec),ft_event(i).value)
                
                % set column types because Net Station evt files can vary
                ns_evt_cols = {};
                for ns = 1:size(ns_evt,2)
                  ns_evt_cols = cat(1,ns_evt_cols,ns_evt{ns}(ec));
                end
                cols.(phaseName).isExp = find(strcmp(ns_evt_cols,'expt'));
                cols.(phaseName).phase = find(strcmp(ns_evt_cols,'exid'));
                if isempty(cols.(phaseName).isExp) || isempty(cols.(phaseName).phase)
                  keyboard
                end
                
                if strcmpi(ns_evt{cols.(phaseName).isExp+1}(ec),'true') &&...
                    strcmp(ns_evt{cols.(phaseName).phase+1}(ec),phaseName)
                  
                 % set column types because Net Station evt files can vary
%                 cols.(phaseName).phaseCount = find(strcmp(ns_evt_cols,'pcou'));
                  cols.(phaseName).trial = find(strcmp(ns_evt_cols,'trl#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rsp#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'eval'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rtim'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'argu'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pos#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'obs#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rsp+'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cel#'));
                  if isempty(cols.(phaseName).phaseCount) || isempty(cols.(phaseName).trial) || isempty(cols.(phaseName).type)
                    keyboard
                  end            
                
                  % original set column
                  cols.(phaseName).phaseCount = find(strcmp(ns_evt_cols,'pcou'));
                  cols.(phaseName).trial = find(strcmp(ns_evt_cols,'trln'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'type'));
                  if isempty(cols.(phaseName).phaseCount) || isempty(cols.(phaseName).trial) || isempty(cols.(phaseName).type)
                    keyboard
                  end
                  
                  % Critical: set up the stimulus type, as well as the
                  % event string to match eventValues
                  if strcmp(ns_evt{cols.(phaseName).type+1}(ec),'11')
                    stimType = 'ao_standard';
                    evVal = 'cel#';
                  elseif strcmp(ns_evt{cols.(phaseName).type+1}(ec),'12')
                    stimType = 'ao_target';
                    evVal = 'cel#';
                  end
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
                    [events_all.(sesName).(phaseName).data.phaseCount] == str2double(ns_evt{cols.(phaseName).phaseCount+1}(ec)) &...
                    ismember({events_all.(sesName).(phaseName).data.type},{stimType}) &...
                    [events_all.(sesName).(phaseName).data.trial] == str2double(ns_evt{cols.(phaseName).trial+1}(ec)) ...
                    );
                  % exclude  buffers (lag~=-1)
                  %  [events_all.(sesName).(phaseName).data.lag] ~= -1 ...
                  
                  if length(this_event) > 1
                    warning('More than one event found! Fix this script before continuing analysis.')
                    keyboard
                  elseif isempty(this_event)
                    warning('No event found! Fix this script before continuing analysis.')
                    keyboard
                  elseif length(this_event) == 1
                    
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
                    
                    % add it to the trial definition
                    this_trl = trl_ini;
                    
                    % get the time of this event
                    this_sample = ft_event(i).sample;
                    
%                     % if we're using the photodiode DIN and we find one
%                     % within the threshold, replace the current sample time
%                     % with that of the DIN
%                     if cfg.eventinfo.usePhotodiodeDIN
%                       photodiodeDIN_thresholdSamp = round((cfg.eventinfo.photodiodeDIN_thresholdMS / 1000) * ft_hdr.Fs);
%                       if strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str)
%                         % if there is a DIN before and after the stim, pick
%                         % the closer one
%                         preDiff = (ft_event(i-1).sample - this_sample);
%                         postDiff = (ft_event(i+1).sample - this_sample);
%                         
%                         if preDiff < 0 && abs(preDiff) <= photodiodeDIN_thresholdSamp
%                           preFlag = true;
%                         else
%                           preFlag = false;
%                         end
%                         if postDiff <= photodiodeDIN_thresholdSamp
%                           postFlag = true;
%                         else
%                           postFlag = false;
%                         end
%                         
%                         if preFlag && ~postFlag
%                           % only the pre-DIN makes sense
%                           this_sample = ft_event(i-1).sample;
%                         elseif ~preFlag && postFlag
%                           % only the post-DIN makes sense
%                           this_sample = ft_event(i+1).sample;
%                         elseif preFlag && postFlag
%                           % choose the smaller one
%                           if abs(preDiff) < abs(postDiff)
%                             this_sample = ft_event(i-1).sample;
%                           elseif abs(preDiff) > abs(postDiff)
%                             this_sample = ft_event(i+1).sample;
%                           elseif abs(preDiff) == abs(postDiff)
%                             keyboard
%                           end
%                         end
%                       elseif strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && ~strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i+1).sample - this_sample) <= photodiodeDIN_thresholdSamp
%                         this_sample = ft_event(i+1).sample;
%                       elseif strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i-1).sample - this_sample) < 0 && abs(ft_event(i-1).sample - this_sample) <= photodiodeDIN_thresholdSamp
%                         % apparently the ethernet tags can be delayed
%                         % enough that the DIN shows up first
%                         this_sample = ft_event(i-1).sample;
%                       end
%                     end
%                     
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
                    
                  end % check the event struct
                end % check the evt event
                
              else
                % the count is off? EEG and evt events don't line up
                keyboard
              end
              
              % put ec back in its prior state
              ec = ec + ec_add;
              
            case 'bgin'
              
            case 'stm+'
              
            case 'resp'
             
          end
          
        end
      end % for
      fprintf('\n');
      
    case {'TC_NEMO_fN400test'}
      
      %% process the test phase
      
      fprintf('Processing %s...\n',phaseName);
      
      % keep track of how many real evt events we have counted
      ec = 0;
      
      fprintf('%s NS event flag count: %s',phaseName,repmat(' ',1,length(num2str(length(ft_event)))));
      
      for i = 1:length(ft_event)
        fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
        
        if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
          % found a trigger in the evt file; increment index if value is correct.
          
          %if ~ismember(event(i).value,{'epoc'})
          if ismember(ft_event(i).value,triggers)
            ec = ec + 1;
          end
          
          switch ft_event(i).value
            case 'TRSP'
              
              % hack: 2 special cases: evt events occurred at the same
              % sample (but possibly at different MS). Since the evt
              % records MS and FieldTrip events records samples, this can
              % cause some screwy things to happen. Both events always
              % exist in the evt file; however, when FT reads events, it
              % seems to respect events with different codes, but it
              % ignores one of the two with the same code. In the former
              % case, sometimes they are in a slightly different order in
              % the evt file compared to the events that FT reads due to
              % two events having the same sample time but different MS
              % times, so ec needs to get reset to its previous state. In
              % the latter case, since FT completely skips duplicate events
              % at the same sample, we simply need to increment ec by 1.
              ec_add = 0;
              if ~strcmp(ns_evt{1}(ec),ft_event(i).value)% && (strcmp(ns_evt{1}(ec-1),ft_event(i).value) || strcmp(ns_evt{1}(ec+1),ft_event(i).value))
                this_time_ms_str = ns_evt{5}(ec);
                this_time_ms = (str2double(this_time_ms_str{1}(2:3)) * 60 * 60 * 1000) + (str2double(this_time_ms_str{1}(5:6)) * 60 * 1000) + (str2double(this_time_ms_str{1}(8:9)) * 1000) + (str2double(this_time_ms_str{1}(11:13)));
                this_time_samp = fix((this_time_ms / 1000) * ft_hdr.Fs);
                prev_time_ms_str = ns_evt{5}(ec-1);
                prev_time_ms = (str2double(prev_time_ms_str{1}(2:3)) * 60 * 60 * 1000) + (str2double(prev_time_ms_str{1}(5:6)) * 60 * 1000) + (str2double(prev_time_ms_str{1}(8:9)) * 1000) + (str2double(prev_time_ms_str{1}(11:13)));
                prev_time_samp = fix((prev_time_ms / 1000) * ft_hdr.Fs);
                
                if this_time_samp == prev_time_samp || abs(this_time_ms - prev_time_ms) == 1
                  % events occurred at the same sample or one ms apart
                  
                  if strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
                    % don't put ec back in its prior state if the event
                    % codes were the same
                    ec = ec + 1;
                  elseif ~strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
                    % put ec back in its prior state if the event codes
                    % were not the same
                    if strcmp(ns_evt{1}(ec-1),ft_event(i).value)
                      ec = ec - 1;
                      ec_add = 1;
                    elseif strcmp(ns_evt{1}(ec+1),ft_event(i).value)
                      ec = ec + 1;
                      ec_add = -1;
                    end
                  end
                end
              end
              
              if strcmp(ns_evt{1}(ec),ft_event(i).value)
                
                % set column types because Net Station evt files can vary
                ns_evt_cols = {};
                for ns = 1:size(ns_evt,2)
                  ns_evt_cols = cat(1,ns_evt_cols,ns_evt{ns}(ec));
                end
                cols.(phaseName).isExp = find(strcmp(ns_evt_cols,'expt'));
                cols.(phaseName).phase = find(strcmp(ns_evt_cols,'phas'));
                if isempty(cols.(phaseName).isExp) || isempty(cols.(phaseName).phase)
                  keyboard
                end
                
                if strcmpi(ns_evt{cols.(phaseName).isExp+1}(ec),'true') &&...
                    strcmp(ns_evt{cols.(phaseName).phase+1}(ec),phaseName)
                  
                  % set column types because Net Station evt files can vary
                ns_evt_cols = {};
                for ns = 1:size(ns_evt,2)
                  ns_evt_cols = cat(1,ns_evt_cols,ns_evt{ns}(ec));
                end
                cols.(phaseName).isExp = find(strcmp(ns_evt_cols,'expt'));
                cols.(phaseName).phase = find(strcmp(ns_evt_cols,'exid'));     %string but the only phase identifier
                if isempty(cols.(phaseName).isExp) || isempty(cols.(phaseName).phase)
                  keyboard
                end
                
                if strcmp(ns_evt{cols.(phaseName).phase+1}(ec),phaseName) &&...
                    strcmpi(ns_evt{cols.(phaseName).isExp+1}(ec),'true')
                              
                 % set column types because Net Station evt files can vary
                  %cols.(phaseName).phaseCount = find(strcmp(ns_evt_cols,'pcou'));
                  cols.(phaseName).trial = find(strcmp(ns_evt_cols,'trl#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cel#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'obs#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pos#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rsp#'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'eval'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rtim'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'prim')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'targ')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'fstg'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'bstg'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'relt')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cnct'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'frqt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'lent'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phot'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'sylt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'ornt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phnt'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'ccat')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'aoat'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'bfqt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'famt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'imgt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cmnt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pmnt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'tfrt'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'valt')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'rspt')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'widt'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cncp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'frqp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'lenp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phop'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'sylp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'ornp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'phnp'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'ccap')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'aoap'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'bfqp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'famp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'imgp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'cmnp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pmnp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'tfrp'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'valp')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'rspp')); %string
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'widp'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'pair'));
                  cols.(phaseName).type = find(strcmp(ns_evt_cols,'rnum'));
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'mody')); %string
                  %cols.(phaseName).type = find(strcmp(ns_evt_cols,'task')); %string
                  if isempty(cols.(phaseName).phaseCount) || isempty(cols.(phaseName).trial)
                    keyboard
                  end
                  
                  % Critical: set up the stimulus type, as well as the
                  % event string to match eventValues
                  stimType = 'test_target';
                  evVal = 'widt';
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
                    [events_all.(sesName).(phaseName).data.phaseCount] == str2double(ns_evt{cols.(phaseName).phaseCount+1}(ec)) &...
                    ismember({events_all.(sesName).(phaseName).data.type},{stimType}) &...
                    [events_all.(sesName).(phaseName).data.trial] == str2double(ns_evt{cols.(phaseName).trial+1}(ec)) ...
                    );
                  
                  if length(this_event) > 1
                    warning('More than one event found! Fix this script before continuing analysis.')
                    keyboard
                  elseif isempty(this_event)
                    warning('No event found! Fix this script before continuing analysis.')
                    keyboard
                  elseif length(this_event) == 1
                    
                    phaseCount = this_event.phaseCount;
                    trial = this_event.trial;
                    response = this_event.resp;
                    acc = this_event.acc;
                    rt = this_event.rt;
                    
                    % add it to the trial definition
                    this_trl = trl_ini;
                    
                    % get the time of this event
                    this_sample = ft_event(i).sample;
                    
%                     % if we're using the photodiode DIN and we find one
%                     % within the threshold, replace the current sample time
%                     % with that of the DIN
%                     if cfg.eventinfo.usePhotodiodeDIN
%                       photodiodeDIN_thresholdSamp = round((cfg.eventinfo.photodiodeDIN_thresholdMS / 1000) * ft_hdr.Fs);
%                       if strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str)
%                         % if there is a DIN before and after the stim, pick
%                         % the closer one
%                         preDiff = (ft_event(i-1).sample - this_sample);
%                         postDiff = (ft_event(i+1).sample - this_sample);
%                         
%                         if preDiff < 0 && abs(preDiff) <= photodiodeDIN_thresholdSamp
%                           preFlag = true;
%                         else
%                           preFlag = false;
%                         end
%                         if postDiff <= photodiodeDIN_thresholdSamp
%                           postFlag = true;
%                         else
%                           postFlag = false;
%                         end
%                         
%                         if preFlag && ~postFlag
%                           % only the pre-DIN makes sense
%                           this_sample = ft_event(i-1).sample;
%                         elseif ~preFlag && postFlag
%                           % only the post-DIN makes sense
%                           this_sample = ft_event(i+1).sample;
%                         elseif preFlag && postFlag
%                           % choose the smaller one
%                           if abs(preDiff) < abs(postDiff)
%                             this_sample = ft_event(i-1).sample;
%                           elseif abs(preDiff) > abs(postDiff)
%                             this_sample = ft_event(i+1).sample;
%                           elseif abs(preDiff) == abs(postDiff)
%                             keyboard
%                           end
%                         end
%                       elseif strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && ~strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i+1).sample - this_sample) <= photodiodeDIN_thresholdSamp
%                         this_sample = ft_event(i+1).sample;
%                       elseif strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i-1).sample - this_sample) < 0 && abs(ft_event(i-1).sample - this_sample) <= photodiodeDIN_thresholdSamp
%                         % apparently the ethernet tags can be delayed
%                         % enough that the DIN shows up first
%                         this_sample = ft_event(i-1).sample;
%                       end
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
                    
                  end % check the event struct
                end % check the evt event
                
              else
                % the count is off? EEG and evt events don't line up
                keyboard
              end
              
              % put ec back in its prior state
              ec = ec + ec_add;
              
            case 'bgin'
              
            case 'trg+'
              
            case 'fix+'
              
            case 'resp'
             
          end
          
        end
      end % for
      fprintf('\n');
      
  end % switch phase
end % pha
% end % ses
