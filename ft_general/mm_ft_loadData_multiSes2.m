function [data,exper] = mm_ft_loadData_multiSes2(cfg,exper,dirs,ana,data_evoked)
%MM_FT_LOADDATA_MULTISES2 Load subject data into a full struct
%
% Load in power data, z-transform relative to all trials, average if wanted
%
% [data,exper] = mm_ft_loadData_multiSes(cfg,exper,ana,dirs,ana,data_evoked)
%
% Input:
%   cfg         = configuration; include normalization and baseline params
%   exper       = the exper struct
%   dirs        = the dirs struct, with the field saveDirProc (or saveDirRaw)
%   ana         = ana struct; contains a cell of the event values to load
%                 (ana.eventValues)
%
% cfg fields:
%   cfg.rmPreviousCfg=true; will remove data.cfg.previous with the intent
%                           to save memory
%
% Output:
%   data        = the data you requested
%   exper       = exper struct
%
% TODO:
%  - for cfg.baseline_type='condition', use ft_math
%  - many other things...
%


% no longer providing support for equatetrials
%   exper       = exper struct. If you used cfg.equatetrials='yes', a
%                 'randTrials' field gets added containing the randomly
%                 selected trial numbers for each subejct, session,
%                 condition type, and condition. You can feed this into
%                 this function again to choose the same trials (e.g., for
%                 getting power and phase/coherence of the same trials).


% OLD:
% ftype       = string included in the filename to load (e.g., 'tla' in
%               'data_tla_CR.mat'); can be 'raw' to load the raw data.
% keeptrials  = optional; default=1. If loading powspctrm data created with
%               keeptrials='yes', set to 0 to return averaged data.
%               NB: uses ft_freqdescriptives

% % make sure eventValues is set up correctly
% if isfield(exper,'eventValuesExtra') && isfield(exper.eventValuesExtra,'newValue') && ~isempty(exper.eventValuesExtra.newValue)
%   if isfield(exper.eventValuesExtra,'onlyKeepExtras') && exper.eventValuesExtra.onlyKeepExtras == 1
%     exper.eventValues = cat(2,exper.eventValuesExtra.newValue{:});
%   else
%     for nVal = 1:length(exper.eventValuesExtra.newValue)
%       if ~ismember(exper.eventValuesExtra.newValue{nVal},exper.eventValues)
%         exper.eventValues = cat(2,exper.eventValues,exper.eventValuesExtra.newValue{nVal});
%       else
%         fprintf('%s is already in the event value list!\n',exper.eventValuesExtra.newValue{nVal}{1});
%       end
%     end
%   end
%   exper.eventValues = sort(exper.eventValues);
% end

%% set some field names for FT data

cfg.fourierparam = 'fourierspctrm';
cfg.powparam = 'powspctrm';
%cohparam = 'cohspctrm';
% use 'powspctrm' for now because that's what FT expects
cfg.cohparam = 'powspctrm';
cfg.phaseparam = 'powspctrm';

%% check on congifuration details

if ~isfield(cfg,'loadMethod')
  error('Must specify cfg.loadMethod (''seg'' or ''trialinfo''');
end

% keep individual trials?
if ~isfield(cfg,'keeptrials')
  cfg.keeptrials = 'no';
  fprintf('Not keeping individual trials.\n');
elseif isfield(cfg,'keeptrials')
  if ~strcmp(cfg.keeptrials,'yes') && ~strcmp(cfg.keeptrials,'no')
    error('cfg.keeptrials must be set to ''yes'' or ''no''.');
  end
end

if ~isfield(cfg,'rmPreviousCfg')
  %warning('Upon loading of data, the data.cfg.previous field will NOT be removed! This may cause increased memory usage! See ''help %s'' for more info.',mfilename);
  %cfg.rmPreviousCfg = false;
  warning('Upon loading of data, the data.cfg.previous field WILL be removed to save memory! See ''help %s'' for more info.',mfilename);
  cfg.rmPreviousCfg = true;
end
if ~cfg.rmPreviousCfg
  warning('The data.cfg.previous field will NOT be removed! This may cause increased memory usage! See ''help %s'' for more info.',mfilename);
end

% % equate trial counts?
% if ~isfield(cfg,'equatetrials')
%   cfg.equatetrials = 'no';
%   fprintf('Not equating trials across conditions.\n');
% elseif isfield(cfg,'equatetrials')
%   if ~strcmp(cfg.equatetrials,'yes') && ~strcmp(cfg.equatetrials,'no')
%     error('cfg.equatetrials must be set to ''yes'' or ''no''.');
%   end
% end

% remove evoked power from whole power?
if ~isfield(cfg,'rmevoked')
  cfg.rmevoked = 'no';
  cfg.rmevokedfourier = 'no';
  cfg.rmevokedpow = 'no';
elseif isfield(cfg,'rmevoked')
  if ~strcmp(cfg.rmevoked,'yes') && ~strcmp(cfg.rmevoked,'no')
    error('cfg.rmevoked must be set to ''yes'' or ''no''.');
  elseif strcmp(cfg.rmevoked,'yes')
    % set the defaults
    if ~isfield(cfg,'rmevokedfourier')
      cfg.rmevokedfourier = 'yes';
    elseif isfield(cfg,'rmevokedfourier') && isempty(cfg.rmevokedfourier)
      cfg.rmevokedfourier = 'yes';
    end
    if ~isfield(cfg,'rmevokedpow')
      cfg.rmevokedpow = 'no';
    elseif isfield(cfg,'rmevokedpow') && isempty(cfg.rmevokedpow)
      cfg.rmevokedfourier = 'no';
    end
    
    if strcmp(cfg.rmevokedfourier,'yes') && strcmp(cfg.rmevokedpow,'yes')
      error('Cannot remove both evoked fourier and power');
    elseif strcmp(cfg.rmevokedfourier,'no') && strcmp(cfg.rmevokedpow,'no')
      % set the defaults
      warning([mfilename,':rmevoked method'],'Both remove options set to ''no'', setting defaults (cfg.rmevokedfourier = ''yes'').');
      cfg.rmevokedfourier = 'yes';
      cfg.rmevokedpow = 'no';
    end
  elseif strcmp(cfg.rmevoked,'no')
    cfg.rmevokedfourier = 'no';
    cfg.rmevokedpow = 'no';
  end
end

if ~exist('data_evoked','var') || isempty(data_evoked)
  data_evoked = [];
  if strcmp(cfg.rmevoked,'yes')
    warning([mfilename,':rmevoked'],'Setting cfg.rmevoked=''no'' because data_evoked was not provided.\n');
    cfg.rmevoked = 'no';
  end
end

if ~isfield(cfg,'zthresh')
  cfg.zthresh = [];
else
  if ~isnumeric(cfg.zthresh)
    error('cfg.zthresh not correctly set');
  end
end

%% check on input and output format

if ~isfield(cfg,'ftype')
  error('Must set cfg.ftype; used in the filename to load (e.g., ''fourier'' in ''data_fourier_CR.mat'');');
elseif strcmp(cfg.ftype,'pow')
  powflg = true;
  csdflg = false;
  fftflg = false;
  if ~isfield(cfg,'output')
    cfg.output = 'pow';
  end
elseif strcmp(cfg.ftype,'powandcsd')
  powflg = false;
  csdflg = true;
  fftflg = false;
  if ~isfield(cfg,'output')
    cfg.output = 'powandcsd';
  end
elseif strcmp(cfg.ftype,'fourier')
  powflg = false;
  csdflg = false;
  fftflg = true;
  if ~isfield(cfg,'output')
    cfg.output = 'fourier';
  elseif strcmp(cfg.output,'coh') && strcmp(cfg.keeptrials,'yes')
    warning([mfilename,':cohNoTrials'],'Setting cfg.keeptrials=''no'' because coherence is calcualted across trials.\n');
    cfg.keeptrials = 'no';
  elseif strcmp(cfg.output,'phase')
    if strcmp(cfg.keeptrials,'no')
      warning([mfilename,':cohNoTrials'],'Setting cfg.keeptrials=''yes'' because phase angle is calcualted for individual trials.\n');
      cfg.keeptrials = 'yes';
    end
    
    % TODO: plotphase is an unused field
    if ~isfield(cfg,'plotphase')
      cfg.plotphase = 'yes';
    end
    
    % define frequencies to use for phase calculation
    if ~isfield(cfg,'phasefreq')
      cfg.phasefreq = [];
    end
    
    % define channels to use for phase calculation
    if ~isfield(cfg,'phaseroi')
      cfg.phaseroi = {{'E55'}};
    elseif ~iscell(cfg.phaseroi)
      cfg.phaseroi = {{cfg.phaseroi}};
    elseif ~iscell(cfg.phaseroi{1})
      cfg.phaseroi = {cfg.phaseroi};
    end
  end
else
  % we haven't defined additional cases
  error('We haven''t defined additional cases.');
  %   powflg = false;
  %   csdflg = false;
  %   fftflg = false;
  %   if ~isfield(cfg,'output')
  %     cfg.output = '???';
  %   end
end

%% make sure some fields are set (correctly)

if ~isfield(cfg,'transform')
  cfg.transform = '';
elseif ~isempty(cfg.transform) && strcmp(cfg.output,'coh')
  error('Transformation for cfg.output=''coh'' is not allowed.');
end

if ~isfield(cfg,'baseline_type')
  cfg.baseline_type = '';
  warning('Default: Setting cfg.baseline_type to an empty string. Not baselining data!');
end

if ~isempty(cfg.baseline_type)
  if ~isfield(cfg,'baseline_events') || (isfield(cfg,'baseline_events') && isempty(cfg.baseline_events))
    error('No events specified for baseline in cfg.baseline_events!');
  end
end

% % no longer having 'db' be a reasonable transformation option
% if strcmp(cfg.transform,'db') && ~isempty(cfg.baseline_type) && ~strcmp(cfg.baseline_type,'relative')
%   warning([mfilename,':dbNormBL'],'db transformation uses a relative baseline type (division). Setting cfg.baseline_type=''relative''.\n');
%   cfg.baseline_type = 'relative';
% end

% if iscell(ana.eventValues)
%   if ~iscell(ana.eventValues{1})
%     ana.eventValues = {ana.eventValues};
%   end
% elseif ~iscell(ana.eventValues)
%   ana.eventValues = {{ana.eventValues}};
% end

% %% if we're equating trials, find out the number of trials to choose
%
% % TODO: see if randTrials is set up for multiple sessions
%
% if strcmp(cfg.equatetrials,'yes')
%   fprintf('Will be equating trials across conditions within a type.\n');
%   if isfield(exper,'randTrials')
%     fprintf('exper struct already contains randTrials field. Using that for sub-sampling trials.\n');
%   else
%     % initialize random number generator
%     rng('shuffle','twister');
%
%     for sub = 1:length(exper.subjects)
%       for ses = 1:length(exper.sesStr)
%         % initialize to find the lowest number of trials within each subject
%         lowEvNum = Inf(length(exper.subjects),length(ana.eventValues{ses}));
%         % initialize to store the randomly chosen events
%         exper.(exper.sesStr{ses}).randTrials = cell(length(exper.subjects),length(ana.eventValues{ses}));
%         for typ = 1:length(ana.eventValues{ses})
%           for evVal = 1:length(ana.eventValues{ses}{typ})
%             if exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub) < lowEvNum(sub,typ) && exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub) > 0
%               lowEvNum(sub,typ) = exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub);
%             end
%           end
%           exper.(exper.sesStr{ses}).randTrials{sub,typ} = nan(lowEvNum(sub,typ),length(ana.eventValues{ses}{typ}));
%           for evVal = 1:length(ana.eventValues{ses}{typ})
%             if exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub) ~= 0
%               exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal) = sort(randperm(exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub),lowEvNum(sub,typ)));
%             else
%               exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal) = [];
%             end
%           end % evVal
%         end % typ
%       end % ses
%     end % sub
%   end
%
%   % print out the trial counts for each subject
%   for typ = 1:length(ana.eventValues{ses})
%     fprintf('Condition type consisting of: %s\n',sprintf(repmat('%s ',1,length(ana.eventValues{ses}{typ})),ana.eventValues{ses}{typ}{:}));
%     if length(exper.subjects{1}) > 7
%       tabchar_sub = '\t';
%     else
%       tabchar_sub = '';
%     end
%     fprintf('Subject%s%s\n',sprintf(tabchar_sub),sprintf(repmat('\t%s',1,length(ana.eventValues{ses}{typ})),ana.eventValues{ses}{typ}{:}));
%     for ses = 1:length(exper.sesStr)
%       for sub = 1:length(exper.subjects)
%         subStr = exper.subjects{sub};
%         for evVal = 1:length(ana.eventValues{ses}{typ})
%           if length(ana.eventValues{ses}{typ}{evVal}) > 7
%             tabchar_ev = '\t';
%           else
%             tabchar_ev = '';
%           end
%           subStr = cat(2,subStr,sprintf('\t%d%s',exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub),sprintf(tabchar_ev)));
%         end
%         fprintf('%s\n',subStr);
%       end
%     end
%     fprintf('\n');
%   end
% end

%% process the data

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    %     % turn the session name into a string for easier printing
    %     if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
    %       exper.sesStr{ses} = exper.sessions{ses}{1};
    %       for i = 2:length(exper.sessions{ses})
    %         exper.sesStr{ses} = cat(2,exper.sesStr{ses},'_',exper.sessions{ses}{i});
    %       end
    %     elseif ~iscell(exper.sessions{ses}) || (iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1)
    %       exper.sesStr{ses} = exper.sessions{ses};
    %     end
    
    sesStr = exper.sesStr{ses};
    
    if strcmp(cfg.ftype,'raw')
      saveFileDir = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
    else
      saveFileDir = fullfile(dirs.saveDirProc,exper.subjects{sub},sesStr);
    end
    
    %for typ = 1:length(ana.eventValues{ses})
    for evVal = 1:length(ana.eventValues{ses})
      
      inputfile = fullfile(saveFileDir,sprintf('data_%s_%s.mat',cfg.ftype,ana.eventValues{ses}{evVal}));
      if exist(inputfile,'file')
        fprintf('Loading %s %s %s: %s...',exper.subjects{sub},sesStr,ana.eventValues{ses}{evVal},inputfile);
        % load the data
        subSesEvData = load(inputfile);
        fprintf('Done.\n');
        % get the name of the field
        fn = fieldnames(subSesEvData);
        
        % rename the field to 'data'
        if length(fn) == 1
          data_fn = cell2mat(fn);
        else
          error('More than one field in data struct! There should only be one.\n');
        end
        
        % save space by removing a potentially large 'previous' field
        if cfg.rmPreviousCfg
          if isfield(subSesEvData.(data_fn).cfg,'previous')
            fprintf('Removing ''previous'' field from data.cfg\n');
            subSesEvData.(data_fn).cfg = rmfield(subSesEvData.(data_fn).cfg,'previous');
          end
        end
        
        if strcmp(cfg.loadMethod,'seg')
          % if we're just loading the segmented files
          eventValue = ana.eventValues{ses}{evVal};
          
          if strcmp(cfg.rmevoked,'yes') && (strcmp(cfg.rmevokedfourier,'yes') || strcmp(cfg.rmevokedpow,'yes'))
            thisData_evoked = data_evoked.(sesStr).(eventValue).sub(sub).data;
          else
            thisData_evoked = [];
          end
          
          %data.(sesStr).(eventValue).sub(sub).data = processData(subSesEvData.(data_fn),exper.subjects{sub},sesStr,eventValue,cfg,dirs,ana,thisData_evoked,powflg,csdflg,fftflg);
          
        elseif strcmp(cfg.loadMethod,'trialinfo')
          % if we're using trialinfo loading method
          trl_order = ana.trl_order.(ana.eventValues{ses}{evVal});
          
          for es = 1:length(ana.eventValuesSplit{ses}{evVal})
            eventValue = ana.eventValuesSplit{ses}{evVal}{es};
            
            fprintf('\nSelecting %s trials...\n',eventValue);
            
            expr = ana.trl_expr{ses}{evVal}{es};
            if isfield(exper,'trialinfo_allEv')
              expr_allEv = ana.trl_expr{ses}{evVal}{es};
            end
            
            for to = 1:length(trl_order)
              % replace the field name in the logical expression
              % ana.trl_expr{ses}{es} with the column number found in
              % ana.trl_order.(ana.eventValues{ses}{evVal})
              
              % find the column number
              fieldNum = find(ismember(trl_order,trl_order{to}));
              
              % set up the string to be replaced
              r_exp = ['\<' trl_order{to} '\>'];
              
              % set up the replacement string
              r_str = sprintf('subSesEvData.%s.trialinfo(:,%d)',data_fn, fieldNum);
              expr = regexprep(expr,r_exp,r_str);
              
              if isfield(exper,'trialinfo_allEv')
                % set up the replacement string
                r_str_allEv = sprintf('exper.trialinfo_allEv.%s.%s{%d}(:,%d)',sesStr, ana.eventValues{ses}{evVal}, sub, fieldNum);
                expr_allEv = regexprep(expr_allEv,r_exp,r_str_allEv);
              end
            end
            
            % set up data selection
            cfg_sel = [];
            cfg_sel.trials = eval(expr);
            if sum(cfg_sel.trials) == 0
              error('no trials found for this setup: %s',expr);
            end
            if isfield(cfg,'latency')
              cfg_sel.latency = cfg.latency;
            else
              cfg_sel.latency = 'all';
            end
            if isfield(cfg,'frequency')
              if ischar(cfg.frequency)
                cfg_sel.frequency = cfg.frequency;
              elseif isnumeric(cfg.frequency) && length(cfg.frequency) == 2
                cfg_sel.foilim = cfg.frequency;
              end
            else
              cfg_sel.frequency = 'all';
            end
            %if isfield(cfg,'foilim')
            %  cfg_sel.foilim = cfg.foilim;
            %else
            %  cfg_sel.foilim = cfg.frequency;
            %end
            if isfield(cfg,'channel')
              cfg_sel.channel = cfg.channel;
            else
              cfg_sel.channel = 'all';
            end
            
            % select these data
            %selData = ft_selectdata(cfg_sel,subSesEvData.(data_fn));
            data.(sesStr).(eventValue).sub(sub).data = ft_selectdata(cfg_sel,subSesEvData.(data_fn));
            
            %             % select these data - use the old selection function
            %             if strcmp(cfg.ftype,'fourier')
            %               thisParam = cfg.fourierparam;
            %             elseif strcmp(cfg.ftype,'pow') || strcmp(cfg.ftype,'powandcsd')
            %               thisParam = cfg.powparam;
            %             end
            %             selData = ft_selectdata_old(...
            %               subSesEvData.(data_fn),'param',thisParam,'foilim',cfg_sel.foilim,'toilim',cfg_sel.latency,'rpt',cfg_sel.trials,'channel',cfg_sel.channel,...
            %               'avgoverchan','no','avgoverfreq','no','avgovertime','no','avgoverrpt','no');
            
            % put in the trial counts
            exper.nTrials.(sesStr).(eventValue)(sub,1) = sum(cfg_sel.trials);
            
            % put in the artifact counts
            if isfield(exper,'artifacts') && ~isempty(exper.artifacts.(sesStr).(ana.eventValues{ses}{evVal}))
              cfg_allEv = [];
              cfg_allEv.trials = eval(expr_allEv);
              theseArt = exper.artifacts.(sesStr).(ana.eventValues{ses}{evVal});
              artTypes = fieldnames(theseArt);
              totalArtCount = [];
              for at = 1:length(artTypes)
                if ~isempty(theseArt.(artTypes{at}){sub})
                  exper.artifacts.(sesStr).(eventValue).(artTypes{at}){sub,1} = theseArt.(artTypes{at}){sub}(cfg_allEv.trials);
                  totalArtCount = cat(2,totalArtCount,theseArt.(artTypes{at}){sub}(cfg_allEv.trials));
                  % % the below sum method does not take into account
                  % % that one trial can be marked for different kinds of
                  % % artifacts
                  % exper.nArtifacts.(sesStr).(eventValue).(artTypes{at})(sub,1) = sum(theseArt.(artTypes{at}){sub}(cfg_allEv.trials));
                else
                  exper.artifacts.(sesStr).(eventValue).(artTypes{at}){sub,1} = [];
                end
              end
              exper.nArtifacts.(sesStr).(eventValue)(sub,1) = sum(sign(sum(totalArtCount,2)));
            end
            
            %             if strcmp(cfg.rmevoked,'yes') && (strcmp(cfg.rmevokedfourier,'yes') || strcmp(cfg.rmevokedpow,'yes'))
            %               thisData_evoked = data_evoked.(sesStr).(eventValue).sub(sub).data;
            %             else
            %               thisData_evoked = [];
            %             end
            %
            %             data.(sesStr).(eventValue).sub(sub).data = processData(selData,exper.subjects{sub},sesStr,eventValue,cfg,dirs,ana,thisData_evoked,powflg,csdflg,fftflg);
          end % es
        end
        
        % % old method
        % data.(sesStr).(ana.eventValues{ses}{evVal}).sub(sub) = load(inputfile);
        fprintf('Done.\n\n');
      else
        warning([mfilename,':noFileFound'],'NOT FOUND: %s\n',inputfile);
        fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n\n');
        data.(sesStr).(ana.eventValues{ses}{evVal}).sub(sub).data.cfg.trl = [];
      end % if
    end % for evVal
    %end % for typ
    
    if ~isempty(cfg.baseline_type)
      % use mean and standard deviation across events to z-transform
      allEvBaseline = [];
      
      % get the baseline time indices
      blt = false(size(data.(sesStr).(ana.eventValuesSplit{ses}{1}{1}).sub(sub).data.time));
      blt_tbeg = nearest(data.(sesStr).(ana.eventValuesSplit{ses}{1}{1}).sub(sub).data.time,cfg.baseline_time(1));
      blt_tend = nearest(data.(sesStr).(ana.eventValuesSplit{ses}{1}{1}).sub(sub).data.time,cfg.baseline_time(2));
      blt(blt_tbeg:blt_tend) = true;
      
      % collect events for baseline
      if ischar(cfg.baseline_events)
        if strcmp(cfg.baseline_events,'all')
          % use all events listed in ana.eventValuesSplit
          for evVal = 1:length(ana.eventValuesSplit{ses})
            for es = 1:length(ana.eventValuesSplit{ses}{evVal})
              eventValue = ana.eventValuesSplit{ses}{evVal}{es};
              allEvBaseline = cat(1,allEvBaseline,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam)(:,:,:,blt));
            end
          end
        else
          % only one event, listed as a string
          eventTypes = fieldnames(data.(sesStr));
          if ismember(cfg.baseline_events,eventTypes)
            eventValue = cfg.baseline_events;
            allEvBaseline = cat(1,allEvBaseline,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam)(:,:,:,blt));
          else
            error('%s is not an event type stored in data struct!',cfg.baseline_events);
          end
        end
      elseif iscell(cfg.baseline_events)
        % events for baseline specified as a cell array
        for evVal = 1:length(cfg.baseline_events)
          eventValue = cfg.baseline_events{evVal};
          if isfield(data.(sesStr),eventValue)
            allEvBaseline = cat(1,allEvBaseline,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam)(:,:,:,blt));
          else
            error('%s is not an event type stored in data struct!',eventValue);
          end
        end
      end
      
      % mean across baseline period, then mean across events
      blm_orig = squeeze(nanmean(nanmean(allEvBaseline,4),1));
      
      if strcmp(cfg.baseline_type,'zscore')
        % std across time, then avg across events (lower freqs often get smaller std)
        blstd_orig = squeeze(nanmean(nanstd(allEvBaseline,0,4),1));
        % % avg across time, then std across events (higher freqs often get smaller std)
        % blstd = nanstd(nanmean(allEvBaseline,4),0,1);
        
        % % concatenate all times of all events and std (equivalent std across freqs)
        %blstd = shiftdim(squeeze(std(double(reshape(shiftdim(allEvBaseline,3),...
        %  size(allEvBaseline,1)*size(allEvBaseline,4),...
        %  size(allEvBaseline,2),size(allEvBaseline,3))),0,1)),-1);
      end
    end
    
    for evVal = 1:length(ana.eventValuesSplit{ses})
      for es = 1:length(ana.eventValuesSplit{ses}{evVal})
        eventValue = ana.eventValuesSplit{ses}{evVal}{es};
        
        if ~isempty(cfg.baseline_type)
          blm = shiftdim(repmat(blm_orig,[1,1,size(data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),1)]),2);
          
          if strcmp(cfg.baseline_type,'zscore')
            fprintf('\tZ-transforming data relative to mean([%.2f %.2f]).\n',cfg.baseline_time(1),cfg.baseline_time(2));
            blstd = shiftdim(repmat(blstd_orig,[1,1,size(data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),1)]),2);
            
            data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam) = bsxfun(@rdivide,bsxfun(@minus,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),blm),blstd);
            
          elseif strcmp(cfg.baseline_type,'db')
            fprintf('\tConverting to db relative to mean([%.2f %.2f]).\n',cfg.baseline_time(1),cfg.baseline_time(2));
            % divide by baseline mean and convert to db
            
            %if strcmp(cfg.norm_trials,'average')
            data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam) = 10*log10(bsxfun(@rdivide,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),blm));
            %elseif strcmp(cfg.norm_trials,'single')
            %  data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam) = bsxfun(@rdivide,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),blm);
            %end
            
            % TODO: mean(10*log10) is not the same as 10*log10(mean)
            
            %data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam) = bsxfun(@rdivide,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),blm);
            
          elseif strcmp(cfg.baseline_type,'absolute')
            fprintf('\tSubtracting mean([%.2f %.2f]) power from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
            data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam) = bsxfun(@minus,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),blm);
            
          elseif strcmp(cfg.baseline_type,'relative')
            fprintf('\tDividing entire trial power by mean([%.2f %.2f]) power.\n',cfg.baseline_time(1),cfg.baseline_time(2));
            data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam) = bsxfun(@rdivide,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),blm);
            
          elseif strcmp(cfg.baseline_type,'relchange')
            fprintf('\tSubtracting mean([%.2f %.2f]) power from entire trial and dividing by that mean.\n',cfg.baseline_time(1),cfg.baseline_time(2));
            data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam) = bsxfun(@rdivide,bsxfun(@minus,data.(sesStr).(eventValue).sub(sub).data.(cfg.powparam),blm),blm);
          end
        end
        
        % average if needed
        if strcmp(cfg.keeptrials,'no')
          cfg_sel = [];
          cfg_sel.avgoverrpt = 'yes';
          data.(sesStr).(eventValue).sub(sub).data = ft_selectdata(cfg_sel,data.(sesStr).(eventValue).sub(sub).data);
        end
        
      end
    end
    
  end % for ses
end % for sub

end % function mm_ft_loadData_multiSes

% function [subSesEvData] = processData(subSesEvData,subject,sesStr,eventValue,cfg,dirs,ana,data_evoked,powflg,csdflg,fftflg)
% 
% % can't have values exactly equal to zero (why?)
% if powflg
%   param = 'powspctrm';
%   subSesEvData.(param)(subSesEvData.(param) == 0) = eps(0);
% end
% 
% if csdflg
%   param = 'crsspctrm';
% end
% 
% % convert complex fourier-spectra to power or coherence
% if fftflg
%   
%   %               % TODO: do I need to update something in the FT structure now
%   %               % that I get the subset of trials at the start? maybe use
%   %               % ft_selectdata?
%   %               %
%   %               % equate the trial counts
%   %               if strcmp(cfg.equatetrials,'yes')
%   %                 % get the subset of trials
%   %                 if size(subSesEvData.(cfg.fourierparam),1) ~= length(exper.(sesStr).randTrials{sub,typ}(:,evVal))
%   %                   fprintf('Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.(sesStr).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(cfg.fourierparam),1));
%   %                   %subSesEvData.(cfg.fourierparam) = subSesEvData.(cfg.fourierparam)(exper.(sesStr).randTrials{sub,typ}(:,evVal),:,:,:);
%   %                   subSesEvData = ft_selectdata_old(subSesEvData,'param',cfg.fourierparam,'rpt',exper.(sesStr).randTrials{sub,typ}(:,evVal));
%   %                 end
%   %               end
%   
%   fprintf('\tConverting %s to %s.\n',cfg.ftype,cfg.output);
%   if strcmp(cfg.output,'pow')
%     %% get the modulus of the whole data and evoked data
%     param = 'modulus';
%     
%     % get the modulus for whole data (i.e., abs(whole fourier))
%     subSesEvData.(param) = abs(subSesEvData.(cfg.fourierparam));
%     
%     % get the modulus for evoked data (i.e., abs(evoked fourier))
%     if strcmp(cfg.rmevoked,'yes')
%       data_evoked.(param) = abs(data_evoked.(cfg.fourierparam));
%     end
%     
%     if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
%       % option 1: subtract evoked modulus from whole modulus;
%       % this will make evoked=positive and induced=negative
%       
%       %fprintf('Subtraction: (abs(whole fourier) - abs(evoked fourier)).^2 = induced power.\n');
%       fprintf('\tSubtraction: ((whole modulus) - (evoked modulus)).^2 = induced power.\n');
%       
%       % subtract evoked from single trials, then square
%       % % evoked_data = repmat(data_evoked.(cfg.fourierparam),[size(subSesEvData.(cfg.fourierparam),1),1,1,1]);
%       % % subSesEvData.(cfg.powparam) = (abs(subSesEvData.(cfg.fourierparam)) - abs(evoked_data)).^2;
%       evoked_data = repmat(data_evoked.(param),[size(subSesEvData.(param),1),1,1,1]);
%       subSesEvData.(cfg.powparam) = (subSesEvData.(param) - evoked_data).^2;
%       
%       % % subtract means of moduluses, then square
%       % % subSesEvData.(cfg.powparam) = (mean(abs(subSesEvData.(cfg.fourierparam)),1) - abs(data_evoked.(cfg.fourierparam))).^2;
%       % subSesEvData.(cfg.powparam) = (mean(subSesEvData.(param),1) - data_evoked.(param)).^2;
%       
%       % % square means of moduluses, then subtract
%       % % subSesEvData.(cfg.powparam) = mean(abs(subSesEvData.(cfg.fourierparam)),1).^2 - abs(data_evoked.(cfg.fourierparam)).^2;
%       % subSesEvData.(cfg.powparam) = mean(subSesEvData.(param),1).^2 - data_evoked.(param).^2;
%       
%     elseif strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
%       % option 2: set up evoked fourier to subtract as power
%       if isfield(data_evoked,cfg.powparam)
%         data_evoked = rmfield(data_evoked,cfg.powparam);
%       end
%       
%       % get the evoked power from evoked fourier data
%       data_evoked.(cfg.powparam) = data_evoked.(cfg.fourierparam).^2;
%       % no zeros
%       data_evoked.(cfg.powparam)(data_evoked.(cfg.powparam) == 0) = eps(0);
%       % no longer need the fourierspctrm or modulus fields
%       data_evoked = rmfield(data_evoked,{cfg.fourierparam param});
%       
%       % square the whole modulus to get power
%       subSesEvData.(cfg.powparam) = subSesEvData.(param).^2;
%     elseif strcmp(cfg.rmevoked,'no')
%       %% otherwise just square the modulus to get power
%       
%       subSesEvData.(cfg.powparam) = subSesEvData.(param).^2;
%     end
%     
%     % no zeros
%     subSesEvData.(cfg.powparam)(subSesEvData.(cfg.powparam) == 0) = eps(0);
%     % no longer need the fourierspctrm or modulus fields
%     subSesEvData = rmfield(subSesEvData,{cfg.fourierparam param});
%     
%     param = cfg.powparam;
%     
%   elseif strcmp(cfg.output,'coh')
%     %% inter-trial phase coherence
%     
%     %if strcmp(cfg.equatetrials,'yes')
%     %  % get the subset of trials
%     %  fprintf('Inter-trial coherence data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.(sesStr).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(cfg.fourierparam),1));
%     %  subSesEvData.(cfg.fourierparam) = subSesEvData.(cfg.fourierparam)(exper.(sesStr).randTrials{sub,typ}(:,evVal),:,:,:);
%     %end
%     
%     subSesEvData.(cfg.cohparam) = subSesEvData.(cfg.fourierparam) ./ abs(subSesEvData.(cfg.fourierparam));
%     subSesEvData.(cfg.cohparam) = abs(squeeze(nanmean(subSesEvData.(cfg.cohparam),1)));
%     subSesEvData = rmfield(subSesEvData,cfg.fourierparam);
%     param = cfg.cohparam;
%     % change dimord
%     % TODO: can I automate this (e.g., with ft_checkdata)? or
%     % maybe ft_freqdescriptives (keeptrials='no')
%     subSesEvData.dimord = 'chan_freq_time';
%     subSesEvData = rmfield(subSesEvData,'trialinfo');
%   elseif strcmp(cfg.output,'phase')
%     %% single-trial phase
%     
%     %if strcmp(cfg.equatetrials,'yes')
%     %  % get the subset of trials
%     %  fprintf('Phase data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.(sesStr).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(cfg.fourierparam),1));
%     %  subSesEvData.(cfg.fourierparam) = subSesEvData.(cfg.fourierparam)(exper.(sesStr).randTrials{sub,typ}(:,evVal),:,:,:);
%     %end
%     
%     % use all the frequencies
%     if isempty(cfg.phasefreq)
%       cfg.phasefreq = subSesEvData.freq;
%     end
%     
%     if ~isvector(cfg.phasefreq) || (size(cfg.phasefreq,1) == 1 && size(cfg.phasefreq,2) == 2)
%       meanfourier = nan(size(subSesEvData.(cfg.fourierparam),1),length(cfg.phaseroi),size(cfg.phasefreq,1),size(subSesEvData.(cfg.fourierparam),4));
%     elseif iscolumn(cfg.phasefreq)
%       meanfourier = nan(size(subSesEvData.(cfg.fourierparam),1),length(cfg.phaseroi),length(cfg.phasefreq),size(subSesEvData.(cfg.fourierparam),4));
%     end
%     
%     for r = 1:length(cfg.phaseroi)
%       cfg.roi = cfg.phaseroi{r};
%       % set the channel information
%       if ismember(cfg.roi,ana.elecGroupsStr)
%         % if it's in the predefined ROIs, get the channel numbers
%         cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
%       else
%         % otherwise it should be the channel number(s) or 'all'
%         cfg.channel = cfg.roi;
%       end
%       chansel = ismember(subSesEvData.label,cfg.channel);
%       
%       if ~isvector(cfg.phasefreq) || (size(cfg.phasefreq,1) == 1 && size(cfg.phasefreq,2) == 2)
%         for f = 1:size(cfg.phasefreq,1)
%           freqsel = subSesEvData.freq >= cfg.phasefreq(f,1) & subSesEvData.freq <= cfg.phasefreq(f,2);
%           meanfourier(:,r,f,:) = nanmean(nanmean(subSesEvData.(cfg.fourierparam)(:,chansel,freqsel,:),3),2);
%         end
%       elseif iscolumn(cfg.phasefreq)
%         meanfourier(:,r,:,:) = nanmean(subSesEvData.(cfg.fourierparam)(:,chansel,:,:),2);
%       end
%     end
%     subSesEvData.(cfg.fourierparam) = meanfourier;
%     
%     % calculate phase angle
%     subSesEvData.(cfg.phaseparam) = angle(subSesEvData.(cfg.fourierparam));
%     subSesEvData = rmfield(subSesEvData,cfg.fourierparam);
%     
%     % print it
%     for r = 1:length(cfg.phaseroi)
%       cfg.roi = cfg.phaseroi{r};
%       % set the channel information
%       if ismember(cfg.roi,ana.elecGroupsStr)
%         % if it's in the predefined ROIs, get the channel numbers
%         cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
%         % set the string for the filename
%         chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
%       else
%         % otherwise it should be the channel number(s) or 'all'
%         cfg.channel = cfg.roi;
%         
%         % set the string for the filename
%         if isfield(cfg,'cohrefchannel')
%           chan_str = [cfg.cohrefchannel,'-',sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:})];
%         else
%           chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
%         end
%       end
%       
%       for f = 1:size(subSesEvData.(cfg.phaseparam),3)
%         if ~isvector(cfg.phasefreq) || length(cfg.phasefreq) == 2
%           freq_str = sprintf('%.1f-%.1f',cfg.phasefreq(f,1),cfg.phasefreq(f,2));
%         else
%           freq_str = sprintf('%.1f',cfg.phasefreq(f));
%         end
%         figure
%         surf(subSesEvData.time,1:size(subSesEvData.(cfg.phaseparam),1),squeeze(subSesEvData.(cfg.phaseparam)(:,r,f,:)));
%         shading interp;view([0,90]);axis tight;
%         title(sprintf('%s, %s, %s, %s, %s Hz',subject,strrep(sesStr,'_',''),eventValue,strrep(chan_str,'_',' '),freq_str));
%         xlabel('Time (s)');
%         ylabel('Trial number');
%         publishfig(gcf,0);
%         figfilename = sprintf('phase_%s_%s_%s_%s%sHz.png',subject,sesStr,eventValue,chan_str,freq_str);
%         dirs.saveDirFigsPhase = fullfile(dirs.saveDirFigs,'tfr_phase');
%         if ~exist(dirs.saveDirFigsPhase,'dir')
%           mkdir(dirs.saveDirFigsPhase)
%         end
%         print(gcf,'-dpng','-r150',fullfile(dirs.saveDirFigsPhase,figfilename));
%         close all
%       end % f
%     end % r
%     
%   end % cfg.output
% end % fftflg
% 
% %% transform the data
% if ~isempty(cfg.transform)
%   fprintf('\tDoing a %s transformation.\n',cfg.transform);
%   
%   if strcmp(cfg.transform,'vec')
%     for tr = 1:size(subSesEvData.(param),1)
%       for ch = 1:size(subSesEvData.(param),2)
%         % doesn't ignore NaNs, but is way faster
%         subSesEvData.(param)(tr,ch,:,:) = normr(squeeze(subSesEvData.(param)(tr,ch,:,:)));
%         % % try to ignore NaNs
%         %for fq = 1:size(subSesEvData.(param),3)
%         %  thisTrl = squeeze(subSesEvData.(param)(tr,ch,fq,:));
%         %  subSesEvData.(param)(tr,ch,fq,~isnan(thisTrl)) = normc(thisTrl(~isnan(thisTrl)));
%         %end
%         % % trying another method, but it doesn't work
%         %thisTrl = squeeze(subSesEvData.(param)(tr,ch,:,:));
%         %subSesEvData.(param)(tr,ch,:,~isnan(thisTrl)) = normr(thisTrl(~isnan(thisTrl)));
%       end
%     end
%   else
%     
%     if ~strcmp(cfg.transform,'log10') && ~strcmp(cfg.transform,'log')
%       % just letting you know...
%       warning([mfilename,':unsure_about_cfgTransform'],'Unsure if %s is a reasonable trasnformation.\n',evcfg.transformentValue);
%     end
%     
%     % run using feval (e.g., for log10 normalization)
%     subSesEvData.(param) = feval(cfg.transform,subSesEvData.(param));
%     
%     % option 2: normalize the evoked data
%     if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
%       data_evoked.(cfg.powparam) = feval(cfg.transform,data_evoked.(cfg.powparam));
%     end
%   end
% end
% 
% % option 2: subtract evoked power from whole power; need to
% % do this after log normalization
% if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
%   fprintf('\tSubtracting evoked power from whole power to get induced power.\n');
%   evoked_data = repmat(data_evoked.(cfg.powparam),[size(subSesEvData.(cfg.powparam),1),1,1,1]);
%   subSesEvData.(cfg.powparam) = subSesEvData.(cfg.powparam) - evoked_data;
% end
% 
% %% baseline correct power
% 
% % single trial normalization uses the entire trial for
% % normalization
% if strcmp(cfg.baseline_data,'pow')
%   if ~isempty(cfg.baseline_type) && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow'))
%     
%     if ndims(subSesEvData.(param)) ~= 4
%       error('data (%s) does not have four dimensions',param);
%     end
%     
%     if strcmp(cfg.norm_trials,'single')
%       % Grandchamp & Delorme (2011)
%       
%       % single-trial normalization uses all the time points
%       blt = true(size(subSesEvData.time));
%       
%       % mean across baseline period
%       blm = nanmean(subSesEvData.(param)(:,:,:,blt),4);
%       
%       if strcmp(cfg.baseline_type,'zscore')
%         fprintf('\tZ-transforming each trial using entire trial as baseline.\n');
%         % std across time
%         blstd = nanstd(subSesEvData.(param)(:,:,:,blt),0,4);
%         
%         % z-transform
%         subSesEvData.(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(param),blm),blstd);
%         
%       elseif strcmp(cfg.baseline_type,'db')
%         fprintf('\tConverting each trial to db using entire trial as baseline.\n');
%         % divide by baseline mean and convert to db
%         subSesEvData.(param) = 10*log10(bsxfun(@rdivide,subSesEvData.(param),blm));
%         
%         % mean(10*log10) is not the same as 10*log10(mean)
%         
%         % debug
%         %figure;imagesc(squeeze(mean(subSesEvData.(param)(:,62,:,:),1)));axis xy;colorbar
%         
%       elseif strcmp(cfg.baseline_type,'absolute')
%         fprintf('\tSubtracting mean(entire trial) power from entire trial.\n');
%         subSesEvData.(param) = bsxfun(@minus,subSesEvData.(param),blm);
%         
%       elseif strcmp(cfg.baseline_type,'relative')
%         fprintf('\tDividing entire trial power by mean(entire trial).\n');
%         subSesEvData.(param) = bsxfun(@rdivide,subSesEvData.(param),blm);
%         
%       elseif strcmp(cfg.baseline_type,'relchange')
%         fprintf('\tSubtracting mean(entire trial) power from entire trial and dividing by that mean.\n');
%         subSesEvData.(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(param),blm),blm);
%       end
%       
%     end
%     
%     % now compare it to the baseline activity
%     if strcmp(cfg.norm_trials,'average') || strcmp(cfg.norm_trials,'single')
%       % TODO: can I combine 'single' with comparison to baseline?
%       
%       % get the baseline time indices
%       blt = subSesEvData.time >= cfg.baseline_time(1) & subSesEvData.time <= cfg.baseline_time(2);
%       
%       % mean across baseline period
%       blm = nanmean(subSesEvData.(param)(:,:,:,blt),4);
%       
%       if strcmp(cfg.baseline_type,'zscore')
%         fprintf('\tZ-transforming data relative to mean([%.2f %.2f]).\n',cfg.baseline_time(1),cfg.baseline_time(2));
%         % std across time, then avg across events (lower freqs often get smaller std)
%         blstd = nanmean(nanstd(subSesEvData.(param)(:,:,:,blt),0,4),1);
%         
%         % % avg across time, then std across events (higher freqs often get smaller std)
%         % blstd = nanstd(nanmean(subSesEvData.(param)(:,:,:,blt),4),0,1);
%         
%         % % concatenate all times of all events and std (equivalent std across freqs)
%         %blstd = shiftdim(squeeze(std(double(reshape(shiftdim(subSesEvData.(param)(:,:,:,blt),3),...
%         %  size(subSesEvData.(param)(:,:,:,blt),1)*size(subSesEvData.(param)(:,:,:,blt),4),...
%         %  size(subSesEvData.(param)(:,:,:,blt),2),size(subSesEvData.(param)(:,:,:,blt),3))),0,1)),-1);
%         
%         subSesEvData.(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(param),blm),blstd);
%         
%       elseif strcmp(cfg.baseline_type,'db')
%         fprintf('\tConverting to db relative to mean([%.2f %.2f]).\n',cfg.baseline_time(1),cfg.baseline_time(2));
%         % divide by baseline mean and convert to db
%         
%         if strcmp(cfg.norm_trials,'average')
%           subSesEvData.(param) = 10*log10(bsxfun(@rdivide,subSesEvData.(param),blm));
%         elseif strcmp(cfg.norm_trials,'single')
%           subSesEvData.(param) = bsxfun(@rdivide,subSesEvData.(param),blm);
%         end
%         
%         % TODO: mean(10*log10) is not the same as 10*log10(mean)
%         
%         %subSesEvData.(param) = bsxfun(@rdivide,subSesEvData.(param),blm);
%         
%       elseif strcmp(cfg.baseline_type,'absolute')
%         fprintf('\tSubtracting mean([%.2f %.2f]) power from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
%         subSesEvData.(param) = bsxfun(@minus,subSesEvData.(param),blm);
%         
%       elseif strcmp(cfg.baseline_type,'relative')
%         fprintf('\tDividing entire trial power by mean([%.2f %.2f]) power.\n',cfg.baseline_time(1),cfg.baseline_time(2));
%         subSesEvData.(param) = bsxfun(@rdivide,subSesEvData.(param),blm);
%         
%       elseif strcmp(cfg.baseline_type,'relchange')
%         fprintf('\tSubtracting mean([%.2f %.2f]) power from entire trial and dividing by that mean.\n',cfg.baseline_time(1),cfg.baseline_time(2));
%         subSesEvData.(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(param),blm),blm);
%       end
%       
%     end % norm_trials
%   end % baseline_type
% end % baseline_data
% 
% %             if strcmp(cfg.baseline_data,'pow')
% %               if ~isempty(cfg.baseline_type) && ~strcmp(cfg.transform,'db')
% %
% %                 % get the baseline time indices
% %                 blt = subSesEvData.time >= cfg.baseline_time(1) & subSesEvData.time <= cfg.baseline_time(2);
% %
% %                 if strcmp(cfg.baseline_type,'zscore') && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(param)) == 4
% %                   fprintf('Z-transforming data relative to mean([%.2f %.2f]).\n',cfg.baseline_time(1),cfg.baseline_time(2));
% %
% %                   nTrl = size(subSesEvData.(param),1);
% %                   %nFrq = size(subSesEvData.(param),3);
% %                   nSmp = size(subSesEvData.(param),4);
% %
% %                   blm = repmat(nanmean(subSesEvData.(param)(:,:,:,blt),4),[1,1,1,nSmp]);
% %                   blstd = repmat(nanmean(nanstd(subSesEvData.(param)(:,:,:,blt),0,1),4),[nTrl,1,1,nSmp]);
% %                   subSesEvData.(param) = (subSesEvData.(param) - blm) ./ blstd;
% %
% %                 elseif (strcmp(cfg.baseline_type,'absolute') || strcmp(cfg.baseline_type,'relative') || strcmp(cfg.baseline_type,'relchange')) && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(param)) == 4
% %                   %blt = subSesEvData.time >= cfg.baseline_time(1) & subSesEvData.time <= cfg.baseline_time(2);
% %
% %                   nSmp = size(subSesEvData.(param),4);
% %
% %                   blm = repmat(nanmean(subSesEvData.(param)(:,:,:,blt),4),[1,1,1,nSmp]);
% %
% %                   if strcmp(cfg.baseline_type,'absolute')
% %                     fprintf('Subtracting mean([%.2f %.2f]) power from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
% %                     subSesEvData.(param) = subSesEvData.(param) - blm;
% %                   elseif strcmp(cfg.baseline_type,'relative')
% %                     error('relative baseline is not working');
% %                     %fprintf('Dividing power by mean([%.2f %.2f]).\n',cfg.baseline_time(1),cfg.baseline_time(2));
% %                     %subSesEvData.(param) = subSesEvData.(param) ./ blm;
% %                   elseif strcmp(cfg.baseline_type,'relchange')
% %                     error('relchange baseline is not working');
% %                     %fprintf('Subtracting mean([%.2f %.2f]) power from entire trial and dividing by the mean.\n',cfg.baseline_time(1),cfg.baseline_time(2));
% %                     %subSesEvData.(param) = (subSesEvData.(param) - blm) ./ blm;
% %                   end
% %
% %                 elseif strcmp(cfg.baseline_type,'absolute') && strcmp(cfg.output,'coh') && ndims(subSesEvData.(param)) == 3
% %                   % not sure that coherence gets baseline correction???
% %
% %                   fprintf('Subtracting mean([%.2f %.2f]) coherence from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
% %                   %blt = subSesEvData.time >= cfg.baseline_time(1) & subSesEvData.time <= cfg.baseline_time(2);
% %
% %                   nSmp = size(subSesEvData.(param),3);
% %
% %                   blm = repmat(nanmean(subSesEvData.(param)(:,:,blt),3),[1,1,nSmp]);
% %                   subSesEvData.(param) = (subSesEvData.(param) - blm);
% %                 end
% %               end
% %             end
% 
% %% debug
% %             chan = 62;
% %             chan = 25;
% %             clim = [-1 1];
% %             %clim = [-3 3];
% %             if strcmp(cfg.transform,'db')
% %               clim = clim .* 10;
% %             end
% %             figure;surf(subSesEvData.time,subSesEvData.freq,squeeze(mean(subSesEvData.(cfg.powparam)(:,chan,:,:),1)));shading interp;view([0,90]);axis tight;colorbar;
% %             title(sprintf('chan %d: rmevokedfourier: %s; rmevokedpow: %s',chan,cfg.rmevokedfourier,cfg.rmevokedpow));
% %             caxis(clim)
% %             figure;imagesc(subSesEvData.time,subSesEvData.freq,squeeze(mean(subSesEvData.(cfg.powparam)(:,chan,:,:),1)),clim);axis xy;colorbar
% %             title(sprintf('chan %d: rmevokedfourier: %s; rmevokedpow: %s',chan,cfg.rmevokedfourier,cfg.rmevokedpow));
% %
% %             figure
% %             for tr = 1:size(subSesEvData.(cfg.powparam),1)
% %               clf
% %               surf(subSesEvData.time,subSesEvData.freq,squeeze(subSesEvData.(cfg.powparam)(tr,chan,:,:)));shading interp;view([0,90]);axis tight;colorbar;
% %               title(sprintf('chan %d: trial %d/%d',chan,tr,size(subSesEvData.(cfg.powparam),1)));
% %               caxis([-3 3])
% %               keyboard
% %             end
% 
% %% save the data in a container struct
% 
% fprintf('\t%s %s %s: ',subject,sesStr,eventValue);
% if (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && isfield(subSesEvData,cfg.powparam) && ndims(subSesEvData.(cfg.powparam)) == 4
%   fprintf('Power data: ');
%   %             elseif strcmp(cfg.keeptrials,'no') && (~strcmp(cfg.ftype,'pow') || ~strcmp(cfg.output,'pow')) && isfield(subSesEvData,'powspctrm') && ndims(subSesEvData.(cfg.powparam)) == 4
%   %               error('\n%s %s %s: Can only keep trials for cfg.ftype=''pow''. You set it to ''%s''.\n',subject,sesStr,eventValue,cfg.ftype);
%   %             elseif strcmp(cfg.keeptrials,'no') && ~isfield(subSesEvData,'powspctrm')
%   %               error('\n%s %s %s: Can only keep trials with ''powspctrm'' field. Please examine your data.\n',subject,sesStr,eventValue);
%   %             elseif strcmp(cfg.keeptrials,'no') && isfield(subSesEvData,'powspctrm') && ndims(subSesEvData.(cfg.powparam)) ~= 4
%   %               error('\n%s %s %s: Can only keep trials for ndims(powspctrm)==4. This data has ndims=%d.\n',subject,sesStr,eventValue,ndims(subSesEvData.(cfg.powparam)));
%   
%   % set the configuration of ft_freqdescriptives
%   cfg_fd = [];
%   cfg_fd.keeptrials = cfg.keeptrials;
%   
%   if strcmp(cfg.keeptrials,'no')
%     fprintf('Averaging over individual trials.\n');
%   elseif strcmp(cfg.keeptrials,'yes')
%     fprintf('Keeping individual trials.\n');
%   end
%   
%   %               if strcmp(cfg.equatetrials,'yes')
%   %                 fprintf(' Equating trial counts across all conditions (within a type).\n');
%   %                 fprintf('Randomly selecting %d of %d trials.\n',length(exper.(sesStr).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(cfg.powparam),1));
%   %                 % get the subset of trials
%   %                 cfg_fd.trials = exper.(sesStr).randTrials{sub,typ}(:,evVal);
%   %               elseif strcmp(cfg.equatetrials,'no')
%   %                 fprintf(' Using all trials.\n');
%   %                 % use all of the trials
%   %                 cfg_fd.trials = 'all';
%   %               end
%   
%   % Do a very coarse rejection for trials with artifacts
%   if ~isempty(cfg.zthresh) && strcmp(cfg.baseline_type,'zscore')
%     % reject trials with huge zscores
%     cfg_fd.trials = true(1,size(subSesEvData.(param),1));
%     [tr,ch] = find(squeeze(nanmean(nanmean(subSesEvData.(param),4),3)) > cfg.zthresh);
%     if ~isempty(tr)
%       fprintf('\tRejecting %d of %d trials due to zscore > %.2f in %d channels (NB: this is a very coarse rejection).\n',length(unique(tr)),length(cfg_fd.trials),cfg.zthresh,length(unique(ch)));
%       cfg_fd.trials(unique(tr)) = false;
%       fprintf('\tKeeping %d of %d trials.\n',sum(ismember(cfg_fd.trials,1)),length(cfg_fd.trials));
%     else
%       fprintf('\tKeeping all %d trials.\n',size(subSesEvData.(param),1));
%       cfg_fd.trials = 'all';
%     end
%   else
%     % TODO: do I need to update something in the FT structure now
%     % that I get the subset of trials at the start?
%     cfg_fd.trials = 'all';
%   end
%   
%   % make sure there are trials left to average
%   if (isvector(cfg_fd.trials) && ismember(1,cfg_fd.trials)) || (ischar(cfg_fd.trials) && strcmp(cfg_fd.trials,'all'))
%     % use ft_freqdescriptives to average over individual
%     % trials, if desired and the data is appropriate
%     % subSesEvData = ft_freqdescriptives(cfg_fd,subSesEvData);
%     % % data.(sesStr).(eventValue).sub(sub).data = ft_freqdescriptives(cfg_fd,subSesEvData);
%     
%     if strcmp(cfg_fd.keeptrials,'yes')
%       cfg_fd.avgoverrpt = 'no';
%     elseif strcmp(cfg_fd.keeptrials,'no')
%       cfg_fd.avgoverrpt = 'yes';
%     end
%     if isfield(cfg_fd,'keeptrials')
%       cfg_fd = rmfield(cfg_fd,'keeptrials');
%     end
%     subSesEvData = ft_selectdata(cfg_fd,subSesEvData);
%   else %if (isvector(cfg_fd.trials) && length(unique(cfg_fd.trials)) == 1 && unique(cfg_fd.trials) == 0)
%     warning([mfilename,':allTrialsRejected'],'ALL TRIALS REJECTED: %s\n',eventValue);
%     fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n\n');
%     subSesEvData.cfg.trl = [];
%     %data.(sesStr).(eventValue).sub(sub).data.cfg.trl = [];
%   end % if
%   
%   
% elseif strcmp(cfg.output,'coh') && isfield(subSesEvData,cfg.cohparam) && ndims(subSesEvData.(cfg.cohparam)) == 3
%   fprintf('\tPhase coherence data: Individual trials are lost.\n');
%   %               if strcmp(cfg.equatetrials,'no')
%   %                 fprintf(' Using all trials. NB: this could affect the coherence calculation if trial counts are unequal!!\n');
%   %               %elseif strcmp(cfg.equatetrials,'yes')
%   %               %  fprintf(' Using randomly selected subsets of trials so each condition (within a type) has the same number of trials.\n');
%   %               end
%   
%   %data.(sesStr).(eventValue).sub(sub).data = subSesEvData;
% elseif strcmp(cfg.output,'phase') && isfield(subSesEvData,cfg.phaseparam) && ndims(subSesEvData.(cfg.phaseparam)) == 4
%   fprintf('\tOutputting phase data.\n');
%   %data.(sesStr).(eventValue).sub(sub).data = subSesEvData;
% else
%   fprintf('\tNeed to figure out what to do for this case.\n');
%   keyboard
% end
% 
% end % function processData
