function [data,exper] = mm_ft_loadData_multiSes(cfg,exper,dirs,ana,data_evoked)
%MM_FT_LOADDATA_MULTISES Load subject data into a full struct
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
% Output:
%   data        = the data you requested
%   exper       = exper struct. If you used cfg.equatetrials='yes', a
%                 'randTrials' field gets added containing the randomly
%                 selected trial numbers for each subejct, session,
%                 condition type, and condition. You can feed this into
%                 this function again to choose the same trials (e.g., for
%                 getting power and phase/coherence of the same trials).
%
% TODO:
%  - for cfg.baseline_type='condition', use ft_freqcomparison
%  - many other things...
%



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

fourierparam = 'fourierspctrm';
powparam = 'powspctrm';
%cohparam = 'cohspctrm';
% use 'powspctrm' for now because that's what FT expects
cohparam = 'powspctrm';
phaseparam = 'powspctrm';

%% check on congifuration details

% keep individual trials?
if ~isfield(cfg,'keeptrials')
  cfg.keeptrials = 'no';
  fprintf('Not keeping individual trials.\n');
elseif isfield(cfg,'keeptrials')
  if ~strcmp(cfg.keeptrials,'yes') && ~strcmp(cfg.keeptrials,'no')
    error('cfg.keeptrials must be set to ''yes'' or ''no''.');
  end
end

% equate trial counts?
if ~isfield(cfg,'equatetrials')
  cfg.equatetrials = 'no';
  fprintf('Not equating trials across conditions.\n');
elseif isfield(cfg,'equatetrials')
  if ~strcmp(cfg.equatetrials,'yes') && ~strcmp(cfg.equatetrials,'no')
    error('cfg.equatetrials must be set to ''yes'' or ''no''.');
  end
end

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
      warning([mfilename,':rmevoked method'],'Both remove options set to ''no'', setting defaults.');
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
  error('Normalizing cfg.output=''coh'' is not allowed.');
end

if ~isfield(cfg,'baseline_type')
  cfg.baseline_type = '';
end

if strcmp(cfg.transform,'dB') && ~isempty(cfg.baseline_type) && ~strcmp(cfg.baseline_type,'relative')
  warning([mfilename,':dbNormBL'],'dB normalization uses a relative baseline type (division). Setting cfg.baseline_type=''relative''.\n');
  cfg.baseline_type = 'relative';
end

% if iscell(ana.eventValues)
%   if ~iscell(ana.eventValues{1})
%     ana.eventValues = {ana.eventValues};
%   end
% elseif ~iscell(ana.eventValues)
%   ana.eventValues = {{ana.eventValues}};
% end

%% if we're equating trials, find out the number of trials to choose

% TODO: see if randTrials is set up for multiple sessions

if strcmp(cfg.equatetrials,'yes')
  fprintf('Will be equating trials across conditions within a type.\n');
  if isfield(exper,'randTrials')
    fprintf('exper struct already contains randTrials field. Using that for sub-sampling trials.\n');
  else
    % initialize random number generator
    rng('shuffle','twister');
    
    for sub = 1:length(exper.subjects)
      for ses = 1:length(exper.sesStr)
        % initialize to find the lowest number of trials within each subject
        lowEvNum = Inf(length(exper.subjects),length(ana.eventValues{ses}));
        % initialize to store the randomly chosen events
        exper.(exper.sesStr{ses}).randTrials = cell(length(exper.subjects),length(ana.eventValues{ses}));
        for typ = 1:length(ana.eventValues{ses})
          for evVal = 1:length(ana.eventValues{ses}{typ})
            if exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub) < lowEvNum(sub,typ) && exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub) > 0
              lowEvNum(sub,typ) = exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub);
            end
          end
          exper.(exper.sesStr{ses}).randTrials{sub,typ} = nan(lowEvNum(sub,typ),length(ana.eventValues{ses}{typ}));
          for evVal = 1:length(ana.eventValues{ses}{typ})
            if exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub) ~= 0
              exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal) = sort(randperm(exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub),lowEvNum(sub,typ)));
            else
              exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal) = [];
            end
          end % evVal
        end % typ
      end % ses
    end % sub
  end
  
  % print out the trial counts for each subject
  for typ = 1:length(ana.eventValues{ses})
    fprintf('Condition type consisting of: %s\n',sprintf(repmat('%s ',1,length(ana.eventValues{ses}{typ})),ana.eventValues{ses}{typ}{:}));
    if length(exper.subjects{1}) > 7
      tabchar_sub = '\t';
    else
      tabchar_sub = '';
    end
    fprintf('Subject%s%s\n',sprintf(tabchar_sub),sprintf(repmat('\t%s',1,length(ana.eventValues{ses}{typ})),ana.eventValues{ses}{typ}{:}));
    for ses = 1:length(exper.sesStr)
      for sub = 1:length(exper.subjects)
        subStr = exper.subjects{sub};
        for evVal = 1:length(ana.eventValues{ses}{typ})
          if length(ana.eventValues{ses}{typ}{evVal}) > 7
            tabchar_ev = '\t';
          else
            tabchar_ev = '';
          end
          subStr = cat(2,subStr,sprintf('\t%d%s',exper.nTrials.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal})(sub),sprintf(tabchar_ev)));
        end
        fprintf('%s\n',subStr);
      end
    end
    fprintf('\n');
  end
end

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
    
    if strcmp(cfg.ftype,'raw')
      saveFileDir = fullfile(dirs.saveDirRaw,exper.subjects{sub},exper.sesStr{ses});
    else
      saveFileDir = fullfile(dirs.saveDirProc,exper.subjects{sub},exper.sesStr{ses});
    end
    
    for typ = 1:length(ana.eventValues{ses})
      for evVal = 1:length(ana.eventValues{ses}{typ})
        
        inputfile = fullfile(saveFileDir,sprintf('data_%s_%s.mat',cfg.ftype,ana.eventValues{ses}{typ}{evVal}));
        if exist(inputfile,'file')
          fprintf('Loading %s %s %s: %s...\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal},inputfile);
          % load the data
          subSesEvData = load(inputfile);
          % get the name of the field
          fn = fieldnames(subSesEvData);
          % rename the field to 'data'
          if length(fn) == 1
            
            % can't have values exactly equal to zero (why?)
            if powflg
              param = 'powspctrm';
              subSesEvData.(fn{1}).(param)(subSesEvData.(fn{1}).(param) == 0) = eps(0);
            end
            
            if csdflg
              param = 'crsspctrm';
            end
            
            % convert complex fourier-spectra to power or coherence
            if fftflg
              
              % TODO: do I need to update something in the FT structure now
              % that I get the subset of trials at the start? maybe use
              % ft_selectdata?
              %
              % equate the trial counts
              if strcmp(cfg.equatetrials,'yes')
                % get the subset of trials
                if size(subSesEvData.(fn{1}).(fourierparam),1) ~= length(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal))
                  fprintf('Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(fn{1}).(fourierparam),1));
                  %subSesEvData.(fn{1}).(fourierparam) = subSesEvData.(fn{1}).(fourierparam)(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal),:,:,:);
                  subSesEvData.(fn{1}) = ft_selectdata_old(subSesEvData.(fn{1}),'param',fourierparam,'rpt',exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal));
                end
              end
                
              fprintf('Converting %s to %s.\n',cfg.ftype,cfg.output);
              if strcmp(cfg.output,'pow')
                %% get the modulus of the whole data and evoked data
                param = 'modulus';
                
                % get the modulus for whole data (i.e., abs(whole fourier))
                subSesEvData.(fn{1}).(param) = abs(subSesEvData.(fn{1}).(fourierparam));
                
                % get the modulus for evoked data (i.e., abs(evoked fourier))
                if strcmp(cfg.rmevoked,'yes')
                  data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(param) = abs(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(fourierparam));
                end
                
                if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
                  % option 1: subtract evoked modulus from whole modulus;
                  % this will make evoked=positive and induced=negative
                  
                  %fprintf('Subtraction: (abs(whole fourier) - abs(evoked fourier)).^2 = induced power.\n');
                  fprintf('Subtraction: ((whole modulus) - (evoked modulus)).^2 = induced power.\n');
                  
                  % subtract evoked from single trials, then square
                  % % evoked_data = repmat(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(fourierparam),[size(subSesEvData.(fn{1}).(fourierparam),1),1,1,1]);
                  % % subSesEvData.(fn{1}).(powparam) = (abs(subSesEvData.(fn{1}).(fourierparam)) - abs(evoked_data)).^2;
                  evoked_data = repmat(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(param),[size(subSesEvData.(fn{1}).(param),1),1,1,1]);
                  subSesEvData.(fn{1}).(powparam) = (subSesEvData.(fn{1}).(param) - evoked_data).^2;
                  
                  % % subtract means of moduluses, then square
                  % % subSesEvData.(fn{1}).(powparam) = (mean(abs(subSesEvData.(fn{1}).(fourierparam)),1) - abs(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(fourierparam))).^2;
                  % subSesEvData.(fn{1}).(powparam) = (mean(subSesEvData.(fn{1}).(param),1) - data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(param)).^2;
                  
                  % % square means of moduluses, then subtract
                  % % subSesEvData.(fn{1}).(powparam) = mean(abs(subSesEvData.(fn{1}).(fourierparam)),1).^2 - abs(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(fourierparam)).^2;
                  % subSesEvData.(fn{1}).(powparam) = mean(subSesEvData.(fn{1}).(param),1).^2 - data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(param).^2;
                  
                elseif strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
                  % option 2: set up evoked fourier to subtract as power
                  if isfield(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,powparam)
                    data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = rmfield(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,powparam);
                  end
                  
                  % get the evoked power from evoked fourier data
                  data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(powparam) = data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(fourierparam).^2;
                  % no zeros
                  data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(powparam)(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(powparam) == 0) = eps(0);
                  % no longer need the fourierspctrm or modulus fields
                  data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = rmfield(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,{fourierparam param});
                  
                  % square the whole modulus to get power
                  subSesEvData.(fn{1}).(powparam) = subSesEvData.(fn{1}).(param).^2;
                elseif strcmp(cfg.rmevoked,'no')
                  %% otherwise just square the modulus to get power
                  
                  subSesEvData.(fn{1}).(powparam) = subSesEvData.(fn{1}).(param).^2;
                end
                
                % no zeros
                subSesEvData.(fn{1}).(powparam)(subSesEvData.(fn{1}).(powparam) == 0) = eps(0);
                % no longer need the fourierspctrm or modulus fields
                subSesEvData.(fn{1}) = rmfield(subSesEvData.(fn{1}),{fourierparam param});
                
                param = powparam;
                
              elseif strcmp(cfg.output,'coh')
                %% inter-trial phase coherence
                
                %if strcmp(cfg.equatetrials,'yes')
                %  % get the subset of trials
                %  fprintf('Inter-trial coherence data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(fn{1}).(fourierparam),1));
                %  subSesEvData.(fn{1}).(fourierparam) = subSesEvData.(fn{1}).(fourierparam)(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal),:,:,:);
                %end
                
                subSesEvData.(fn{1}).(cohparam) = subSesEvData.(fn{1}).(fourierparam) ./ abs(subSesEvData.(fn{1}).(fourierparam));
                subSesEvData.(fn{1}).(cohparam) = abs(squeeze(nanmean(subSesEvData.(fn{1}).(cohparam),1)));
                subSesEvData.(fn{1}) = rmfield(subSesEvData.(fn{1}),fourierparam);
                param = cohparam;
                % change dimord
                % TODO: can I automate this (e.g., with ft_checkdata)? or
                % maybe ft_freqdescriptives (keeptrials='no')
                subSesEvData.(fn{1}).dimord = 'chan_freq_time';
                subSesEvData.(fn{1}) = rmfield(subSesEvData.(fn{1}),'trialinfo');
              elseif strcmp(cfg.output,'phase')
                %% single-trial phase
                
                %if strcmp(cfg.equatetrials,'yes')
                %  % get the subset of trials
                %  fprintf('Phase data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(fn{1}).(fourierparam),1));
                %  subSesEvData.(fn{1}).(fourierparam) = subSesEvData.(fn{1}).(fourierparam)(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal),:,:,:);
                %end
                
                % use all the frequencies
                if isempty(cfg.phasefreq)
                  cfg.phasefreq = subSesEvData.(fn{1}).freq;
                end
                
                if ~isvector(cfg.phasefreq) || (size(cfg.phasefreq,1) == 1 && size(cfg.phasefreq,2) == 2)
                  meanfourier = nan(size(subSesEvData.(fn{1}).(fourierparam),1),length(cfg.phaseroi),size(cfg.phasefreq,1),size(subSesEvData.(fn{1}).(fourierparam),4));
                elseif iscolumn(cfg.phasefreq)
                  meanfourier = nan(size(subSesEvData.(fn{1}).(fourierparam),1),length(cfg.phaseroi),length(cfg.phasefreq),size(subSesEvData.(fn{1}).(fourierparam),4));
                end
                  
                for r = 1:length(cfg.phaseroi)
                  cfg.roi = cfg.phaseroi{r};
                  % set the channel information
                  if ismember(cfg.roi,ana.elecGroupsStr)
                    % if it's in the predefined ROIs, get the channel numbers
                    cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
                  else
                    % otherwise it should be the channel number(s) or 'all'
                    cfg.channel = cfg.roi;
                  end
                  chansel = ismember(subSesEvData.(fn{1}).label,cfg.channel);
                  
                  if ~isvector(cfg.phasefreq) || (size(cfg.phasefreq,1) == 1 && size(cfg.phasefreq,2) == 2)
                    for f = 1:size(cfg.phasefreq,1)
                      freqsel = subSesEvData.(fn{1}).freq >= cfg.phasefreq(f,1) & subSesEvData.(fn{1}).freq <= cfg.phasefreq(f,2);
                      meanfourier(:,r,f,:) = nanmean(nanmean(subSesEvData.(fn{1}).(fourierparam)(:,chansel,freqsel,:),3),2);
                    end
                  elseif iscolumn(cfg.phasefreq)
                    meanfourier(:,r,:,:) = nanmean(subSesEvData.(fn{1}).(fourierparam)(:,chansel,:,:),2);
                  end
                end
                subSesEvData.(fn{1}).(fourierparam) = meanfourier;
                
                % calculate phase angle
                subSesEvData.(fn{1}).(phaseparam) = angle(subSesEvData.(fn{1}).(fourierparam));
                subSesEvData.(fn{1}) = rmfield(subSesEvData.(fn{1}),fourierparam);
                
                % print it
                for r = 1:length(cfg.phaseroi)
                  cfg.roi = cfg.phaseroi{r};
                  % set the channel information
                  if ismember(cfg.roi,ana.elecGroupsStr)
                    % if it's in the predefined ROIs, get the channel numbers
                    cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
                    % set the string for the filename
                    chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
                  else
                    % otherwise it should be the channel number(s) or 'all'
                    cfg.channel = cfg.roi;
                    
                    % set the string for the filename
                    if isfield(cfg,'cohrefchannel')
                      chan_str = [cfg.cohrefchannel,'-',sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:})];
                    else
                      chan_str = sprintf(repmat('%s_',1,length(cfg.roi)),cfg.roi{:});
                    end
                  end
                  
                  for f = 1:size(subSesEvData.(fn{1}).(phaseparam),3)
                    if ~isvector(cfg.phasefreq) || length(cfg.phasefreq) == 2
                      freq_str = sprintf('%.1f-%.1f',cfg.phasefreq(f,1),cfg.phasefreq(f,2));
                    else
                      freq_str = sprintf('%.1f',cfg.phasefreq(f));
                    end
                    figure
                    surf(subSesEvData.(fn{1}).time,1:size(subSesEvData.(fn{1}).(phaseparam),1),squeeze(subSesEvData.(fn{1}).(phaseparam)(:,r,f,:)));
                    shading interp;view([0,90]);axis tight;
                    title(sprintf('%s, %s, %s, %s, %s Hz',exper.subjects{sub},strrep(exper.sessions{ses},'_',''),ana.eventValues{ses}{typ}{evVal},strrep(chan_str,'_',' '),freq_str));
                    xlabel('Time (s)');
                    ylabel('Trial number');
                    publishfig(gcf,0);
                    figfilename = sprintf('phase_%s_%s_%s_%s%sHz.png',exper.subjects{sub},exper.sessions{ses},ana.eventValues{ses}{typ}{evVal},chan_str,freq_str);
                    dirs.saveDirFigsPhase = fullfile(dirs.saveDirFigs,'tfr_phase');
                    if ~exist(dirs.saveDirFigsPhase,'dir')
                      mkdir(dirs.saveDirFigsPhase)
                    end
                    print(gcf,'-dpng','-r150',fullfile(dirs.saveDirFigsPhase,figfilename));
                    close all
                  end % f
                end % r
                
              end % cfg.output
            end % fftflg
            
            %% transform the data
            if ~isempty(cfg.transform)
              fprintf('Doing a %s transformation.\n',cfg.transform);
              
              if strcmp(cfg.transform,'vec')
                for tr = 1:size(subSesEvData.(fn{1}).(param),1)
                  for ch = 1:size(subSesEvData.(fn{1}).(param),2)
                    % doesn't ignore NaNs, but is way faster
                    subSesEvData.(fn{1}).(param)(tr,ch,:,:) = normr(squeeze(subSesEvData.(fn{1}).(param)(tr,ch,:,:)));
                    % % try to ignore NaNs
                    %for fq = 1:size(subSesEvData.(fn{1}).(param),3)
                    %  thisTrl = squeeze(subSesEvData.(fn{1}).(param)(tr,ch,fq,:));
                    %  subSesEvData.(fn{1}).(param)(tr,ch,fq,~isnan(thisTrl)) = normc(thisTrl(~isnan(thisTrl)));
                    %end
                    % % trying another method, but it doesn't work
                    %thisTrl = squeeze(subSesEvData.(fn{1}).(param)(tr,ch,:,:));
                    %subSesEvData.(fn{1}).(param)(tr,ch,:,~isnan(thisTrl)) = normr(thisTrl(~isnan(thisTrl)));
                  end
                end
              else
                % run using feval (e.g., for log10 normalization)
                subSesEvData.(fn{1}).(param) = feval(cfg.transform,subSesEvData.(fn{1}).(param));
                
                % option 2: normalize the evoked data
                if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
                  data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(powparam) = feval(cfg.transform,data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(powparam));
                end
              end
            end
            
            % option 2: subtract evoked power from whole power; need to
            % do this after log normalization
            if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
              fprintf('Subtracting evoked power from whole power to get induced power.\n');
              evoked_data = repmat(data_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.(powparam),[size(subSesEvData.(fn{1}).(powparam),1),1,1,1]);
              subSesEvData.(fn{1}).(powparam) = subSesEvData.(fn{1}).(powparam) - evoked_data;
            end
            
            %% baseline correct power
            
            % single trial normalization uses the entire trial for
            % normalization
            if strcmp(cfg.baseline_data,'pow')
              if ~isempty(cfg.baseline_type) && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow'))
                
                if ndims(subSesEvData.(fn{1}).(param)) ~= 4
                  error('data (%s) does not have four dimensions',param);
                end
                
                if strcmp(cfg.norm_trials,'single')
                  
                  % single-trial normalization uses all the time points
                  blt = true(size(subSesEvData.(fn{1}).time));
                  
                  % mean across baseline period
                  blm = nanmean(subSesEvData.(fn{1}).(param)(:,:,:,blt),4);
                  
                  if strcmp(cfg.baseline_type,'zscore')
                    fprintf('Z-transforming each trial using entire trial as baseline.\n');
                    % std across time
                    blstd = nanstd(subSesEvData.(fn{1}).(param)(:,:,:,blt),0,4);
                    
                    % z-transform
                    subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(fn{1}).(param),blm),blstd);
                    
                  elseif strcmp(cfg.baseline_type,'dB')
                    fprintf('Converting each trial to dB using entire trial as baseline.\n');
                    % divide by baseline mean and convert to dB
                    subSesEvData.(fn{1}).(param) = 10*log10(bsxfun(@rdivide,subSesEvData.(fn{1}).(param),blm));
                    
                    % mean(10*log10) is not the same as 10*log10(mean)
                    
                    % debug
                    %figure;imagesc(squeeze(mean(subSesEvData.(fn{1}).(param)(:,62,:,:),1)));axis xy;colorbar
                    
                  elseif strcmp(cfg.baseline_type,'absolute')
                    fprintf('Subtracting mean(entire trial) power from entire trial.\n');
                    subSesEvData.(fn{1}).(param) = bsxfun(@minus,subSesEvData.(fn{1}).(param),blm);
                    
                  elseif strcmp(cfg.baseline_type,'relative')
                    fprintf('Dividing entire trial power by mean(entire trial).\n');
                    subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,subSesEvData.(fn{1}).(param),blm);
                    
                  elseif strcmp(cfg.baseline_type,'relchange')
                    fprintf('Subtracting mean(entire trial) power from entire trial and dividing by that mean.\n');
                    subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(fn{1}).(param),blm),blm);
                  end
                  
                end
                
                % now compare it to the baseline activity
                if strcmp(cfg.norm_trials,'average') || strcmp(cfg.norm_trials,'single')
                  
                  % get the baseline time indices
                  blt = subSesEvData.(fn{1}).time >= cfg.baseline_time(1) & subSesEvData.(fn{1}).time <= cfg.baseline_time(2);
                  
                  % mean across baseline period
                  blm = nanmean(subSesEvData.(fn{1}).(param)(:,:,:,blt),4);
                  
                  if strcmp(cfg.baseline_type,'zscore')
                    fprintf('Z-transforming data relative to mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline_time(1),cfg.baseline_time(2));
                    % std across time, then avg across events (lower freqs often get smaller std)
                    blstd = nanmean(nanstd(subSesEvData.(fn{1}).(param)(:,:,:,blt),0,4),1);
                    
                    % % avg across time, then std across events (higher freqs often get smaller std)
                    % blstd = nanstd(nanmean(subSesEvData.(fn{1}).(param)(:,:,:,blt),4),0,1);
                    
                    % % concatenate all times of all events and std (equivalent std across freqs)
                    %blstd = shiftdim(squeeze(std(double(reshape(shiftdim(subSesEvData.(fn{1}).(param)(:,:,:,blt),3),...
                    %  size(subSesEvData.(fn{1}).(param)(:,:,:,blt),1)*size(subSesEvData.(fn{1}).(param)(:,:,:,blt),4),...
                    %  size(subSesEvData.(fn{1}).(param)(:,:,:,blt),2),size(subSesEvData.(fn{1}).(param)(:,:,:,blt),3))),0,1)),-1);
                    
                    subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(fn{1}).(param),blm),blstd);
                    
                  elseif strcmp(cfg.baseline_type,'dB')
                    fprintf('Converting to dB relative to mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline_time(1),cfg.baseline_time(2));
                    % divide by baseline mean and convert to dB
                    
                    subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,subSesEvData.(fn{1}).(param),blm);
                    if ~strcmp(cfg.norm_trials,'single')
                      subSesEvData.(fn{1}).(param) = 10*log10(subSesEvData.(fn{1}).(param));
                    end
                    
                    % TODO: mean(10*log10) is not the same as 10*log10(mean)
                    
                    %subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,subSesEvData.(fn{1}).(param),blm);
                    
                  elseif strcmp(cfg.baseline_type,'absolute')
                    fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
                    subSesEvData.(fn{1}).(param) = bsxfun(@minus,subSesEvData.(fn{1}).(param),blm);
                    
                  elseif strcmp(cfg.baseline_type,'relative')
                    fprintf('Dividing entire trial power by mean([%.2f %.2f] pre-stimulus) power.\n',cfg.baseline_time(1),cfg.baseline_time(2));
                    subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,subSesEvData.(fn{1}).(param),blm);
                    
                  elseif strcmp(cfg.baseline_type,'relchange')
                    fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial and dividing by that mean.\n',cfg.baseline_time(1),cfg.baseline_time(2));
                    subSesEvData.(fn{1}).(param) = bsxfun(@rdivide,bsxfun(@minus,subSesEvData.(fn{1}).(param),blm),blm);
                  end

                end % norm_trials
              end % baseline_type
            end % baseline_data
            
%             if strcmp(cfg.baseline_data,'pow')
%               if ~isempty(cfg.baseline_type) && ~strcmp(cfg.transform,'dB')
%                 
%                 % get the baseline time indices
%                 blt = subSesEvData.(fn{1}).time >= cfg.baseline_time(1) & subSesEvData.(fn{1}).time <= cfg.baseline_time(2);
%                 
%                 if strcmp(cfg.baseline_type,'zscore') && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(fn{1}).(param)) == 4
%                   fprintf('Z-transforming data relative to mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline_time(1),cfg.baseline_time(2));
%                   
%                   nTrl = size(subSesEvData.(fn{1}).(param),1);
%                   %nFrq = size(subSesEvData.(fn{1}).(param),3);
%                   nSmp = size(subSesEvData.(fn{1}).(param),4);
%                   
%                   blm = repmat(nanmean(subSesEvData.(fn{1}).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
%                   blstd = repmat(nanmean(nanstd(subSesEvData.(fn{1}).(param)(:,:,:,blt),0,1),4),[nTrl,1,1,nSmp]);
%                   subSesEvData.(fn{1}).(param) = (subSesEvData.(fn{1}).(param) - blm) ./ blstd;
%                   
%                 elseif (strcmp(cfg.baseline_type,'absolute') || strcmp(cfg.baseline_type,'relative') || strcmp(cfg.baseline_type,'relchange')) && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(fn{1}).(param)) == 4
%                   %blt = subSesEvData.(fn{1}).time >= cfg.baseline_time(1) & subSesEvData.(fn{1}).time <= cfg.baseline_time(2);
%                   
%                   nSmp = size(subSesEvData.(fn{1}).(param),4);
%                   
%                   blm = repmat(nanmean(subSesEvData.(fn{1}).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
%                   
%                   if strcmp(cfg.baseline_type,'absolute')
%                     fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
%                     subSesEvData.(fn{1}).(param) = subSesEvData.(fn{1}).(param) - blm;
%                   elseif strcmp(cfg.baseline_type,'relative')
%                     error('relative baseline is not working');
%                     %fprintf('Dividing power by mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline_time(1),cfg.baseline_time(2));
%                     %subSesEvData.(fn{1}).(param) = subSesEvData.(fn{1}).(param) ./ blm;
%                   elseif strcmp(cfg.baseline_type,'relchange')
%                     error('relchange baseline is not working');
%                     %fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial and dividing by the mean.\n',cfg.baseline_time(1),cfg.baseline_time(2));
%                     %subSesEvData.(fn{1}).(param) = (subSesEvData.(fn{1}).(param) - blm) ./ blm;
%                   end
%                   
%                 elseif strcmp(cfg.baseline_type,'absolute') && strcmp(cfg.output,'coh') && ndims(subSesEvData.(fn{1}).(param)) == 3
%                   % not sure that coherence gets baseline correction???
%                   
%                   fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) coherence from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
%                   %blt = subSesEvData.(fn{1}).time >= cfg.baseline_time(1) & subSesEvData.(fn{1}).time <= cfg.baseline_time(2);
%                   
%                   nSmp = size(subSesEvData.(fn{1}).(param),3);
%                   
%                   blm = repmat(nanmean(subSesEvData.(fn{1}).(param)(:,:,blt),3),[1,1,nSmp]);
%                   subSesEvData.(fn{1}).(param) = (subSesEvData.(fn{1}).(param) - blm);
%                 end
%               end
%             end
            
            %% debug
%             chan = 62;
%             chan = 25;
%             clim = [-1 1];
%             %clim = [-3 3];
%             if strcmp(cfg.transform,'dB')
%               clim = clim .* 10;
%             end
%             figure;surf(subSesEvData.(fn{1}).time,subSesEvData.(fn{1}).freq,squeeze(mean(subSesEvData.(fn{1}).(powparam)(:,chan,:,:),1)));shading interp;view([0,90]);axis tight;colorbar;
%             title(sprintf('chan %d: rmevokedfourier: %s; rmevokedpow: %s',chan,cfg.rmevokedfourier,cfg.rmevokedpow));
%             caxis(clim)
%             figure;imagesc(subSesEvData.(fn{1}).time,subSesEvData.(fn{1}).freq,squeeze(mean(subSesEvData.(fn{1}).(powparam)(:,chan,:,:),1)),clim);axis xy;colorbar
%             title(sprintf('chan %d: rmevokedfourier: %s; rmevokedpow: %s',chan,cfg.rmevokedfourier,cfg.rmevokedpow));
%             
%             figure
%             for tr = 1:size(subSesEvData.(fn{1}).(powparam),1)
%               clf
%               surf(subSesEvData.(fn{1}).time,subSesEvData.(fn{1}).freq,squeeze(subSesEvData.(fn{1}).(powparam)(tr,chan,:,:)));shading interp;view([0,90]);axis tight;colorbar;
%               title(sprintf('chan %d: trial %d/%d',chan,tr,size(subSesEvData.(fn{1}).(powparam),1)));
%               caxis([-3 3])
%               keyboard
%             end
            
            %% save the data in a container struct
            fprintf('%s %s %s: ',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
            if (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && isfield(subSesEvData.(fn{1}),powparam) && ndims(subSesEvData.(fn{1}).(powparam)) == 4
              fprintf('Power data: ');
              %             elseif strcmp(cfg.keeptrials,'no') && (~strcmp(cfg.ftype,'pow') || ~strcmp(cfg.output,'pow')) && isfield(subSesEvData.(fn{1}),'powspctrm') && ndims(subSesEvData.(fn{1}).(powparam)) == 4
              %               error('\n%s %s %s: Can only keep trials for cfg.ftype=''pow''. You set it to ''%s''.\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal},cfg.ftype);
              %             elseif strcmp(cfg.keeptrials,'no') && ~isfield(subSesEvData.(fn{1}),'powspctrm')
              %               error('\n%s %s %s: Can only keep trials with ''powspctrm'' field. Please examine your data.\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
              %             elseif strcmp(cfg.keeptrials,'no') && isfield(subSesEvData.(fn{1}),'powspctrm') && ndims(subSesEvData.(fn{1}).(powparam)) ~= 4
              %               error('\n%s %s %s: Can only keep trials for ndims(powspctrm)==4. This data has ndims=%d.\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal},ndims(subSesEvData.(fn{1}).(powparam)));
              
              % set the configuration of ft_freqdescriptives
              cfg_fd = [];
              cfg_fd.keeptrials = cfg.keeptrials;
              
              if strcmp(cfg.keeptrials,'no')
                fprintf('Averaging over individual trials.\n');
              elseif strcmp(cfg.keeptrials,'yes')
                fprintf('Keeping individual trials.\n');
              end
              
%               if strcmp(cfg.equatetrials,'yes')
%                 fprintf(' Equating trial counts across all conditions (within a type).\n');
%                 fprintf('Randomly selecting %d of %d trials.\n',length(exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal)),size(subSesEvData.(fn{1}).(powparam),1));
%                 % get the subset of trials
%                 cfg_fd.trials = exper.(exper.sesStr{ses}).randTrials{sub,typ}(:,evVal);
%               elseif strcmp(cfg.equatetrials,'no')
%                 fprintf(' Using all trials.\n');
%                 % use all of the trials
%                 cfg_fd.trials = 'all';
%               end
              
              % Do a very coarse rejection for trials with artifacts
              if ~isempty(cfg.zthresh) && strcmp(cfg.baseline_type,'zscore')
                % reject trials with huge zscores
                cfg_fd.trials = true(1,size(subSesEvData.(fn{1}).(param),1));
                [tr,ch] = find(squeeze(nanmean(nanmean(subSesEvData.(fn{1}).(param),4),3)) > cfg.zthresh);
                if ~isempty(tr)
                  fprintf('Rejecting %d of %d trials due to zscore > %.2f in %d channels (NB: this is a very coarse rejection).\n',length(unique(tr)),length(cfg_fd.trials),cfg.zthresh,length(unique(ch)));
                  cfg_fd.trials(unique(tr)) = false;
                  fprintf('Keeping %d of %d trials.\n',sum(ismember(cfg_fd.trials,1)),length(cfg_fd.trials));
                else
                  fprintf('Keeping all %d trials.\n',size(subSesEvData.(fn{1}).(param),1));
                  cfg_fd.trials = 'all';
                end
              else
                % TODO: do I need to update something in the FT structure now
                % that I get the subset of trials at the start?
                cfg_fd.trials = 'all';
              end
              
              % make sure there are trials left to average
              if (isvector(cfg_fd.trials) && ismember(1,cfg_fd.trials)) || (ischar(cfg_fd.trials) && strcmp(cfg_fd.trials,'all'))
                  % use ft_freqdescriptives to average over individual
                  % trials, if desired and the data is appropriate
                  data.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_freqdescriptives(cfg_fd,subSesEvData.(fn{1}));
              else %if (isvector(cfg_fd.trials) && length(unique(cfg_fd.trials)) == 1 && unique(cfg_fd.trials) == 0)
                  warning([mfilename,':allTrialsRejected'],'ALL TRIALS REJECTED: %s\n',ana.eventValues{ses}{typ}{evVal});
                  fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n\n');
                  data.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.cfg.trl = [];
              end % if

              
            elseif strcmp(cfg.output,'coh') && isfield(subSesEvData.(fn{1}),cohparam) && ndims(subSesEvData.(fn{1}).(cohparam)) == 3
              fprintf('Phase coherence data: Individual trials are lost.');
              if strcmp(cfg.equatetrials,'no')
                fprintf(' Using all trials. NB: this could affect the coherence calculation if trial counts are unequal!!\n');
              %elseif strcmp(cfg.equatetrials,'yes')
              %  fprintf(' Using randomly selected subsets of trials so each condition (within a type) has the same number of trials.\n');
              end
              data.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = subSesEvData.(fn{1});
            elseif strcmp(cfg.output,'phase') && isfield(subSesEvData.(fn{1}),phaseparam) && ndims(subSesEvData.(fn{1}).(phaseparam)) == 4
              data.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = subSesEvData.(fn{1});
            else
              fprintf('Need to figure out what to do for this case.\n');
              keyboard
            end
            
          else
            error('More than one field in data struct! There should only be one.\n');
          end
          
          % % old method
          % data.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub) = load(inputfile);
          fprintf('Done.\n\n');
        else
          warning([mfilename,':noFileFound'],'NOT FOUND: %s\n',inputfile);
          fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n\n');
          data.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.cfg.trl = [];
        end % if
      end % for evVal
    end % for typ
  end % for ses
end % for sub

end
