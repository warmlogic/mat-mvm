function [data,exper] = mm_ft_loadData(cfg,exper,dirs,ana,data_evoked)
%MM_FT_LOADDATA Load subject data into a full struct
%
% [data,exper] = mm_ft_loadData(cfg,exper,ana,dirs,ana,data_evoked)
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
%  - for cfg.baselinetype='condition', use ft_freqcomparison
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

if ~isfield(cfg,'normalize')
  cfg.normalize = '';
elseif ~isempty(cfg.normalize) && strcmp(cfg.output,'coh')
  error('Normalizing cfg.output=''coh'' is not allowed.');
end

if ~isfield(cfg,'baselinetype')
  cfg.baselinetype = '';
end

if strcmp(cfg.normalize,'dB') && ~isempty(cfg.baselinetype) && ~strcmp(cfg.baselinetype,'relative')
  warning([mfilename,':dbNormBL'],'dB normalization uses a relative baseline type (division). Setting cfg.baselinetype=''relative''.\n');
  cfg.baselinetype = 'relative';
end

if iscell(ana.eventValues)
  if ~iscell(ana.eventValues{1})
    ana.eventValues = {ana.eventValues};
  end
elseif ~iscell(ana.eventValues)
  ana.eventValues = {{ana.eventValues}};
end

%% if we're equating trials, find out the number of trials to choose

if strcmp(cfg.equatetrials,'yes')
  fprintf('Will be equating trials across conditions within a type.\n');
  if isfield(exper,'randTrials')
    fprintf('exper struct already contains randTrials field. Using that for sub-sampling trials.\n');
  else
    % initialize random number generator
    rng('shuffle','twister');
    % initialize to find the lowest number of trials within each subject
    lowEvNum = Inf(length(exper.subjects),length(exper.sessions),length(ana.eventValues));
    % initialize to store the randomly chosen events
    exper.randTrials = cell(length(exper.subjects),length(exper.sessions),length(ana.eventValues));
    
    for sub = 1:length(exper.subjects)
      for ses = 1:length(exper.sessions)
        for typ = 1:length(ana.eventValues)
          for evVal = 1:length(ana.eventValues{typ})
            if exper.nTrials.(ana.eventValues{typ}{evVal})(sub,ses) < lowEvNum(sub,ses,typ) && exper.nTrials.(ana.eventValues{typ}{evVal})(sub,ses) > 0
              lowEvNum(sub,ses,typ) = exper.nTrials.(ana.eventValues{typ}{evVal})(sub,ses);
            end
          end
          exper.randTrials{sub,ses,typ} = nan(lowEvNum(sub,ses,typ),length(ana.eventValues{typ}));
          for evVal = 1:length(ana.eventValues{typ})
            if exper.nTrials.(ana.eventValues{typ}{evVal})(sub,ses) ~= 0
              exper.randTrials{sub,ses,typ}(:,evVal) = sort(randperm(exper.nTrials.(ana.eventValues{typ}{evVal})(sub,ses),lowEvNum(sub,ses,typ)));
            else
              exper.randTrials{sub,ses,typ}(:,evVal) = [];
            end
          end % evVal
        end % typ
      end % ses
    end % sub
  end
  
  % print out the trial counts for each subject
  for typ = 1:length(ana.eventValues)
    fprintf('Condition type consisting of: %s\n',sprintf(repmat('%s ',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:}));
    if length(exper.subjects{1}) > 7
      tabchar_sub = '\t';
    else
      tabchar_sub = '';
    end
    fprintf('Subject%s%s\n',sprintf(tabchar_sub),sprintf(repmat('\t%s',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:}));
    for ses = 1:length(exper.sessions)
      for sub = 1:length(exper.subjects)
        subStr = exper.subjects{sub};
        for evVal = 1:length(ana.eventValues{typ})
          if length(ana.eventValues{typ}{evVal}) > 7
            tabchar_ev = '\t';
          else
            tabchar_ev = '';
          end
          subStr = cat(2,subStr,sprintf('\t%d%s',exper.nTrials.(ana.eventValues{typ}{evVal})(sub,ses),sprintf(tabchar_ev)));
        end
        fprintf('%s\n',subStr);
      end
    end
    fprintf('\n');
  end
end

%% process the data

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    % turn the session name into a string for easier printing
    if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
      sesStr = exper.sessions{ses}{1};
      for i = 2:length(exper.sessions{ses})
        sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
      end
    elseif ~iscell(exper.sessions{ses}) || (iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1)
      sesStr = exper.sessions{ses};
    end
    
    if strcmp(cfg.ftype,'raw')
      saveFileDir = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
    else
      saveFileDir = fullfile(dirs.saveDirProc,exper.subjects{sub},sesStr);
    end
    
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        
        inputfile = fullfile(saveFileDir,sprintf('data_%s_%s.mat',cfg.ftype,ana.eventValues{typ}{evVal}));
        if exist(inputfile,'file')
          fprintf('Loading %s %s %s: %s...\n',exper.subjects{sub},sesStr,ana.eventValues{typ}{evVal},inputfile);
          % load the data
          subSesEvData = load(inputfile);
          % get the name of the field
          fn = fieldnames(subSesEvData);
          % rename the field to 'data'
          if length(fn) == 1
            
            % can't have values exactly equal to zero (why?)
            if powflg
              param = 'powspctrm';
              subSesEvData.(cell2mat(fn)).(param)(subSesEvData.(cell2mat(fn)).(param) == 0) = eps(0);
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
                fprintf('Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.randTrials{sub,ses,typ}(:,evVal)),size(subSesEvData.(cell2mat(fn)).(fourierparam),1));
                subSesEvData.(cell2mat(fn)).(fourierparam) = subSesEvData.(cell2mat(fn)).(fourierparam)(exper.randTrials{sub,ses,typ}(:,evVal),:,:,:);
              end
                
              fprintf('Converting %s to %s.\n',cfg.ftype,cfg.output);
              if strcmp(cfg.output,'pow')
                %% get the modulus of the whole data and evoked data
                param = 'modulus';
                
                % get the modulus for whole data
                subSesEvData.(cell2mat(fn)).(param) = abs(subSesEvData.(cell2mat(fn)).(fourierparam));
                
                % get the modulus for evoked data
                if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
                  data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = abs(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(fourierparam));
                end
                
                %% baseline correct the modulus
                if strcmp(cfg.baselinedata,'mod')
                  if ~isempty(cfg.baselinetype)
                    % find the baseline time indices
                    blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                    
                    if strcmp(cfg.baselinetype,'zscore')
                      fprintf('Z-transforming data relative to mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline(1),cfg.baseline(2));
                      
                      nTrl = size(subSesEvData.(cell2mat(fn)).(param),1);
                      %nFrq = size(subSesEvData.(cell2mat(fn)).(param),3);
                      nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                      
                      % get baseline mean and standard deviation
                      blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                      blstd = repmat(nanmean(nanstd(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),0,1),4),[nTrl,1,1,nSmp]);
                      
                      % z-transform
                      subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm) ./ blstd;
                      
                      % can't z-transform evoked data: there's only 1 "trial"
                      
                    elseif (strcmp(cfg.baselinetype,'absolute') || strcmp(cfg.baselinetype,'relative') || strcmp(cfg.baselinetype,'relchange')) && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(cell2mat(fn)).(param)) == 4
                      % find the baseline time indices
                      %blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                      
                      nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                      % get baseline mean
                      blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                      
                      % get the details for evoked data
                      if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
                        blt_ev = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.time >= cfg.baseline(1) & data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.time <= cfg.baseline(2);
                        nSmp_ev = size(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param),4);
                        blm_ev = repmat(nanmean(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,:,:,blt_ev),4),[1,1,1,nSmp_ev]);
                      end
                      
                      if strcmp(cfg.baselinetype,'absolute')
                        % subtraction
                        fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial.\n',cfg.baseline(1),cfg.baseline(2));
                        
                        % whole
                        subSesEvData.(cell2mat(fn)).(param) = subSesEvData.(cell2mat(fn)).(param) - blm;
                        
                        % evoked
                        if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
                          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) - blm_ev;
                        end
                      elseif strcmp(cfg.baselinetype,'relative')
                        %error('relative baseline is not working');
                        
                        % division
                        fprintf('Dividing power by mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline(1),cfg.baseline(2));
                        
                        % whole
                        subSesEvData.(cell2mat(fn)).(param) = subSesEvData.(cell2mat(fn)).(param) ./ blm;
                        
                        % evoked
                        if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
                          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) ./ blm_ev;
                        end
                      elseif strcmp(cfg.baselinetype,'relchange')
                        %error('relchange baseline is not working');
                        
                        % subtraction then division
                        fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial and dividing by the mean.\n',cfg.baseline(1),cfg.baseline(2));
                        
                        % whole
                        subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm) ./ blm;
                        
                        % evoked
                        if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
                          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = (data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) - blm_ev) ./ blm_ev;
                        end
                      end
                    end
                  end % if baseline
                end
                
                
                if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedfourier,'yes')
                  % option 1: subtract evoked modulus from whole modulus;
                  % this will make evoked=positive and induced=negative
                  
                  %fprintf('Subtraction: (abs(whole fourier) - abs(evoked fourier)).^2 = induced power.\n');
                  fprintf('Subtraction: ((whole modulus) - (evoked modulus)).^2 = induced power.\n');
                  
                  % subtract evoked from single trials, then square
                  % % evoked_data = repmat(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(fourierparam),[size(subSesEvData.(cell2mat(fn)).(fourierparam),1),1,1,1]);
                  % % subSesEvData.(cell2mat(fn)).(powparam) = (abs(subSesEvData.(cell2mat(fn)).(fourierparam)) - abs(evoked_data)).^2;
                  evoked_data = repmat(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param),[size(subSesEvData.(cell2mat(fn)).(param),1),1,1,1]);
                  subSesEvData.(cell2mat(fn)).(powparam) = (subSesEvData.(cell2mat(fn)).(param) - evoked_data).^2;
                  
                  % % subtract means of moduluses, then square
                  % % subSesEvData.(cell2mat(fn)).(powparam) = (mean(abs(subSesEvData.(cell2mat(fn)).(fourierparam)),1) - abs(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(fourierparam))).^2;
                  % subSesEvData.(cell2mat(fn)).(powparam) = (mean(subSesEvData.(cell2mat(fn)).(param),1) - data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)).^2;
                  
                  % % square means of moduluses, then subtract
                  % % subSesEvData.(cell2mat(fn)).(powparam) = mean(abs(subSesEvData.(cell2mat(fn)).(fourierparam)),1).^2 - abs(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(fourierparam)).^2;
                  % subSesEvData.(cell2mat(fn)).(powparam) = mean(subSesEvData.(cell2mat(fn)).(param),1).^2 - data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param).^2;
                  
                elseif strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
                  % option 2: set up evoked fourier to subtract as power
                  if isfield(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,powparam)
                    data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,powparam);
                  end
                  
                  % get the evoked power from evoked fourier data
                  data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(powparam) = abs(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(fourierparam)).^2;
                  % no zeros
                  data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(powparam)(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(powparam) == 0) = eps(0);
                  % no longer need the fourierspctrm field
                  data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,fourierparam);
                elseif strcmp(cfg.rmevoked,'no')
                  %% otherwise just square the modulus to get power
                  
                  %subSesEvData.(cell2mat(fn)).(powparam) = abs(subSesEvData.(cell2mat(fn)).(fourierparam)).^2;
                  subSesEvData.(cell2mat(fn)).(powparam) = subSesEvData.(cell2mat(fn)).(param).^2;
                end
                
                % no zeros
                subSesEvData.(cell2mat(fn)).(powparam)(subSesEvData.(cell2mat(fn)).(powparam) == 0) = eps(0);
                % no longer need the fourierspctrm field
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),fourierparam);
                % no longer need the modulus field
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),param);
                
                param = powparam;
                % TODO? tell FT that we have power?
                % subSesEvData.(cell2mat(fn)) = ft_checkdata(subSesEvData.(cell2mat(fn)),'haspow','yes');
              elseif strcmp(cfg.output,'coh')
                %% inter-trial phase coherence
                
                %if strcmp(cfg.equatetrials,'yes')
                %  % get the subset of trials
                %  fprintf('Inter-trial coherence data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.randTrials{sub,ses,typ}(:,evVal)),size(subSesEvData.(cell2mat(fn)).(fourierparam),1));
                %  subSesEvData.(cell2mat(fn)).(fourierparam) = subSesEvData.(cell2mat(fn)).(fourierparam)(exper.randTrials{sub,ses,typ}(:,evVal),:,:,:);
                %end
                
                subSesEvData.(cell2mat(fn)).(cohparam) = subSesEvData.(cell2mat(fn)).(fourierparam) ./ abs(subSesEvData.(cell2mat(fn)).(fourierparam));
                subSesEvData.(cell2mat(fn)).(cohparam) = abs(squeeze(nanmean(subSesEvData.(cell2mat(fn)).(cohparam),1)));
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),fourierparam);
                param = cohparam;
                % change dimord
                % TODO: can I automate this (e.g., with ft_checkdata)? or
                % maybe ft_freqdescriptives (keeptrials='no')
                subSesEvData.(cell2mat(fn)).dimord = 'chan_freq_time';
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),'trialinfo');
              elseif strcmp(cfg.output,'phase')
                %% single-trial phase
                
                %if strcmp(cfg.equatetrials,'yes')
                %  % get the subset of trials
                %  fprintf('Phase data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.randTrials{sub,ses,typ}(:,evVal)),size(subSesEvData.(cell2mat(fn)).(fourierparam),1));
                %  subSesEvData.(cell2mat(fn)).(fourierparam) = subSesEvData.(cell2mat(fn)).(fourierparam)(exper.randTrials{sub,ses,typ}(:,evVal),:,:,:);
                %end
                
                % use all the frequencies
                if isempty(cfg.phasefreq)
                  cfg.phasefreq = subSesEvData.(cell2mat(fn)).freq;
                end
                
                if ~isvector(cfg.phasefreq) || (size(cfg.phasefreq,1) == 1 && size(cfg.phasefreq,2) == 2)
                  meanfourier = nan(size(subSesEvData.(cell2mat(fn)).(fourierparam),1),length(cfg.phaseroi),size(cfg.phasefreq,1),size(subSesEvData.(cell2mat(fn)).(fourierparam),4));
                elseif iscolumn(cfg.phasefreq)
                  meanfourier = nan(size(subSesEvData.(cell2mat(fn)).(fourierparam),1),length(cfg.phaseroi),length(cfg.phasefreq),size(subSesEvData.(cell2mat(fn)).(fourierparam),4));
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
                  chansel = ismember(subSesEvData.(cell2mat(fn)).label,cfg.channel);
                  
                  if ~isvector(cfg.phasefreq) || (size(cfg.phasefreq,1) == 1 && size(cfg.phasefreq,2) == 2)
                    for f = 1:size(cfg.phasefreq,1)
                      freqsel = subSesEvData.(cell2mat(fn)).freq >= cfg.phasefreq(f,1) & subSesEvData.(cell2mat(fn)).freq <= cfg.phasefreq(f,2);
                      meanfourier(:,r,f,:) = nanmean(nanmean(subSesEvData.(cell2mat(fn)).(fourierparam)(:,chansel,freqsel,:),3),2);
                    end
                  elseif iscolumn(cfg.phasefreq)
                    meanfourier(:,r,:,:) = nanmean(subSesEvData.(cell2mat(fn)).(fourierparam)(:,chansel,:,:),2);
                  end
                end
                subSesEvData.(cell2mat(fn)).(fourierparam) = meanfourier;
                
                % calculate phase angle
                subSesEvData.(cell2mat(fn)).(phaseparam) = angle(subSesEvData.(cell2mat(fn)).(fourierparam));
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),fourierparam);
                
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
                  
                  for f = 1:size(subSesEvData.(cell2mat(fn)).(phaseparam),3)
                    if ~isvector(cfg.phasefreq) || length(cfg.phasefreq) == 2
                      freq_str = sprintf('%.1f-%.1f',cfg.phasefreq(f,1),cfg.phasefreq(f,2));
                    else
                      freq_str = sprintf('%.1f',cfg.phasefreq(f));
                    end
                    figure
                    surf(subSesEvData.(cell2mat(fn)).time,1:size(subSesEvData.(cell2mat(fn)).(phaseparam),1),squeeze(subSesEvData.(cell2mat(fn)).(phaseparam)(:,r,f,:)));
                    shading interp;view([0,90]);axis tight;
                    title(sprintf('%s, %s, %s, %s, %s Hz',exper.subjects{sub},strrep(exper.sessions{ses},'_',''),ana.eventValues{typ}{evVal},strrep(chan_str,'_',' '),freq_str));
                    xlabel('Time (s)');
                    ylabel('Trial number');
                    publishfig(gcf,0);
                    figfilename = sprintf('phase_%s_%s_%s_%s%sHz.png',exper.subjects{sub},exper.sessions{ses},ana.eventValues{typ}{evVal},chan_str,freq_str);
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
            
            %% normalize
            if ~isempty(cfg.normalize)
              fprintf('Doing a %s normalization.\n',cfg.normalize);
              if strcmp(cfg.normalize,'dB')
                if strcmp(cfg.baselinedata,'mod')
                  subSesEvData.(cell2mat(fn)).(param) = 10*log10(subSesEvData.(cell2mat(fn)).(param));
                elseif strcmp(cfg.baselinedata,'pow')
                  if strcmp(cfg.baselinetype,'relative') || strcmp(cfg.baselinetype,'relchange')
                    fprintf('Dividing entire trial by mean([%.2f %.2f] pre-stimulus) (baseline correcting) before normalizing.\n',cfg.baseline(1),cfg.baseline(2));
                    nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                    blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                    blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                    if strcmp(cfg.baselinetype,'relative')
                      subSesEvData.(cell2mat(fn)).(param) = 10*log10(subSesEvData.(cell2mat(fn)).(param) ./ blm);
                    elseif strcmp(cfg.baselinetype,'relchange')
                      subSesEvData.(cell2mat(fn)).(param) = 10*log10((subSesEvData.(cell2mat(fn)).(param) - blm) ./ blm);
                    end
                  elseif isempty(cfg.baselinetype)
                    subSesEvData.(cell2mat(fn)).(param) = 10*log10(subSesEvData.(cell2mat(fn)).(param));
                  else
                    % TODO: check if this could occur
                    fprintf('Fix this normalization!\n');
                    keyboard
                  end
                end
              elseif strcmp(cfg.normalize,'vec')
                for tr = 1:size(subSesEvData.(cell2mat(fn)).(param),1)
                  for ch = 1:size(subSesEvData.(cell2mat(fn)).(param),2)
                    % doesn't ignore NaNs, but is way faster
                    subSesEvData.(cell2mat(fn)).(param)(tr,ch,:,:) = normr(squeeze(subSesEvData.(cell2mat(fn)).(param)(tr,ch,:,:)));
                    % % try to ignore NaNs
                    %for fq = 1:size(subSesEvData.(cell2mat(fn)).(param),3)
                    %  thisTrl = squeeze(subSesEvData.(cell2mat(fn)).(param)(tr,ch,fq,:));
                    %  subSesEvData.(cell2mat(fn)).(param)(tr,ch,fq,~isnan(thisTrl)) = normc(thisTrl(~isnan(thisTrl)));
                    %end
                    % % trying another method, but it doesn't work
                    %thisTrl = squeeze(subSesEvData.(cell2mat(fn)).(param)(tr,ch,:,:));
                    %subSesEvData.(cell2mat(fn)).(param)(tr,ch,:,~isnan(thisTrl)) = normr(thisTrl(~isnan(thisTrl)));
                  end
                end
              else
                % run using feval (e.g., for log10 normalization)
                subSesEvData.(cell2mat(fn)).(param) = feval(cfg.normalize,subSesEvData.(cell2mat(fn)).(param));
                
                % option 2: normalize the evoked data
                if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
                  data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(powparam) = feval(cfg.normalize,data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(powparam));
                end
              end
            end
            
            % option 2: subtract evoked power from whole power; need to
            % do this after log normalization
            if strcmp(cfg.rmevoked,'yes') && strcmp(cfg.rmevokedpow,'yes')
              fprintf('Subtracting evoked power from whole power to get induced power.\n');
              evoked_data = repmat(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(powparam),[size(subSesEvData.(cell2mat(fn)).(powparam),1),1,1,1]);
              subSesEvData.(cell2mat(fn)).(powparam) = subSesEvData.(cell2mat(fn)).(powparam) - evoked_data;
            end
            
            %% baseline correct power
            if strcmp(cfg.baselinedata,'pow')
              if ~isempty(cfg.baselinetype) && ~strcmp(cfg.normalize,'dB')
                
                % get the baseline time indices
                blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                
                if strcmp(cfg.baselinetype,'zscore') && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(cell2mat(fn)).(param)) == 4
                  fprintf('Z-transforming data relative to mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline(1),cfg.baseline(2));
                  
                  nTrl = size(subSesEvData.(cell2mat(fn)).(param),1);
                  %nFrq = size(subSesEvData.(cell2mat(fn)).(param),3);
                  nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                  
                  blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                  blstd = repmat(nanmean(nanstd(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),0,1),4),[nTrl,1,1,nSmp]);
                  subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm) ./ blstd;
                  
                elseif (strcmp(cfg.baselinetype,'absolute') || strcmp(cfg.baselinetype,'relative') || strcmp(cfg.baselinetype,'relchange')) && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(cell2mat(fn)).(param)) == 4
                  %blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                  
                  nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                  
                  blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                  
                  if strcmp(cfg.baselinetype,'absolute')
                    fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial.\n',cfg.baseline(1),cfg.baseline(2));
                    subSesEvData.(cell2mat(fn)).(param) = subSesEvData.(cell2mat(fn)).(param) - blm;
                  elseif strcmp(cfg.baselinetype,'relative')
                    error('relative baseline is not working');
                    %fprintf('Dividing power by mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline(1),cfg.baseline(2));
                    %subSesEvData.(cell2mat(fn)).(param) = subSesEvData.(cell2mat(fn)).(param) ./ blm;
                  elseif strcmp(cfg.baselinetype,'relchange')
                    error('relchange baseline is not working');
                    %fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial and dividing by the mean.\n',cfg.baseline(1),cfg.baseline(2));
                    %subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm) ./ blm;
                  end
                  
                elseif strcmp(cfg.baselinetype,'absolute') && strcmp(cfg.output,'coh') && ndims(subSesEvData.(cell2mat(fn)).(param)) == 3
                  fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) coherence from entire trial.\n',cfg.baseline(1),cfg.baseline(2));
                  %blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                  
                  nSmp = size(subSesEvData.(cell2mat(fn)).(param),3);
                  
                  blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,blt),3),[1,1,nSmp]);
                  subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm);
                end
              end
            end
            
            %% debug
%             chan = 62;
%             figure;surf(subSesEvData.(cell2mat(fn)).time,subSesEvData.(cell2mat(fn)).freq,squeeze(mean(subSesEvData.(cell2mat(fn)).(powparam)(:,chan,:,:),1)));shading interp;view([0,90]);axis tight;colorbar;
%             title(sprintf('chan %d: rmevokedfourier: %s; rmevokedpow: %s',chan,cfg.rmevokedfourier,cfg.rmevokedpow));
%             %clim = [-1 1];
%             clim = [-3 3];
%             if strcmp(cfg.normalize,'dB')
%               clim = clim .* 10;
%             end
%             figure;imagesc(subSesEvData.(cell2mat(fn)).time,subSesEvData.(cell2mat(fn)).freq,squeeze(mean(subSesEvData.(cell2mat(fn)).(powparam)(:,chan,:,:),1)),clim);axis xy;colorbar
%             title(sprintf('chan %d: rmevokedfourier: %s; rmevokedpow: %s',chan,cfg.rmevokedfourier,cfg.rmevokedpow));
            
            %% save the data in a container struct
            fprintf('%s %s %s: ',exper.subjects{sub},sesStr,ana.eventValues{typ}{evVal});
            if (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && isfield(subSesEvData.(cell2mat(fn)),powparam) && ndims(subSesEvData.(cell2mat(fn)).(powparam)) == 4
              fprintf('Power data: ');
              %             elseif strcmp(cfg.keeptrials,'no') && (~strcmp(cfg.ftype,'pow') || ~strcmp(cfg.output,'pow')) && isfield(subSesEvData.(cell2mat(fn)),'powspctrm') && ndims(subSesEvData.(cell2mat(fn)).(powparam)) == 4
              %               error('\n%s %s %s: Can only keep trials for cfg.ftype=''pow''. You set it to ''%s''.\n',exper.subjects{sub},sesStr,ana.eventValues{typ}{evVal},cfg.ftype);
              %             elseif strcmp(cfg.keeptrials,'no') && ~isfield(subSesEvData.(cell2mat(fn)),'powspctrm')
              %               error('\n%s %s %s: Can only keep trials with ''powspctrm'' field. Please examine your data.\n',exper.subjects{sub},sesStr,ana.eventValues{typ}{evVal});
              %             elseif strcmp(cfg.keeptrials,'no') && isfield(subSesEvData.(cell2mat(fn)),'powspctrm') && ndims(subSesEvData.(cell2mat(fn)).(powparam)) ~= 4
              %               error('\n%s %s %s: Can only keep trials for ndims(powspctrm)==4. This data has ndims=%d.\n',exper.subjects{sub},sesStr,ana.eventValues{typ}{evVal},ndims(subSesEvData.(cell2mat(fn)).(powparam)));
              
              % set the configuration of ft_freqdescriptives
              cfg_fd = [];
              cfg_fd.keeptrials = cfg.keeptrials;
              
              if strcmp(cfg.keeptrials,'no')
                fprintf('Averaging over individual trials.');
              elseif strcmp(cfg.keeptrials,'yes')
                fprintf('Keeping individual trials.');
              end
              
%               if strcmp(cfg.equatetrials,'yes')
%                 fprintf(' Equating trial counts across all conditions (within a type).\n');
%                 fprintf('Randomly selecting %d of %d trials.\n',length(exper.randTrials{sub,ses,typ}(:,evVal)),size(subSesEvData.(cell2mat(fn)).(powparam),1));
%                 % get the subset of trials
%                 cfg_fd.trials = exper.randTrials{sub,ses,typ}(:,evVal);
%               elseif strcmp(cfg.equatetrials,'no')
%                 fprintf(' Using all trials.\n');
%                 % use all of the trials
%                 cfg_fd.trials = 'all';
%               end
              
              % TODO: do I need to update something in the FT structure now
              % that I get the subset of trials at the start?
              cfg_fd.trials = 'all';
              
                % use ft_freqdescriptives to average over individual
                % trials, if desired and the data is appropriate
              data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg_fd,subSesEvData.(cell2mat(fn)));
            elseif strcmp(cfg.output,'coh') && isfield(subSesEvData.(cell2mat(fn)),cohparam) && ndims(subSesEvData.(cell2mat(fn)).(cohparam)) == 3
              fprintf('Phase coherence data: Individual trials are lost.');
              if strcmp(cfg.equatetrials,'no')
                fprintf(' Using all trials. NB: this could affect the coherence calculation if trial counts are unequal!!\n');
              %elseif strcmp(cfg.equatetrials,'yes')
              %  fprintf(' Using randomly selected subsets of trials so each condition (within a type) has the same number of trials.\n');
              end
              data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = subSesEvData.(cell2mat(fn));
            elseif strcmp(cfg.output,'phase') && isfield(subSesEvData.(cell2mat(fn)),phaseparam) && ndims(subSesEvData.(cell2mat(fn)).(phaseparam)) == 4
              data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = subSesEvData.(cell2mat(fn));
            else
              fprintf('Need to figure out what to do for this case.\n');
              keyboard
            end
            
          else
            error('More than one field in data struct! There should only be one.\n');
          end
          
          % % old method
          % data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses) = load(inputfile);
          fprintf('Done.\n\n');
        else
          warning([mfilename,':noFileFound'],'NOT FOUND: %s\n',inputfile);
          fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n\n');
          data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.cfg.trl = [];
        end % if
      end % for evVal
    end % for typ
  end % for ses
end % for sub

end
