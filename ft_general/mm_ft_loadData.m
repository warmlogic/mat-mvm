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

fourierparam = 'fourierspctrm';
powparam = 'powspctrm';
%cohparam = 'cohspctrm';
% use 'powspctrm' for now because that's what FT expects
cohparam = 'powspctrm';
phaseparam = 'powspctrm';

if ~isfield(cfg,'keeptrials')
  cfg.keeptrials = 'no';
elseif isfield(cfg,'keeptrials')
  if ~strcmp(cfg.keeptrials,'yes') && ~strcmp(cfg.keeptrials,'no')
    error('cfg.keeptrials must be set to ''yes'' or ''no''.');
  end
end

if ~isfield(cfg,'equatetrials')
  cfg.equatetrials = 'no';
elseif isfield(cfg,'equatetrials')
  if ~strcmp(cfg.equatetrials,'yes') && ~strcmp(cfg.equatetrials,'no')
    error('cfg.equatetrials must be set to ''yes'' or ''no''.');
  end
end

if ~isfield(cfg,'rmevoked')
  cfg.rmevoked = 'no';
elseif isfield(cfg,'rmevoked')
  if ~strcmp(cfg.rmevoked,'yes') && ~strcmp(cfg.rmevoked,'no')
    error('cfg.rmevoked must be set to ''yes'' or ''no''.');
  end
end

if ~exist('data_evoked','var') || isempty(data_evoked)
  data_evoked = [];
  if strcmp(cfg.rmevoked,'yes')
    warning([mfilename,':rmevoked'],'Setting cfg.rmevoked=''no'' because data_evoked was not provided.\n');
    cfg.rmevoked = 'no';
  end
end

if ~isfield(cfg,'ftype')
  error('Must set cfg.ftype; used in the filename to load (e.g., ''tla'' in ''data_tla_CR.mat'');');
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
    if ~isfield(cfg,'plotphase')
      cfg.plotphase = 'yes';
    end
    if ~isfield(cfg,'phasefreq')
      cfg.phasefreq = [];
    end
    if ~isfield(cfg,'phaseroi')
      cfg.phaseroi = {{'E55'}};
    elseif ~iscell(cfg.phaseroi)
      cfg.phaseroi = {{cfg.phaseroi}};
    elseif ~iscell(cfg.phaseroi{1})
      cfg.phaseroi = {cfg.phaseroi};
    end
  end
else
  powflg = false;
  csdflg = false;
  fftflg = false;
  if ~isfield(cfg,'output')
    cfg.output = '???';
  end
end

if ~isfield(cfg,'normalize')
  cfg.normalize = '';
elseif ~isempty(cfg.normalize) && strcmp(cfg.output,'coh')
  error('Normalizing cfg.output=''coh'' is not allowed');
end

if ~isfield(cfg,'baselinetype')
  cfg.baselinetype = '';
end

if strcmp(cfg.normalize,'dB') && ~isempty(cfg.baselinetype) && ~strcmp(cfg.baselinetype,'relative')
  warning([mfilename,':dbNormBL'],'dB normalization uses a relative baseline type (division). setting cfg.baselinetype=''relative''.\n');
  cfg.baselinetype = 'relative';
end

if iscell(ana.eventValues)
  if ~iscell(ana.eventValues{1})
    ana.eventValues = {ana.eventValues};
  end
elseif ~iscell(ana.eventValues)
  ana.eventValues = {{ana.eventValues}};
end

% if we're equating trials, find out the number of trials to choose
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
          end
        end
      end
    end
  end
  for typ = 1:length(ana.eventValues)
    fprintf('Condition type consisting of: %s\n',sprintf(repmat('%s ',1,length(ana.eventValues{typ})),ana.eventValues{typ}{:}));
    if length(exper.subjects{1}) > 7
      tabchar_sub = '\t';
    else
      tabchar_sub = '';
    end
    % print out the trial counts for each subject
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

% process the data
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
            
            if powflg
              param = 'powspctrm';
              subSesEvData.(cell2mat(fn)).(param)(subSesEvData.(cell2mat(fn)).(param) == 0) = eps(0);
            end
            
            if csdflg
              param = 'crsspctrm';
            end
            
            % convert complex fourier-spectra to power or coherence
            if fftflg
              fprintf('Converting %s to %s.\n',cfg.ftype,cfg.output);
              if strcmp(cfg.output,'pow')
                if strcmp(cfg.rmevoked,'yes')
                  evoked_data = repmat(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(fourierparam),[size(subSesEvData.(cell2mat(fn)).(fourierparam),1),1,1,1]);
                  subSesEvData.(cell2mat(fn)).(fourierparam) = subSesEvData.(cell2mat(fn)).(fourierparam) - evoked_data;
                end
                subSesEvData.(cell2mat(fn)).(powparam) = abs(subSesEvData.(cell2mat(fn)).(fourierparam)).^2;
                subSesEvData.(cell2mat(fn)).(powparam)(subSesEvData.(cell2mat(fn)).(powparam) == 0) = eps(0);
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),fourierparam);
                param = powparam;
                % TODO? tell FT that we have power?
                % subSesEvData.(cell2mat(fn)) = ft_checkdata(subSesEvData.(cell2mat(fn)),'haspow','yes');
              elseif strcmp(cfg.output,'coh')
                % inter-trial phase coherence
                if strcmp(cfg.equatetrials,'yes')
                  % get the subset of trials
                  fprintf('Inter-trial coherence data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.randTrials{sub,ses,typ}(:,evVal)),size(subSesEvData.(cell2mat(fn)).(fourierparam),1));
                  subSesEvData.(cell2mat(fn)).(fourierparam) = subSesEvData.(cell2mat(fn)).(fourierparam)(exper.randTrials{sub,ses,typ}(:,evVal),:,:,:);
                end
                subSesEvData.(cell2mat(fn)).(cohparam) = subSesEvData.(cell2mat(fn)).(fourierparam) ./ abs(subSesEvData.(cell2mat(fn)).(fourierparam));
                subSesEvData.(cell2mat(fn)).(cohparam) = abs(squeeze(nanmean(subSesEvData.(cell2mat(fn)).(cohparam),1)));
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),fourierparam);
                param = cohparam;
                % change dimord
                % TODO: can I automate this (e.g., with ft_checkdata)?
                subSesEvData.(cell2mat(fn)).dimord = 'chan_freq_time';
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),'trialinfo');
              elseif strcmp(cfg.output,'phase')
                % single-trial phase
                if strcmp(cfg.equatetrials,'yes')
                  % get the subset of trials
                  fprintf('Phase data: Equating trial counts across conditions (within a type). Randomly selecting %d (of %d) trials.\n',length(exper.randTrials{sub,ses,typ}(:,evVal)),size(subSesEvData.(cell2mat(fn)).(fourierparam),1));
                  subSesEvData.(cell2mat(fn)).(fourierparam) = subSesEvData.(cell2mat(fn)).(fourierparam)(exper.randTrials{sub,ses,typ}(:,evVal),:,:,:);
                end
                
                if isempty(cfg.phasefreq)
                  cfg.phasefreq = subSesEvData.(cell2mat(fn)).freq;
                end
                
                if ~isvector(cfg.phasefreq) || (size(cfg.phasefreq,1) == 1 && size(cfg.phasefreq,2) == 2)
                  meanfourier = nan(size(subSesEvData.(cell2mat(fn)).(fourierparam),1),length(cfg.phaseroi),size(cfg.phasefreq,1),size(subSesEvData.(cell2mat(fn)).(fourierparam),4));
                elseif isvector(cfg.phasefreq) && size(cfg.phasefreq,2) == 1
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
                  elseif isvector(cfg.phasefreq) && size(cfg.phasefreq,2) == 1
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
                
              end
            end
            
            % normalize
            if ~isempty(cfg.normalize)
              fprintf('Doing a %s normalization.\n',cfg.normalize);
              if strcmp(cfg.normalize,'dB')
                if strcmp(cfg.baselinetype,'relative')
                  fprintf('Dividing entire trial by mean([%.2f %.2f] pre-stimulus) (baseline correcting) before normalizing.\n',cfg.baseline(1),cfg.baseline(2));
                  nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                  blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                  blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                  subSesEvData.(cell2mat(fn)).(param) = 10*log10(subSesEvData.(cell2mat(fn)).(param) ./ blm);
                elseif isempty(cfg.baselinetype)
                  subSesEvData.(cell2mat(fn)).(param) = 10*log10(subSesEvData.(cell2mat(fn)).(param));
                else
                  % TODO: check if this could occur
                  fprintf('Fix this normalization!\n');
                  keyboard
                end
              elseif strcmp(cfg.normalize,'vector')
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
                % run using feval
                subSesEvData.(cell2mat(fn)).(param) = feval(cfg.normalize,subSesEvData.(cell2mat(fn)).(param));
              end
            end
            
            % baseline
            if ~isempty(cfg.baselinetype) && ~strcmp(cfg.normalize,'dB')
              if strcmp(cfg.baselinetype,'zscore') && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(cell2mat(fn)).(param)) == 4
                fprintf('Z-transforming data relative to mean([%.2f %.2f] pre-stimulus).\n',cfg.baseline(1),cfg.baseline(2));
                blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                
                nTrl = size(subSesEvData.(cell2mat(fn)).(param),1);
                %nFrq = size(subSesEvData.(cell2mat(fn)).(param),3);
                nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                
                blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                blstd = repmat(nanmean(nanstd(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),0,1),4),[nTrl,1,1,nSmp]);
                subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm) ./ blstd;
                
              elseif strcmp(cfg.baselinetype,'absolute') && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(cell2mat(fn)).(param)) == 4
                fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) power from entire trial.\n',cfg.baseline(1),cfg.baseline(2));
                blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= baseline(2);
                
                nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                
                blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm);
                
              elseif strcmp(cfg.baselinetype,'absolute') && strcmp(cfg.output,'coh') && ndims(subSesEvData.(cell2mat(fn)).(param)) == 3
                fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) coherence from entire trial.\n',cfg.baseline(1),cfg.baseline(2));
                blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                
                nSmp = size(subSesEvData.(cell2mat(fn)).(param),3);
                
                blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,blt),3),[1,1,nSmp]);
                subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm);
              end
            end
            
            % save the data in a container struct
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
              
              if strcmp(cfg.equatetrials,'yes')
                fprintf(' Equating trial counts across all conditions (within a type).\n');
                fprintf('Randomly selecting %d of %d trials.\n',length(exper.randTrials{sub,ses,typ}(:,evVal)),size(subSesEvData.(cell2mat(fn)).(powparam),1));
                % get the subset of trials
                cfg_fd.trials = exper.randTrials{sub,ses,typ}(:,evVal);
              elseif strcmp(cfg.equatetrials,'no')
                fprintf(' Using all trials.\n');
                % use all of the trials
                cfg_fd.trials = 'all';
              end
              
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
