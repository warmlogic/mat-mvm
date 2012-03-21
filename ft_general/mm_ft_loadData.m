function [data] = mm_ft_loadData(cfg,exper,dirs,eventValues)
%MM_FT_LOADDATA Load subject data into a full struct
%
% [data] = mm_ft_loadData(cfg,exper,dirs,eventValues)
%
% cfg         = configuration; include normalization and baseline params
% exper       = the exper struct
% dirs        = the dirs struct, with the field saveDirProc (or saveDirRaw)
% eventValues = a cell of the event values to load (e.g., ana.eventValues)
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

powparam = 'powspctrm';
phaseparam = 'plvspctrm';

if ~isfield(cfg,'keeptrials')
  cfg.keeptrials = 'yes';
end

if ~isfield(cfg,'ftype')
  error('Must set cfg.ftype; used in the filename to load (e.g., ''tla'' in ''data_tla_CR.mat'');');
elseif strcmp(cfg.ftype,'pow')
  powflg = true;
  csdflg = false;
  fftflg = false;
  if ~isfield(cfg,'output')
    cfg.output = '';
  end
elseif strcmp(cfg.ftype,'powandcsd')
  powflg = false;
  csdflg = true;
  fftflg = false;
  if ~isfield(cfg,'output')
    cfg.output = '';
  end
elseif strcmp(cfg.ftype,'fourier')
  powflg = false;
  csdflg = false;
  fftflg = true;
  if ~isfield(cfg,'output')
    cfg.output = '';
  elseif strcmp(cfg.output,'phase') && strcmp(cfg.keeptrials,'yes')
    fprintf('Setting cfg.keeptrials=''no'' because phase is calcualted (averaged) across trials.\n');
    cfg.keeptrials = 'no';
  end
else
  powflg = false;
  csdflg = false;
  fftflg = false;
  if ~isfield(cfg,'output')
    cfg.output = '';
  end
end

if ~isfield(cfg,'normalize')
  cfg.normalize = '';
elseif ~isempty(cfg.normalize) && strcmp(cfg.output,'phase')
  error('Normalizing cfg.output=''phase'' is not allowed');
end

if ~isfield(cfg,'baselinetype')
  cfg.baselinetype = '';
end

if strcmp(cfg.normalize,'dB') && ~isempty(cfg.baselinetype) && ~strcmp(cfg.baselinetype,'relative')
  warning([mfilename,':dbNormBL'],'dB normalization uses a relative baseline type (division). setting cfg.baselinetype=''relative''.\n');
  cfg.baselinetype = 'relative';
end

if iscell(eventValues)
  if ~iscell(eventValues{1})
    eventValues = {eventValues};
  end
elseif ~iscell(eventValues)
  eventValues = {{eventValues}};
end

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
    
    for typ = 1:length(eventValues)
      for evVal = 1:length(eventValues{typ})
        
        inputfile = fullfile(saveFileDir,sprintf('data_%s_%s.mat',cfg.ftype,eventValues{typ}{evVal}));
        if exist(inputfile,'file')
          fprintf('Loading %s %s %s: %s...\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal},inputfile);
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
            
            % convert complex fourier-spectra to power or phase
            if fftflg
              param = 'fourierspctrm';
              fprintf('Converting %s to %s.\n',cfg.ftype,cfg.output);
              if strcmp(cfg.output,'pow')
                subSesEvData.(cell2mat(fn)).(powparam) = abs(subSesEvData.(cell2mat(fn)).(param)).^2;
                subSesEvData.(cell2mat(fn)).(powparam)(subSesEvData.(cell2mat(fn)).(powparam) == 0) = eps(0);
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),param);
                param = powparam;
                % TODO? tell FT that we have power?
                % subSesEvData.(cell2mat(fn)) = ft_checkdata(subSesEvData.(cell2mat(fn)),'haspow','yes');
              elseif strcmp(cfg.output,'phase')
                subSesEvData.(cell2mat(fn)).(phaseparam) = subSesEvData.(cell2mat(fn)).(param) ./ abs(subSesEvData.(cell2mat(fn)).(param));
                subSesEvData.(cell2mat(fn)).(phaseparam) = abs(squeeze(nanmean(subSesEvData.(cell2mat(fn)).(phaseparam),1)));
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),param);
                param = phaseparam;
                % change dimord
                % TODO: automate?
                subSesEvData.(cell2mat(fn)).dimord = 'chan_freq_time';
                subSesEvData.(cell2mat(fn)) = rmfield(subSesEvData.(cell2mat(fn)),'trialinfo');
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
                  error('Fix this normalization!');
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
                blstd = repmat(nanmean(std(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),0,1),4),[nTrl,1,1,nSmp]);
                subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm) ./ blstd;
                
              elseif strcmp(cfg.baselinetype,'absolute') && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && ndims(subSesEvData.(cell2mat(fn)).(param)) == 4
                fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) from entire trial.\n',cfg.baseline(1),cfg.baseline(2));
                blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= baseline(2);
                
                nSmp = size(subSesEvData.(cell2mat(fn)).(param),4);
                
                blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,:,blt),4),[1,1,1,nSmp]);
                subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm);
                
              elseif strcmp(cfg.baselinetype,'absolute') && strcmp(cfg.output,'phase') && ndims(subSesEvData.(cell2mat(fn)).(param)) == 3
                fprintf('Subtracting mean([%.2f %.2f] pre-stimulus) from entire trial.\n',cfg.baseline(1),cfg.baseline(2));
                blt = subSesEvData.(cell2mat(fn)).time >= cfg.baseline(1) & subSesEvData.(cell2mat(fn)).time <= cfg.baseline(2);
                
                nSmp = size(subSesEvData.(cell2mat(fn)).(param),3);
                
                blm = repmat(nanmean(subSesEvData.(cell2mat(fn)).(param)(:,:,blt),3),[1,1,nSmp]);
                subSesEvData.(cell2mat(fn)).(param) = (subSesEvData.(cell2mat(fn)).(param) - blm);
              end
            end
            
            % average
            if strcmp(cfg.keeptrials,'no') && (strcmp(cfg.ftype,'pow') || strcmp(cfg.output,'pow')) && isfield(subSesEvData.(cell2mat(fn)),'powspctrm') && ndims(subSesEvData.(cell2mat(fn)).(powparam)) == 4
              % use ft_freqdescriptives to average over individual trials
              % if desired and the data is appropriate
              fprintf('\n%s %s %s has individual trials. Using ft_freqdescriptives with cfg.keeptrials=''no'' to load only the average.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal});
              cfg_fd = [];
              cfg_fd.cfg.keeptrials = 'no';
              data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg_fd,subSesEvData.(cell2mat(fn)));
              %             elseif strcmp(cfg.keeptrials,'no') && (~strcmp(cfg.ftype,'pow') || ~strcmp(cfg.output,'pow')) && isfield(subSesEvData.(cell2mat(fn)),'powspctrm') && ndims(subSesEvData.(cell2mat(fn)).(powparam)) == 4
              %               error('\n%s %s %s: Can only keep trials for cfg.ftype=''pow''. You set it to ''%s''.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal},cfg.ftype);
              %             elseif strcmp(cfg.keeptrials,'no') && ~isfield(subSesEvData.(cell2mat(fn)),'powspctrm')
              %               error('\n%s %s %s: Can only keep trials with ''powspctrm'' field. Please examine your data.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal});
              %             elseif strcmp(cfg.keeptrials,'no') && isfield(subSesEvData.(cell2mat(fn)),'powspctrm') && ndims(subSesEvData.(cell2mat(fn)).(powparam)) ~= 4
              %               error('\n%s %s %s: Can only keep trials for ndims(powspctrm)==4. This data has ndims=%d.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal},ndims(subSesEvData.(cell2mat(fn)).(powparam)));
            else
              data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data = subSesEvData.(cell2mat(fn));
            end
            
          else
            error('More than one field in data struct! There should only be one.\n');
          end
          
          % % old method
          % data.(eventValues{typ}{evVal}).sub(sub).ses(ses) = load(inputfile);
          fprintf('Done.\n\n');
        else
          warning([mfilename,':noFileFound'],'NOT FOUND: %s\n',inputfile);
          fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n');
          data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data.cfg.trl = [];
        end % if
      end % for evVal
    end % for typ
  end % for ses
end % for sub

end
