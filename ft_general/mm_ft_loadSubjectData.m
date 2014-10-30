function [data,exper] = mm_ft_loadSubjectData(exper,dirs,ana,ftype,keeptrials,loadMethod)
%MM_FT_LOADSUBJECTDATA Load subject data into a full struct
%
% NB: THIS FUNCTION IS OUTDATED. USE MM_LOADSUBJECTDATA INSTEAD.
%
% [data,exper] = mm_ft_loadSubjectData(exper,dirs,ana,ftype,keeptrials,loadMethod)
%
% exper       = the exper struct
% dirs        = the dirs struct, with the field saveDirProc (or saveDirRaw)
% ana         = a struct containing a cell of the event values to load (e.g., ana, with ana.eventValues)
% ftype       = string included in the filename to load (e.g., 'tla' in
%               'data_tla_CR.mat'); can be 'raw' to load the raw data.
% keeptrials  = optional; default=1. If loading powspctrm data created with
%               keeptrials='yes', set to 0 to return averaged data.
%               NB: uses ft_freqdescriptives
% loadMethod  = string defining how to load different event types; can be
%               'seg' or 'trialinfo' (segmented or FT's trialinfo); see
%               space_ft_seg_tla.m for an example of trialinfo.

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

% error('use mm_loadSubjectData.m instead of this function (%s).',mfilename);

%warning([mfilename,':oldFxn'],'%s will soon become an old function! Use mm_ft_loadData instead.',mfilename);

if ~exist('keeptrials','var') || isempty(keeptrials)
  warning('Keeping all trials!');
  keeptrials = 1;
end

  loadMethod = 'seg';
if ~exist('loadMethod','var') || isempty(loadMethod)
  error('Need to provide variable ''loadMethod''. See ''help %s'' for details.',mfilename);
end

if isstruct(ana)
  eventValues = ana.eventValues;
elseif iscell(ana)
  warning('The setup for %s has changed, please input ana instead of eventValues!',mfilename);
  eventValues = ana;
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
    
    if strcmp(ftype,'raw')
      saveFileDir = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
    else
      saveFileDir = fullfile(dirs.saveDirProc,exper.subjects{sub},sesStr);
    end
    for typ = 1:length(eventValues)
      for evVal = 1:length(eventValues{typ})
        
        inputfile = fullfile(saveFileDir,sprintf('data_%s_%s.mat',ftype,eventValues{typ}{evVal}));
        if exist(inputfile,'file')
          fprintf('Loading %s %s %s: %s...\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal},inputfile);
          % load the data
          subSesEvData = load(inputfile);
          % get the name of the field
          fn = fieldnames(subSesEvData);
          % rename the field to 'data'
          if length(fn) == 1
            data_fn = cell2mat(fn);
            
            if keeptrials == 0 && strcmp(ftype,'pow') && isfield(subSesEvData.(data_fn),'powspctrm') && ndims(subSesEvData.(data_fn).powspctrm) == 4
              if strcmp(loadMethod,'seg')
                % use ft_freqdescriptives to average over individual trials
                % if desired and the data is appropriate
                fprintf('\n%s %s %s has individual trials. Using ft_freqdescriptives with keeptrials=''no'' to load only the average.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal});
                cfg_fd = [];
                cfg_fd.keeptrials = 'no';
                data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg_fd,subSesEvData.(data_fn));
                
              elseif strcmp(loadMethod,'trialinfo')
                trl_order = ana.trl_order.(eventValues{typ}{evVal});
                
                for es = 1:length(ana.eventValuesSplit{typ})
                  fprintf('Selecting %s trials...\n',ana.eventValuesSplit{typ}{es});
                  
                  expr = ana.trl_expr{typ}{es};
                  
                  for to = 1:length(trl_order)
                    % replace the field name in the logical expression
                    % ana.trl_expr{typ}{es} with the column number found in
                    % ana.trl_order.(eventValues{typ}{evVal})
                    
                    % find the column number
                    fieldNum = find(ismember(trl_order,trl_order{to}));
                    
                    % set up the string to be replaced
                    r_exp = ['\<' trl_order{to} '\>'];
                    % set up the replacement string
                    r_str = sprintf('subSesEvData.%s.trialinfo(:,%d)',data_fn, fieldNum);
                    
                    expr = regexprep(expr,r_exp,r_str);
                  end
                  
                  cfg_fd = [];
                  cfg_fd.trials = eval(expr);
                  cfg_fd.keeptrials = 'no';
                  data.(ana.eventValuesSplit{typ}{es}).sub(sub).ses(ses).data = ft_freqdescriptives(cfg_fd,subSesEvData.(data_fn));
                  % put in the trial counts
                  exper.nTrials.(ana.eventValuesSplit{typ}{es})(sub,ses) = sum(cfg.trials);
                end % es
              end
              
            elseif keeptrials == 0 && ~strcmp(ftype,'pow')
              error('\n%s %s %s: Can only keep trials for ftype=''pow''. You set it to ''%s''.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal},ftype);
            elseif keeptrials == 0 && ~isfield(subSesEvData.(data_fn),'powspctrm')
              error('\n%s %s %s: Can only keep trials with ''powspctrm'' field. Please examine your data.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal});
            elseif keeptrials == 0 && isfield(subSesEvData.(data_fn),'powspctrm') && ndims(subSesEvData.(data_fn).powspctrm) ~= 4
              error('\n%s %s %s: Can only keep trials for ndims(powspctrm)==4. This data has ndims=%d.\n',exper.subjects{sub},sesStr,eventValues{typ}{evVal},ndims(subSesEvData.(data_fn).powspctrm));
            else
              if strcmp(loadMethod,'seg')
                data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data = subSesEvData.(data_fn);
                
              elseif strcmp(loadMethod,'trialinfo')
                trl_order = ana.trl_order.(eventValues{typ}{evVal});
                
                for es = 1:length(ana.eventValuesSplit{typ})
                  fprintf('Selecting %s trials...\n',ana.eventValuesSplit{typ}{es});
                  
                  expr = ana.trl_expr{typ}{es};
                  
                  for to = 1:length(trl_order)
                    % replace the field name in the logical expression
                    % ana.trl_expr{typ}{es} with the column number found in
                    % ana.trl_order.(eventValues{typ}{evVal})
                    
                    % find the column number
                    fieldNum = find(ismember(trl_order,trl_order{to}));
                    
                    % set up the string to be replaced
                    r_exp = ['\<' trl_order{to} '\>'];
                    % set up the replacement string
                    r_str = sprintf('subSesEvData.%s.trialinfo(:,%d)',data_fn, fieldNum);
                    
                    expr = regexprep(expr,r_exp,r_str);
                  end
                  
                  cfg = [];
                  cfg.trials = eval(expr);
                  cfg.keeptrials = 'yes';
                  data.(ana.eventValuesSplit{typ}{es}).sub(sub).ses(ses).data = ft_timelockanalysis(cfg, subSesEvData.(data_fn));
                  %data.(ana.eventValuesSplit{typ}{es}).sub(sub).ses(ses).data = ft_redefinetrial(cfg, subSesEvData.(data_fn));
                  
                  % put in the trial counts
                  exper.nTrials.(ana.eventValuesSplit{typ}{es})(sub,ses) = sum(cfg.trials);
                end % es
                
              end % loadMethod
            end % keeptrials
          else
            error('More than one field in data struct! There should only be one.\n');
          end
          
          % % old method
          % data.(eventValues{typ}{evVal}).sub(sub).ses(ses) = load(inputfile);
          fprintf('Done.\n\n');
        else
          warning([mfilename,':noFileFound'],'NOT FOUND: %s\n',inputfile);
          fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n');
          if strcmp(loadMethod,'seg')
            data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data.cfg.trl = [];
          elseif strcmp(loadMethod,'trialinfo')
            for es = 1:length(ana.eventValuesSplit{typ})
              data.(ana.eventValuesSplit{typ}{es}).sub(sub).ses(ses).data.cfg.trl = [];
            end
          end
        end % if
      end % for evVal
    end % for typ
  end % for ses
end % for sub

end % function
