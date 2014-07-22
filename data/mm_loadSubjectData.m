function [data,exper] = mm_loadSubjectData(exper,dirs,ana,ftype,keeptrials,loadMethod,rmPreviousCfg)
%MM_LOADSUBJECTDATA Load subject data into a full struct
%
% NB: THIS FUNCTION WILL BE REPLACED WITH MM_FT_LOADDATA
%
% [data,exper] = mm_loadSubjectData(exper,dirs,ana,ftype,keeptrials,loadMethod,rmPreviousCfg)
%
% exper       = the exper struct
% dirs        = the dirs struct, with the field saveDirProc (or saveDirRaw)
% ana         = a struct containing a cell of the event values to load
%               (e.g., ana, with ana.eventValues)
% ftype       = string included in the filename to load (e.g., 'tla' in
%               'data_tla_CR.mat'); can be 'raw' to load the raw data.
% keeptrials  = optional; default=1. If loading powspctrm data created with
%               keeptrials='yes', set to 0 to return averaged data.
%               NB: uses ft_freqdescriptives
% loadMethod  = string defining how to load different event types; can be
%               'seg' or 'trialinfo' (segmented or FT's trialinfo); see
%               space_ft_seg_tla.m for an example of trialinfo.
% rmPreviousCfg = will remove data.cfg.previous with the intent to save
%                 memory
%

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

%warning([mfilename,':oldFxn'],'%s will soon become an old function! Use mm_ft_loadData_multiSes instead.',mfilename);

if ~exist('keeptrials','var') || isempty(keeptrials)
  warning('Keeping all trials!');
  keeptrials = true;
end

if ~exist('loadMethod','var') || isempty(loadMethod)
  error('Need to provide variable ''loadMethod''. See ''help %s'' for details.',mfilename);
  %loadMethod = 'seg';
end

if ~exist('rmPreviousCfg','var') || isempty(rmPreviousCfg)
  %warning('Upon loading of data, the data.cfg.previous field will NOT be removed! This may cause increased memory usage! See ''help %s'' for more info.',mfilename);
  %rmPreviousCfg = false;
  warning('Upon loading of data, the data.cfg.previous field WILL be removed to save memory! See ''help %s'' for more info.',mfilename);
  rmPreviousCfg = true;
end
if ~rmPreviousCfg
  warning('The data.cfg.previous field will NOT be removed! This may cause increased memory usage! See ''help %s'' for more info.',mfilename);
end

% eventValues = ana.eventValues;
% if isstruct(ana)
%   eventValues = ana.eventValues;
% elseif iscell(ana)
%   warning('The setup for %s has changed, please input ana instead of eventValues!',mfilename);
%   eventValues = ana;
% end
%
% if iscell(eventValues)
%   if ~iscell(eventValues{1})
%     eventValues = {eventValues};
%   end
% elseif ~iscell(eventValues)
%   eventValues = {{eventValues}};
% end

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    %     % turn the session name into a string for easier printing
    %     if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
    %       sesStr = exper.sessions{ses}{1};
    %       for i = 2:length(exper.sessions{ses})
    %         sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
    %       end
    %     elseif ~iscell(exper.sessions{ses}) || (iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1)
    %       sesStr = exper.sessions{ses};
    %     end
    
    sesStr = exper.sesStr{ses};
    
    if strcmp(ftype,'raw')
      saveFileDir = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
    else
      saveFileDir = fullfile(dirs.saveDirProc,exper.subjects{sub},sesStr);
    end
    for evVal = 1:length(ana.eventValues{ses})
      
      inputfile = fullfile(saveFileDir,sprintf('data_%s_%s.mat',ftype,ana.eventValues{ses}{evVal}));
      if exist(inputfile,'file')
        fprintf('Loading %s %s %s: %s...\n',exper.subjects{sub},sesStr,ana.eventValues{ses}{evVal},inputfile);
        % load the data
        subSesEvData = load(inputfile);
        % get the name of the field
        fn = fieldnames(subSesEvData);
        % rename the field to 'data'
        if length(fn) == 1
          data_fn = cell2mat(fn);
        else
          error('More than one field in data struct! There should only be one.\n');
        end
        
        % save space by removing a potentially large 'previous' field
        if rmPreviousCfg
          if isfield(subSesEvData.(data_fn).cfg,'previous')
            fprintf('Removing ''previous'' field from data.cfg\n');
            subSesEvData.(data_fn).cfg = rmfield(subSesEvData.(data_fn).cfg,'previous');
          end
        end
        
        if strcmp(ftype,'pow')
          warning('You probably want to use mm_ft_loadData_multiSes.m instead of %s for loading time-freq data.',mfilename);
          
          if ~keeptrials
            if isfield(subSesEvData.(data_fn),'powspctrm') && ndims(subSesEvData.(data_fn).powspctrm) == 4
              if strcmp(loadMethod,'seg')
                % use ft_freqdescriptives to average over individual trials
                % if desired and the data is appropriate
                fprintf('\n%s %s %s has individual trials. Using ft_freqdescriptives with keeptrials=''no'' to load only the average.\n',exper.subjects{sub},sesStr,ana.eventValues{ses}{evVal});
                cfg_fd = [];
                %cfg_fd.keeptrials = 'no';
                %data.(sesStr).(ana.eventValues{ses}{evVal}).sub(sub).data = ft_freqdescriptives(cfg_fd,subSesEvData.(data_fn));
                cfg_fd.avgoverrpt = 'yes';
                data.(sesStr).(ana.eventValues{ses}{evVal}).sub(sub).data = ft_selectdata(cfg_fd,subSesEvData.(data_fn));
                
              elseif strcmp(loadMethod,'trialinfo')
                trl_order = ana.trl_order.(ana.eventValues{ses}{evVal});
                
                for es = 1:length(ana.eventValuesSplit{ses}{evVal})
                  eventValue = ana.eventValuesSplit{ses}{evVal}{es};
                  
                  fprintf('Selecting %s trials...\n',eventValue);
                  
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
                  
                  cfg_fd = [];
                  cfg_fd.trials = eval(expr);
                  %cfg_fd.keeptrials = 'no';
                  %data.(sesStr).(eventValue).sub(sub).data = ft_freqdescriptives(cfg_fd,subSesEvData.(data_fn));
                  cfg_fd.avgoverrpt = 'yes';
                  data.(sesStr).(eventValue).sub(sub).data = ft_selectdata(cfg_fd,subSesEvData.(data_fn));
                  
                  % put in the trial counts
                  exper.nTrials.(sesStr).(eventValue)(sub,1) = sum(cfg_fd.trials);
                  
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
                end % es
              end
              %elseif ~strcmp(ftype,'pow')
              %  error('\n%s %s %s: Can only keep trials for ftype=''pow''. You set it to ''%s''.\n',exper.subjects{sub},sesStr,ana.eventValues{ses}{evVal},ftype);
            elseif ~isfield(subSesEvData.(data_fn),'powspctrm')
              error('\n%s %s %s: Need ''powspctrm'' field. Please examine your data.\n',exper.subjects{sub},sesStr,ana.eventValues{ses}{evVal});
            elseif isfield(subSesEvData.(data_fn),'powspctrm') && ndims(subSesEvData.(data_fn).powspctrm) ~= 4
              error('\n%s %s %s: Needed ndims(powspctrm)==4. This data has ndims=%d.\n',exper.subjects{sub},sesStr,ana.eventValues{ses}{evVal},ndims(subSesEvData.(data_fn).powspctrm));
            end
          elseif keeptrials
            error('Use mm_ft_loadData_multiSes.m instead of %s for keeping time-freq trials.',mfilename);
          end
        elseif strcmp(ftype,'tla')
          if strcmp(loadMethod,'seg')
            if ~keeptrials
              error('need to keep trials')
            end
            data.(sesStr).(ana.eventValues{ses}{evVal}).sub(sub).data = subSesEvData.(data_fn);
            
          elseif strcmp(loadMethod,'trialinfo')
            trl_order = ana.trl_order.(ana.eventValues{ses}{evVal});
            
            for es = 1:length(ana.eventValuesSplit{ses}{evVal})
              eventValue = ana.eventValuesSplit{ses}{evVal}{es};
              
              fprintf('Selecting %s trials...\n',eventValue);
              
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
              
              cfg = [];
              cfg.trials = eval(expr);
              if keeptrials
                cfg.keeptrials = 'yes';
              else
                cfg.keeptrials = 'no';
              end
              if any(cfg.trials == 1)
                data.(sesStr).(eventValue).sub(sub).data = ft_timelockanalysis(cfg, subSesEvData.(data_fn));
                %data.(sesStr).(eventValue).sub(sub).data = ft_redefinetrial(cfg, subSesEvData.(data_fn));
              else
                warning('No events found for %s %s %s',exper.subjects{sub},sesStr,eventValue);
                fprintf('\tEvaluated this expression: %s\n',expr);
                data.(sesStr).(eventValue).sub(sub).data.trial = [];
              end
              
              % put in the trial counts
              exper.nTrials.(sesStr).(eventValue)(sub,1) = sum(cfg.trials);
              
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
            end % es
            
          end % loadMethod
        end % ftype
        
        % % old method
        % data.(ana.eventValues{ses}{evVal}).sub(sub).ses(ses) = load(inputfile);
        fprintf('Done.\n\n');
      else
        warning([mfilename,':noFileFound'],'NOT FOUND: %s\n',inputfile);
        fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n');
        if strcmp(loadMethod,'seg')
          data.(sesStr).(ana.eventValues{ses}{evVal}).sub(sub).data.cfg.trl = [];
        elseif strcmp(loadMethod,'trialinfo')
          for es = 1:length(ana.eventValuesSplit{ses}{evVal})
            eventValue = ana.eventValuesSplit{ses}{evVal}{es};
            data.(sesStr).(eventValue).sub(sub).data.cfg.trl = [];
          end
        end
      end % if
    end % for evVal
  end % for ses
end % for sub

end % function
