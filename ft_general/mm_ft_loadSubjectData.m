function [data] = mm_ft_loadSubjectData(exper,dirs,eventValues,ftype)
%MM_FT_LOADSUBJECTDATA Load subject data into a full struct
%
% [data] = mm_ft_loadSubjectData(exper,dirs,ftype)
%
% ftype       = string included in the filename to load (e.g., 'tla' in
%               'data_tla_CR.mat')
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
      ses_str = exper.sessions{ses}{1};
      for i = 2:length(exper.sessions{ses})
        ses_str = cat(2,ses_str,'_',exper.sessions{ses}{i});
      end
    elseif ~iscell(exper.sessions{ses}) || (iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1)
      ses_str = exper.sessions{ses};
    end

    saveFileDir = fullfile(dirs.saveDirProc,exper.subjects{sub},ses_str);
    for typ = 1:length(eventValues)
      for evVal = 1:length(eventValues{typ})
        
        inputfile = fullfile(saveFileDir,sprintf('data_%s_%s.mat',ftype,eventValues{typ}{evVal}));
        if exist(inputfile,'file')
          fprintf('Loading %s...',inputfile);
          % load the data
          subSesEvData = load(inputfile);
          % get the name of the field
          fn = fieldnames(subSesEvData);
          % rename the field to 'data'
          if length(fn) == 1
            data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data = subSesEvData.(cell2mat(fn));
          else
            error('More than one field in data struct! There should only be one.');
          end
          
          % % old method
          % data.(eventValues{typ}{evVal}).sub(sub).ses(ses) = load(inputfile);
          fprintf('Done.\n');
        else
          fprintf('NOT FOUND: %s\n',inputfile);
          fprintf('Setting data.cfg.trl to empty brackets [] for compatibility.\n');
          data.(eventValues{typ}{evVal}).sub(sub).ses(ses).data.cfg.trl = [];
        end % if
      end % for evVal
    end % for typ
  end % for ses
end % for sub

end
