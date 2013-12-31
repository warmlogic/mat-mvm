function process_ft_data(ana,cfg_proc,exper,dirs)
%PROCESS_FT_DATA processes the raw FieldTrip data
%
% process_ft_data(ana,cfg_proc,exper,dirs)
%
%
% ana.ftFxn           = the FieldTrip function used to process the data
%                       (e.g., 'ft_timelockanalysis' or 'ft_freqanalysis')
% ana.ftype           = added to the filenames output from ana.ftFxn
%                       (e.g., 'tla' or 'pow')
% ana.usePeer         = use the peer toolbox (binary; default: 0); see
%                       http://www.ru.nl/neuroimaging/fieldtrip/ for more
%                       info on useing peer distributed computing
%
% ana.overwrite.proc  = to prevent overwriting of processed data. (binary;
%                       default: 1)
%                       DO NOT USE YET, NOT FULLY IMPLEMENTED. (TODO)
%
% This function will overwrite any raw or processed subject files without
% warning.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If an error occurs when running the ftFxn:
%
% If there's an error when processing an event, an error file will be saved
% (err_FTYPE_EVENT.mat) containing the MException object. Load this file
% for more information about what went wrong. It's likely that error is
% beacuse the raw file doesn't exist, resulting from either (1) the event
% values being incorrect (Net Station can be tricky with its event names)
% or (2) there were zero events after rejecting events due to artifacts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also: CREATE_FT_STRUCT, SEG2FT
%

%% make sure some required fields are set

% overwrite by default
if ~isfield(ana,'overwrite')
  ana.overwrite.proc = 1;
elseif isfield(ana,'overwrite') && ~isfield(ana.overwrite,'proc')
  ana.overwrite.proc = 1;
elseif isfield(ana,'overwrite') && isfield(ana.overwrite,'proc') && ana.overwrite.proc == 0
  % TODO
  error('Prevention of overwriting is not fully implemented. You must set overwrite to 1.');
end

% make sure exper.sessions is a cell
if ~iscell(exper.sessions)
  exper.sessions = {exper.sessions};
end

% need a FieldTrip analysis function
if ~isfield(ana,'ftFxn')
  error('Need to set the FieldTrip analysis function (e.g., ''ft_timelockanalysis'' or ''ft_freqanalysis'').');
end

% need a file type
if ~isfield(ana,'ftype')
  error('Need to set a file type in ana.ftype to be part of the outputfile name (e.g., ''tla'', ''pow'', ''powandcsd'', ''fourier'').');
elseif isfield(ana,'ftype') && strcmp(ana.ftype,'raw')
  error('Cannot set ana.ftype to ''raw'' because it is used for raw FieldTrip data.')
end

% don't run the peer toolbox by default
if ~isfield(ana,'usePeer')
  ana.usePeer = 0;
end

%% run the ft_*analysis function on the data, either using the peer toolbox or not

% make sure data_suffix is set up for loading inputfile
data_suffix = [];
if isfield(ana,'otherFxn')
  for i = 1:length(ana.otherFxn)
    cfg = ana.cfg_other{i};
    data_suffix = cat(2,data_suffix,sprintf('_%s',cfg.ftype));
  end
end

% if using the peer toolbox, set up the cfg first, then run
if ana.usePeer
  % initialize to store the cfgs
  cfg_peer = cell(1,length(exper.subjects)*length(exper.sessions)*length(exper.eventValues));
  % initialize a counter for indexing the cfgs
  count = 1;
  
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      
      % turn the session name into a string for easier printing
      if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
        sesStr = exper.sessions{ses}{1};
        for i = 2:length(exper.sessions{ses})
          sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
        end
      elseif ~iscell(exper.sessions{ses})
        sesStr = exper.sessions{ses};
      elseif iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1
        sesStr = cell2mat(exper.sessions{ses});
      end
      
      saveDirRawFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
      saveDirProcFile = fullfile(dirs.saveDirProc,exper.subjects{sub},sesStr);
      
      if ~exist(saveDirProcFile,'dir')
        mkdir(saveDirProcFile);
      end
      
      for evVal = 1:length(exper.eventValues)
        % get the name of this event type
        eventVal = exper.eventValues{evVal};
        
        inputfile = fullfile(saveDirRawFile,sprintf('data_raw%s_%s.mat',data_suffix,eventVal));
        outputfile = fullfile(saveDirProcFile,sprintf('data_%s_%s.mat',ana.ftype,eventVal));
        
        if (exist(outputfile,'file') && ana.overwrite.proc) || ~exist(outputfile,'file')
          if exist(inputfile,'file')
            cfg_peer{count} = cfg_proc;
            
            cfg_peer{count}.inputfile = inputfile;
            cfg_peer{count}.outputfile = outputfile;
            
            count = count + 1;
            
            % TODO: what happens when the inputfile doesn't exist because
            % of no trials or missing eventValues?
          end
        elseif (exist(outputfile,'file') && ~ana.overwrite.proc)
          fprintf('%s already exists and ana.overwrite.proc=0. Skipping processing of %s.\n',outputfile,eventVal);
        end
      end
    end
  end
  
  % run the analysis using the peer toolbox
  peercellfun(str2func(ana.ftFxn),cfg_peer);
  
else
  % run the analysis script once for each sub, ses, evVal
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      
      % turn the session name into a string for easier printing
      if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
        sesStr = exper.sessions{ses}{1};
        for i = 2:length(exper.sessions{ses})
          sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
        end
      elseif ~iscell(exper.sessions{ses})
        sesStr = exper.sessions{ses};
      elseif iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1
        sesStr = cell2mat(exper.sessions{ses});
      end
      
      saveDirRawFile = fullfile(dirs.saveDirRaw,exper.subjects{sub},sesStr);
      saveDirProcFile = fullfile(dirs.saveDirProc,exper.subjects{sub},sesStr);
      
      if ~exist(saveDirProcFile,'dir')
        mkdir(saveDirProcFile);
      end
      
      for evVal = 1:length(exper.eventValues)
        % get the name of this event type
        eventVal = exper.eventValues{evVal};
        
        inputfile = fullfile(saveDirRawFile,sprintf('data_raw%s_%s.mat',data_suffix,eventVal));
        outputfile = fullfile(saveDirProcFile,sprintf('data_%s_%s.mat',ana.ftype,eventVal));
        
        if (exist(outputfile,'file') && ana.overwrite.proc) || ~exist(outputfile,'file')
          if exist(inputfile,'file')
            
            cfg_proc.inputfile = inputfile;
            cfg_proc.outputfile = outputfile;
            
            fprintf('Running %s on %s, %s, %s...\n',ana.ftFxn,exper.subjects{sub},sesStr,eventVal);
            try
              feval(str2func(ana.ftFxn),cfg_proc);
              fprintf('Done with %s, %s, %s.\n',exper.subjects{sub},sesStr,eventVal);
            catch ME
              disp(ME)
              warning([mfilename,':ftFxn_problem'],'Something is wrong with %s, %s, %s.',exper.subjects{sub},sesStr,eventVal);
              errFile = fullfile(saveDirProcFile,sprintf('err_%s_%s.mat',ana.ftype,eventVal));
              fprintf('Saving error information in %s.\n',errFile);
              save(errFile,'ME');
            end % try
            
          end
        elseif (exist(outputfile,'file') && ~ana.overwrite.proc)
          fprintf('%s already exists and ana.overwrite.proc=0. Skipping processing of %s.\n',outputfile,eventVal);
        end
      end % for evVal
    end % ses
  end % sub
end % if ana.usePeer

end

