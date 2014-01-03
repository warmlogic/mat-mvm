function [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_loadAD(procDir,subjects,sesNames,eventValues,replaceDataroot)
%MM_LOADAD Load in the analysis details by aggregating subject details.
%Modify the file location if necessary (e.g., loading files on a different
%computer than the processing occurred).
%
% [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_loadAD(procDir,subjects,sesNames,eventValues,replaceDataroot)
%
% NB: uses the first subject as the basis for analysis details variables
% 
% procDir         = root path to where subject/session folders are that
%                   contain processed data and subjectDetails.mat
% subjects        = which subjects to load
% sesNames        = cell of strings for session names
% eventValues     = cell of cells (one per session) for which events to
%                   load
% replaceDataroot = replace the dataroot stored in analysisDetails.mat
%                   with the one applicable now. true or false (default = false).
%                   OR
%                   {'old_dataroot','new_dataroot'}
%
%                   If using the cell option, dataroot usually includes
%                   everything before 'ft_data', something like:
%                   /path/to/ExpName/eeg/-1000_2000
%
%                   but also may simply be an old string to be replaced by
%                   a new string!
% 

if ~exist('replaceDataroot','var') || isempty(replaceDataroot)
  replaceDataroot = false;
end

if ~islogical(replaceDataroot) && ~iscell(replaceDataroot)
  error('replace_dataroot must either be logical (true or false) or a cell containing {''old_dataroot'',''new_dataroot''}');
end

% make sure replace_dataroot is in the right format
do_the_replace = false;
if islogical(replaceDataroot)
  if replaceDataroot
    do_the_replace = true;
  end
elseif ~islogical(replaceDataroot)
  if length(replaceDataroot) ~= 2
    error('replace_dataroot must either be logical (true or false) or a cell containing {''old_dataroot'',''new_dataroot''}');
  elseif length(replaceDataroot) == 2
    do_the_replace = true;
  end
end

% use the first subject as the basis for analysis details variables
sub=1;
for ses = 1:length(sesNames)
  sdFile = fullfile(procDir,subjects{sub},sesNames{ses},'subjectDetails.mat');
  if exist(sdFile,'file')
    % load in analysisDetails.mat: exper, ana, dirs, files, cfg_proc
    sd = load(sdFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
    
    if ~exist('exper','var')
      % initialize using the first subject
      if isfield(sd,'exper')
        exper = sd.exper;
      else
        error('exper variable not found');
      end
      if isfield(sd,'ana')
        ana = sd.ana;
      else
        error('ana variable not found');
      end
      if isfield(sd,'dirs')
        dirs = sd.dirs;
      else
        error('dirs variable not found');
      end
      if isfield(sd,'files')
        files = sd.files;
      else
        error('files variable not found');
      end
      if isfield(sd,'cfg_proc')
        cfg_proc = sd.cfg_proc;
      else
        error('cfg_proc variable not found');
      end
      if isfield(sd,'cfg_pp')
        cfg_pp = sd.cfg_pp;
      else
        error('cfg_pp variable not found');
      end
      
      % only include the sessions and relevant info for the input sesNames
      keepSesInd = ismember(exper.sesNames,sesNames{ses});
      if ~any(keepSesInd)
        exper.sessions = exper.sessions(keepSesInd);
        exper.prepost = exper.prepost(keepSesInd);
        
        for k = 1:length(keepSesInd)
          if ~keepSesInd(k)
            exper.badChan = rmfield(exper.badChan,exper.sesNames{k});
            exper.nTrials = rmfield(exper.nTrials,exper.sesNames{k});
            exper.badEv = rmfield(exper.badEv,exper.sesNames{k});
          end
        end
      end
      
      % only include nTrials and badEv for the input eventValues
      for evVal = 1:length(exper.eventValues{ses})
        if ~ismember(exper.eventValues{ses}{evVal},eventValues{ses})
          exper.nTrials.(sesNames{ses}) = rmfield(exper.nTrials.(sesNames{ses}),exper.eventValues{ses}{evVal});
          exper.badEv.(sesNames{ses}) = rmfield(exper.badEv.(sesNames{ses}),exper.eventValues{ses}{evVal});
        end
      end
      
      % overwrite these fields using the input variables
      exper.subjects = subjects;
      exper.sesNames = sesNames;
      exper.eventValues = eventValues;
      
    else
      
      % concatenate: sessions
      exper.sessions = cat(2,exper.sessions,sd.exper.sessions);
      
      % add for each subject and session: badChan, badEv, nTrials
      exper.badChan.(exper.sesStr{ses}) = sd.exper.badChan.(exper.sesStr{ses})(sub);
      for evVal = 1:length(exper.eventValues{ses})
        exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = sd.exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal})(sub);
        exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = sd.exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal})(sub);
      end
      
    end
  else
    error('%s does not exist',sdFile);
  end
end

if length(subjects) > 1
  for sub = 2:length(subjects)
    for ses = 1:length(sesNames)
      sdFile = fullfile(procDir,subjects{sub},sesNames{ses},'subjectDetails.mat');
      if exist(sdFile,'file')
        
        % load in analysisDetails.mat: exper, ana, dirs, files, cfg_proc
        sd = load(sdFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
        
        % add for each subject and session: badChan, badEv, nTrials
        exper.badChan.(exper.sesStr{ses}) = cat(1,exper.badChan.(exper.sesStr{ses}),sd.exper.badChan.(exper.sesStr{ses}));
        for evVal = 1:length(exper.eventValues{ses})
          exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = cat(1,exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}),sd.exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}));
          exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = cat(1,exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}),sd.exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}));
        end
        
      else
        error('%s does not exist',sdFile);
      end
    end
  end
end

% add in the electrode information, if necessary
if ~isfield(ana,'elec')
  ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);
end

% % turn the session name into a string
% if ~isfield(exper,'sesStr')
%   % initialize for the sesStr
%   exper.sesStr = cell(size(exper.sessions));
%   
%   for ses = 1:length(exper.sessions)
%     
%     % turn the session name into a string for easier printing
%     if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
%       sesStr = exper.sessions{ses}{1};
%       for i = 2:length(exper.sessions{ses})
%         sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
%       end
%     elseif ~iscell(exper.sessions{ses})
%       sesStr = exper.sessions{ses};
%     elseif iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1
%       sesStr = cell2mat(exper.sessions{ses});
%     end
%     
%     % store the sesStr
%     exper.sesStr{ses} = sesStr;
%   end
% end

%% process the dataroot string

if islogical(replaceDataroot) && replaceDataroot
  % get the old dataroot from the data we loaded in
  old_dataroot = dirs.dataroot;
  
  % make sure directories are set up correctly for all possible computers
  dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
  dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
  dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
  dirs.localDir = fullfile(getenv('HOME'),'data');
  
  % pick the current dataroot
  if isfield(dirs,'serverDir') && exist(dirs.serverDir,'dir')
    new_dataroot = dirs.serverDir;
  elseif isfield(dirs,'serverLocalDir') && exist(dirs.serverLocalDir,'dir')
    new_dataroot = dirs.serverLocalDir;
  elseif isfield(dirs,'dreamDir') && exist(dirs.dreamDir,'dir')
    new_dataroot = dirs.dreamDir;
  elseif isfield(dirs,'localDir') && exist(dirs.localDir,'dir')
    new_dataroot = dirs.localDir;
  else
    error('Data directory not found.');
  end
  
  if strcmp(dirs.dataroot,new_dataroot)
    fprintf('Saved dirs.dataroot and dataroot for this location are the same (%s), nothing will be changed.\n',dirs.dataroot);
    do_the_replace = false;
  end
elseif iscell(replaceDataroot)
  old_dataroot = replaceDataroot{1};
  new_dataroot = replaceDataroot{2};
  
  if isempty(old_dataroot)
    error('First replacement string cannot be empty!')
  end
  
  % make sure the strings will replace each other correctly
  if strcmp(old_dataroot(end),filesep) && ~isempty(new_dataroot) && ~strcmp(new_dataroot(end),filesep)
    old_dataroot = old_dataroot(1:end-1);
  elseif ~strcmp(old_dataroot(end),filesep) && ~isempty(new_dataroot) && strcmp(new_dataroot(end),filesep)
    new_dataroot = new_dataroot(1:end-1);
  end
  
%   if strcmp(old_dataroot(end),filesep) && ~strcmp(new_dataroot(end),filesep)
%     old_dataroot = old_dataroot(1:end-1);
%   elseif ~strcmp(old_dataroot(end),filesep) && strcmp(new_dataroot(end),filesep)
%     new_dataroot = new_dataroot(1:end-1);
%   end
end

% do the replacement
if do_the_replace
  %fprintf('Replacing saved dirs.dataroot (%s) with dataroot for this location (%s)...',old_dataroot,new_dataroot);
  fprintf('Replacing strings in the dirs struct containing ''%s'' with ''%s''...',old_dataroot,new_dataroot);
  
  % backwards compatibility, we no longer make the "saveDir" field
  if ~isfield(dirs,'saveDirProc') && isfield(dirs,'saveDir')
    dirs.saveDirProc = dirs.saveDir;
    dirs = rmfield(dirs,'saveDir');
    % and put in a raw field for good measure
    if ~isfield(dirs,'saveDirRaw')
      seps = strfind(dirs.saveDirProc,filesep);
      dirs.saveDirRaw = strrep(dirs.saveDirProc,dirs.saveDirProc(seps(end)+1:end),'ft_raw');
    end
  end
  
  % simply replace the old string with the new string in all fields
  fn_dir = fieldnames(dirs);
  for i = 1:length(fn_dir)
    dirs.(fn_dir{i}) = strrep(dirs.(fn_dir{i}),old_dataroot,new_dataroot);
  end
  
%   % replace the old directories with the new
%   dirs.saveDirProc = strrep(dirs.saveDirProc,old_dataroot,new_dataroot);
%   dirs.saveDirRaw = strrep(dirs.saveDirRaw,old_dataroot,new_dataroot);
%   dirs.saveDirFigs = strrep(dirs.saveDirFigs,old_dataroot,new_dataroot);
%   % setting the new dataroot must go last
%   dirs.dataroot = new_dataroot;
  fprintf('Done.\n');
  
  if ~exist(dirs.dataroot,'dir')
    warning([mfilename,':newDatarootDoesNotExist'],'The new dataroot does not exist on this system. This may be OK, e.g. if you are setting up files for another system.');
  end
end

% %% dynamic renaming of cfg_*
% ad = load(filename);
% fn = fieldnames(ad);
% 
% % pull the struct names out into variable space; rename the cfg_* to cfg_proc
% for i = 1:length(fn)
%   if strfind(fn{i},'cfg_proc')
%     cfg_proc = eval(sprintf('ad.%s',fn{i}));
%   else
%     eval(sprintf('%s = ad.%s;',fn{i},fn{i}));
%   end
% end

%% set some default options for figures

% do we want to save them?
if ~isfield(files,'saveFigs')
  files.saveFigs = 1;
end

% set default figure file format to color encapsulated PS (vector graphics)
if ~isfield(files,'figPrintFormat')
  files.figPrintFormat = 'epsc2';
  %files.figPrintFormat = 'png';
else
  if strcmp(files.figPrintFormat(1:2),'-d')
    files.figPrintFormat = files.figPrintFormat(3:end);
  end
end

% default resolution for printed figures (DPI)
if ~isfield(files,'figPrintRes')
  files.figPrintRes = 150; % 150 is the Matlab default
end

% default font for figures
if ~isfield(files,'figFontName')
  files.figFontName = 'Helvetica'; % Helvetica is the Matlab default
end

end
