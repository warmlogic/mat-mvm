function [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_loadAD(procDir,subjects,sesStr,replaceDataroot,replaceDatatype,subDir)
% MM_LOADAD Load in the analysis details by aggregating subject details.
% All event values are loaded!
%
% Can modify the file location if necessary (e.g., loading files on a
% different computer than the processing occurred).
%
% [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot)
%
% NB: uses the first subject as the basis for analysis details variables
% 
% procDir         = root path to where subject/session folders are that
%                   contain processed data and subjectDetails.mat
% subjects        = which subjects to load
% sesStr          = cell of strings for session names (same name(s) as
%                   directories where processed data was saved)
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
% replaceDatatype = similar to replaceDataroot but replaces the type of
%                   data in procDir (the final directory) with a different
%                   type; e.g., change tla to pow (tla needs to be in
%                   procDir because this is where subjectDetails.mat are
%                   stored) e.g., {'tla','pow'}
%
% subDir          = extra directory under ftpp
%
% TODO: save subjectDetails.mat in their own directory instead of in the
%       first data type to be processed
%

if ~exist('subDir','var') || isempty(subDir)
  subDir = '';
end
if ischar(subDir)
  dirs.subDir = subDir;
else
  error('subDir needs to be a character string');
end

if ~exist('replaceDatatype','var') || isempty(replaceDatatype)
  replaceDatatype = {};
end

if ischar(replaceDatatype)
  replaceDatatype = {replaceDatatype};
end

datatype_replace = false;
if ~isempty(replaceDatatype)
  datatype_replace = true;
end

if ~exist('replaceDataroot','var') || isempty(replaceDataroot)
  replaceDataroot = false;
end

if ~islogical(replaceDataroot) && ~iscell(replaceDataroot)
  error('replace_dataroot must either be logical (true or false) or a cell containing {''old_dataroot'',''new_dataroot''}');
end

% make sure replace_dataroot is in the right format
dataroot_replace = false;
if islogical(replaceDataroot)
  if replaceDataroot
    dataroot_replace = true;
  end
elseif ~islogical(replaceDataroot)
  if length(replaceDataroot) ~= 2
    error('replace_dataroot must either be logical (true or false) or a cell containing {''old_dataroot'',''new_dataroot''}');
  elseif length(replaceDataroot) == 2
    dataroot_replace = true;
  end
end

% use the first subject as the basis for analysis details variables
sub=1;
artTypes_all = {};
for ses = 1:length(sesStr)
  sdFile = fullfile(procDir,subjects{sub},sesStr{ses},'subjectDetails.mat');
  if exist(sdFile,'file')
    % load in analysisDetails.mat: exper, ana, dirs, files, cfg_proc
    sd = load(sdFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
    
    if ~exist('exper','var')
      % ses == 1, because we haven't set exper=sd.exper yet
      
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
      keepSesInd = ismember(exper.sesStr,sesStr{ses});
      if any(keepSesInd)
        % % % do not modify exper.sessions
        % % exper.sessions = exper.sessions(keepSesInd);
        exper.prepost = exper.prepost(keepSesInd);
        exper.sesStr = exper.sesStr(keepSesInd);
        exper.eventValues = exper.eventValues(keepSesInd);
        
        % remove extraneous session fields
        fieldsToCheck = {'badChan','nTrials','badEv','artifacts','trialinfo_allEv'};
        for ftc = 1:length(fieldsToCheck)
          if isfield(exper,fieldsToCheck{ftc})
            fn = fieldnames(exper.(fieldsToCheck{ftc}));
            exper.(fieldsToCheck{ftc}) = rmfield(exper.(fieldsToCheck{ftc}),fn(~keepSesInd));
          end
        end
        
        % collect all the artifact types
        if isfield(exper,'artifacts')
          kCount = 0;
          for k = 1:length(keepSesInd)
            if keepSesInd(k)
              kCount = kCount + 1;
              for evVal = 1:length(exper.eventValues{kCount})
                if ~isempty(exper.artifacts.(exper.sesStr{kCount}).(exper.eventValues{kCount}{evVal}))
                  artTypes = fieldnames(exper.artifacts.(exper.sesStr{kCount}).(exper.eventValues{kCount}{evVal}));
                  for at = 1:length(artTypes)
                    if ~ismember(artTypes{at},artTypes_all)
                      artTypes_all = cat(1,artTypes_all,artTypes{at});
                    end
                  end
                end
              end
            end
          end
        end
        
      else
        fprintf('sesStr does not contain any sessions that are in %s\n',sdFile);
        keyboard
      end
      
%       % only include nTrials and badEv for the input eventValues
%       for evVal = 1:length(exper.eventValues{ses})
%         if ~ismember(exper.eventValues{ses}{evVal},eventValues{ses})
%           exper.nTrials.(sesNames{ses}) = rmfield(exper.nTrials.(sesNames{ses}),exper.eventValues{ses}{evVal});
%           exper.badEv.(sesNames{ses}) = rmfield(exper.badEv.(sesNames{ses}),exper.eventValues{ses}{evVal});
%           if isfield(exper,'artifacts')
%             exper.artifacts.(sesNames{ses}) = rmfield(exper.artifacts.(sesNames{ses}),exper.eventValues{ses}{evVal});
%           end
%           if isfield(exper,'trialinfo_allEv')
%             exper.trialinfo_allEv.(sesNames{ses}) = rmfield(exper.trialinfo_allEv.(sesNames{ses}),exper.eventValues{ses}{evVal});
%           end
%         end
%       end
      
      % overwrite the subject field using the input variable
      exper.subjects = subjects;
      
      % % do not overwrite the eventValues field!! We need this for proper
      % % indexing of eventNumber
      % exper.eventValues = eventValues;
      
    else
      % ses > 1, because we already set exper=sd.exper
      
      % find out where in the sd struct this data is coming from
      thisSesInd = ismember(sd.exper.sesStr,sesStr{ses});
      
      % concatenate: sessions (length of sd.exper.sessions should be 1)
      exper.sessions = cat(2,exper.sessions,sd.exper.sessions);
      exper.sesStr = cat(2,exper.sesStr,sesStr{ses});
      exper.eventValues = cat(2,exper.eventValues,sd.exper.eventValues(thisSesInd));
      
      exper.prepost = cat(2,exper.prepost,sd.exper.prepost{thisSesInd});
      
      % add for each subject and session: badChan, badEv, nTrials
      exper.badChan.(sesStr{ses}) = sd.exper.badChan.(sesStr{ses})(sub);
      for evVal = 1:length(exper.eventValues{ses})
        exper.nTrials.(sesStr{ses}).(exper.eventValues{ses}{evVal}) = sd.exper.nTrials.(sesStr{ses}).(exper.eventValues{ses}{evVal})(sub);
        exper.badEv.(sesStr{ses}).(exper.eventValues{ses}{evVal}) = sd.exper.badEv.(sesStr{ses}).(exper.eventValues{ses}{evVal})(sub);
        if isfield(sd.exper,'artifacts')
          artTypes = fieldnames(sd.exper.artifacts.(sesStr{ses}).(exper.eventValues{ses}{evVal}));
          for at = 1:length(artTypes)
            exper.artifacts.(sesStr{ses}).(exper.eventValues{ses}{evVal}).(artTypes{at}) = sd.exper.artifacts.(sesStr{ses}).(exper.eventValues{ses}{evVal}).(artTypes{at})(sub);
            if ~ismember(artTypes{at},artTypes_all)
              artTypes_all = cat(1,artTypes_all,artTypes{at});
            end
          end
        end
        if isfield(sd.exper,'trialinfo_allEv')
          exper.trialinfo_allEv.(sesStr{ses}).(exper.eventValues{ses}{evVal}) = sd.exper.trialinfo_allEv.(sesStr{ses}).(exper.eventValues{ses}{evVal})(sub);
        end
      end
      
    end
  else
    error('%s does not exist',sdFile);
  end
end

if length(subjects) > 1
  for sub = 2:length(subjects)
    for ses = 1:length(sesStr)
      sdFile = fullfile(procDir,subjects{sub},sesStr{ses},'subjectDetails.mat');
      if exist(sdFile,'file')
        
        % load in analysisDetails.mat: exper, ana, dirs, files, cfg_proc
        sd = load(sdFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
        
        % add for each subject and session: badChan, badEv, nTrials
        exper.badChan.(exper.sesStr{ses}) = cat(1,exper.badChan.(exper.sesStr{ses}),sd.exper.badChan.(exper.sesStr{ses}));
        for evVal = 1:length(exper.eventValues{ses})
          exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = cat(1,exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}),sd.exper.nTrials.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}));
          exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = cat(1,exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}),sd.exper.badEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}));
          if isfield(sd.exper,'artifacts')
            artTypes = fieldnames(sd.exper.artifacts.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}));
            for at = 1:length(artTypes)
              if ~ismember(artTypes{at},artTypes_all)
                artTypes_all = cat(1,artTypes_all,artTypes{at});
              end
            end
          end
          if isfield(sd.exper,'trialinfo_allEv')
            exper.trialinfo_allEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}) = cat(1,exper.trialinfo_allEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}),sd.exper.trialinfo_allEv.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}));
          end
        end
        
      else
        error('%s does not exist',sdFile);
      end
    end
  end
end

% add all the artifact information
if isfield(exper,'artifacts')
  for at = 1:length(artTypes_all)
    for ses = 1:length(sesStr)
      for evVal = 1:length(exper.eventValues{ses})
        exper.artifacts.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}).(artTypes_all{at}) = cell(length(subjects),1);
      end
    end
  end
  for sub = 1:length(subjects)
    for ses = 1:length(sesStr)
      sdFile = fullfile(procDir,subjects{sub},sesStr{ses},'subjectDetails.mat');
      if exist(sdFile,'file')
        % load in analysisDetails.mat: exper, ana, dirs, files, cfg_proc
        sd = load(sdFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
        for evVal = 1:length(exper.eventValues{ses})
          if isfield(sd.exper,'artifacts')
            for at = 1:length(artTypes_all)
              if isfield(sd.exper.artifacts.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}),artTypes_all{at})
                exper.artifacts.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}).(artTypes_all{at}){sub} = sd.exper.artifacts.(exper.sesStr{ses}).(exper.eventValues{ses}{evVal}).(artTypes_all{at}){1};
              end
            end
          end
        end
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
    dataroot_replace = false;
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
if dataroot_replace
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

if datatype_replace
  if length(replaceDatatype) == 1
    newDatatype = replaceDatatype{1};
    
    % guessing at oldDatatype
    dirSep = strfind(dirs.saveDirProc,filesep);
    oldDatatype = dirs.saveDirProc(dirSep(end)+1:end);
    
    if ~strcmp(oldDatatype,newDatatype)
      fprintf('Guessing that old datatype is ''%s''... replacing with ''%s''.\n',oldDatatype,newDatatype);
      
      dirs.saveDirProc = dirs.saveDirProc(1:dirSep(end));
      dirs.saveDirProc = sprintf('%s%s',dirs.saveDirProc,newDatatype);
      
      dirs.saveDirFigs = sprintf('%s%sfigs',dirs.saveDirProc,filesep);
    else
      fprintf('Old datatype (''%s'') is the same as new datatype (''%s''). Not replacing anything.\n',oldDatatype,newDatatype);
    end
      
  elseif length(replaceDatatype) == 2
    oldDatatype = replaceDatatype{1};
    newDatatype = replaceDatatype{2};
    
    if ~strcmp(oldDatatype,newDatatype)
      fprintf('Replacing old datatype ''%s'' with new datatype ''%s''.\n',oldDatatype,newDatatype);
      
      dirSep = strfind(dirs.saveDirProc,filesep);
      if strcmp(oldDatatype,dirs.saveDirProc(dirSep(end)+1:end));
        dirs.saveDirProc = dirs.saveDirProc(1:dirSep(end));
        dirs.saveDirProc = sprintf('%s%s',dirs.saveDirProc,newDatatype);
        
        dirs.saveDirFigs = sprintf('%s%sfigs',dirs.saveDirProc,filesep);
      else
        error('old datatype ''%s'' does not match what is at the end of procDir (''%s'').',oldDatatype,dirs.saveDirProc(dirSep(end)+1:end));
      end
    else
      fprintf('Old datatype (''%s'') is the same as new datatype (''%s''). Not replacing anything.\n',oldDatatype,newDatatype);
    end
    
  else
    error('incorrect setup of replaceDatatype');
  end
end

if ~isempty(subDir)
  dirs.dataDir = strrep(dirs.dataDir,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirStem = strrep(dirs.saveDirStem,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirRaw = strrep(dirs.saveDirRaw,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirProc = strrep(dirs.saveDirProc,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirFigs = strrep(dirs.saveDirFigs,'ftpp',fullfile('ftpp',subDir));
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
