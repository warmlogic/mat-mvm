function [dirs] = mm_checkDirs(dirs,replaceDatatype,subDir)

% function [dirs] = mm_checkDirs(dirs,replaceDatatype,subDir)
%
% Changes dataroot to match the local setup. Can also change datatype.
%
% dirs = dirs struct
%
% replaceDatatype = replaces the type of data in procDir (the final
%                   directory) with a different type; e.g., change tla to
%                   pow (tla needs to be in procDir because this is where
%                   subjectDetails.mat are stored) e.g., {'tla','pow'}. Set
%                   to {} if not using.
%
% subDir          = character string. an extra directory under ftpp.
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
  %dataroot_replace = false;
  return
else
  dataroot_replace = true;
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

%% datatype

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

%% check on subDir
if ~isempty(subDir)
  dirs.dataDir = strrep(dirs.dataDir,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirStem = strrep(dirs.saveDirStem,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirRaw = strrep(dirs.saveDirRaw,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirProc = strrep(dirs.saveDirProc,'ftpp',fullfile('ftpp',subDir));
  dirs.saveDirFigs = strrep(dirs.saveDirFigs,'ftpp',fullfile('ftpp',subDir));
end
