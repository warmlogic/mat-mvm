function [exper,ana,dirs,files,cfg_proc] = mm_ft_loadAD(filename,replace_dataroot)
%MM_FT_LOADAD Load in the analysis details. Modify the file location if
%necessary (e.g., loading files on a different computer than the processing
%occurred).
%
% [exper,ana,dirs,files,cfg_proc] = mm_ft_loadAD(filename,replace_dataroot)
%
% filename         = full path to analysisDetails.mat
% replace_dataroot = replace the dataroot stored in analysisDetails.mat
%                    with the one applicable now. 1 or 0 (default = 0).
%                    OR
%                    {'old_dataroot','new_dataroot'}
%                    If using the cell option, dataroot should include
%                    everything before 'ft_data', something like:
%                    /path/to/ExpName/eeg/-1000_2000
% 

if nargin < 2
  replace_dataroot = 0;
end

% make sure replace_dataroot is in the right format
if isnumeric(replace_dataroot)
  if replace_dataroot ~= 1 && replace_dataroot ~= 0
    error('replace_dataroot must either be logical (1 or 0) or a cell containing {''old_dataroot'',''new_dataroot''}');
  elseif replace_dataroot == 1
    do_the_replace = 1;
  elseif replace_dataroot == 0
    do_the_replace = 0;
  end
elseif ~isnumeric(replace_dataroot)
  if ~iscell(replace_dataroot)
    error('replace_dataroot must either be logical (1 or 0) or a cell containing {''old_dataroot'',''new_dataroot''}');
  elseif iscell(replace_dataroot)
    if length(replace_dataroot) ~= 2
      error('replace_dataroot must either be logical (1 or 0) or a cell containing {''old_dataroot'',''new_dataroot''}');
    elseif length(replace_dataroot) == 2
      do_the_replace = 1;
    end
  end
end

% load in analysisDetails.mat: exper, ana, dirs, files, cfg_proc
load(filename);

% add in the electrode information, if necessary
if ~isfield(ana,'elec')
  ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);
end

if isnumeric(replace_dataroot) && replace_dataroot == 1
  % get the old dataroot from the data we loaded in
  old_dataroot = dirs.dataroot;
  
  % make sure directories are set up correctly for all possible computers
  dirs.serverDir = fullfile('/Volumes','curranlab','Data');
  dirs.serverLocalDir = fullfile('/Volumes','RAID','curranlab','Data');
  dirs.dreamDir = fullfile('/data/projects/curranlab');
  dirs.localDir = fullfile(getenv('HOME'),'data');
  
  % pick the current dataroot
  if exist(dirs.serverDir,'dir')
    new_dataroot = dirs.serverDir;
  elseif exist(dirs.serverLocalDir,'dir')
    new_dataroot = dirs.serverLocalDir;
  elseif exist(dirs.dreamDir,'dir')
    new_dataroot = dirs.dreamDir;
  elseif exist(dirs.localDir,'dir')
    new_dataroot = dirs.localDir;
  else
    error('Data directory not found.');
  end
  
  if strcmp(dirs.dataroot,new_dataroot)
    fprintf('dirs.dataroot and the dataroot for this location are the same, nothing will be changed.\n');
    do_the_replace = 0;
  end
elseif iscell(replace_dataroot)
  old_dataroot = replace_dataroot{1};
  new_dataroot = replace_dataroot{2};
  
  % make sure the strings will replace each other correctly
  if strcmp(old_dataroot(end),filesep) && ~strcmp(new_dataroot(end),filesep)
    old_dataroot = old_dataroot(1:end-1);
  elseif ~strcmp(old_dataroot(end),filesep) && strcmp(new_dataroot(end),filesep)
    new_dataroot = new_dataroot(1:end-1);
  end
end

% do the replacement
if do_the_replace
  fprintf('Replacing %s with %s...',old_dataroot,new_dataroot);
  
  % replace the old directories with the new
  dirs.saveDir = strrep(dirs.saveDir,old_dataroot,new_dataroot);
  dirs.saveDirFigs = strrep(dirs.saveDirFigs,old_dataroot,new_dataroot);
  % setting the new dataroot must go last
  dirs.dataroot = new_dataroot;
  fprintf('Done.\n');
  
  if ~exist(dirs.dataroot,'dir')
    warning([mfilename,':newDatarootDoesNotExist'],'The new dataroot does not exist on this system. This may be OK if you are setting up files for another system.');
  end
end

end

% % dynamic renaming of cfg_*
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
