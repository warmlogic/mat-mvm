function [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(filename,replace_dataroot)
%MM_FT_LOADAD Load in the analysis details. Modify the file location if
%necessary (e.g., loading files on a different computer than the processing
%occurred).
%
% [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(filename,replace_dataroot)
%
% filename         = full path to analysisDetails.mat
% replace_dataroot = replace the dataroot stored in analysisDetails.mat
%                    with the one applicable now. true or false (default = false).
%                    OR
%                    {'old_dataroot','new_dataroot'}
%                    If using the cell option, dataroot should include
%                    everything before 'ft_data', something like:
%                    /path/to/ExpName/eeg/-1000_2000
% 

if nargin < 2
  replace_dataroot = false;
end

if ~islogical(replace_dataroot) && ~iscell(replace_dataroot)
  error('replace_dataroot must either be logical (true or false) or a cell containing {''old_dataroot'',''new_dataroot''}');
end

% make sure replace_dataroot is in the right format
do_the_replace = false;
if islogical(replace_dataroot)
  if replace_dataroot
    do_the_replace = true;
  end
elseif ~islogical(replace_dataroot)
  if length(replace_dataroot) ~= 2
    error('replace_dataroot must either be logical (true or false) or a cell containing {''old_dataroot'',''new_dataroot''}');
  elseif length(replace_dataroot) == 2
    do_the_replace = true;
  end
end

% load in analysisDetails.mat: exper, ana, dirs, files, cfg_proc
load(filename,'exper','ana','dirs','files','cfg_proc','cfg_pp');

% add in the electrode information, if necessary
if ~isfield(ana,'elec')
  ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);
end

% turn the session name into a string
if ~isfield(exper,'sesStr')
  % initialize for the sesStr
  exper.sesStr = cell(size(exper.sessions));
  
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
    % store the sesStr
    exper.sesStr{ses} = sesStr;
  end
end

%% process the dataroot string

if islogical(replace_dataroot) && replace_dataroot
  % get the old dataroot from the data we loaded in
  old_dataroot = dirs.dataroot;
  
  % make sure directories are set up correctly for all possible computers
  dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
  dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
  dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
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
    fprintf('Saved dirs.dataroot and dataroot for this location are the same (%s), nothing will be changed.\n',dirs.dataroot);
    do_the_replace = false;
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
  fprintf('Replacing saved dirs.dataroot (%s) with dataroot for this location (%s)...',old_dataroot,new_dataroot);
  
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
  
  % replace the old directories with the new
  dirs.saveDirProc = strrep(dirs.saveDirProc,old_dataroot,new_dataroot);
  dirs.saveDirRaw = strrep(dirs.saveDirRaw,old_dataroot,new_dataroot);
  dirs.saveDirFigs = strrep(dirs.saveDirFigs,old_dataroot,new_dataroot);
  % setting the new dataroot must go last
  dirs.dataroot = new_dataroot;
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
