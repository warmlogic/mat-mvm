%% initialize

myMatlabDir = fullfile(getenv('HOME'),'Documents','MATLAB');

%% set up eegToolbox path
eegToolboxDir = dir(fullfile(myMatlabDir,'eegToolbox*'));
if ~isempty(eegToolboxDir)
  eegToolboxDir = fullfile(myMatlabDir,eegToolboxDir.name);
  % add top folder and all subfolders
  addpath(genpath(eegToolboxDir));
end

%% set up eeglab path
eeglabDir = dir(fullfile(myMatlabDir,'eeglab*'));
if ~isempty(eeglabDir)
  eeglabDir = fullfile(myMatlabDir,eeglabDir.name);
  % add top folder and all subfolders
  addpath(genpath(eeglabDir));

  % remove eeglab's external directory if it was added
  eeglabExtDir = fullfile(eeglabDir,'external');
  if ~isempty(eeglabExtDir)
    rmpath(genpath(eeglabExtDir));
  end
  
  % % remove eeglab's fieldtrip directory if it was added
  % eeglabFtDir = dir(fullfile(eeglabDir,'external','fieldtrip*'));
  % if ~isempty(eeglabFtDir)
  %   eeglabFtDir = fullfile(myMatlabDir,eeglabFtDir.name);
  %   rmpath(genpath(eeglabFtDir));
  % end
end

%% set up MVPA toolbox path
mvpaDir = fullfile(myMatlabDir,'mvpa');
if exist(mvpaDir,'dir')
  % add only the top folder
  addpath(mvpaDir);
  
  % add the subdirectories that MVPA needs
  mvpa_add_paths;
end

%% set up fieldtrip path
ftDir = dir(fullfile(myMatlabDir,'fieldtrip-*'));
if ~isempty(ftDir)
  ftDir = fullfile(myMatlabDir,ftDir.name);
  % add only the top folder
  addpath(ftDir);
  % add the subdirectories that FT needs
  ft_defaults;
  
  % add the peer directory
  addpath(fullfile(ftDir,'peer'));
  
  % add the multivariate directory (outdated; use DMLT)
  %addpath(genpath(fullfile(ftDir,'multivariate')));
  
  % add the SPM directory
  %addpath(fullfile(ftDir,'external','spm8'));
  
  % % remove fieldtrip's external directory
  % ftExtDir = fullfile(ftDir,'external');
  % if ~isempty(ftExtDir)
  %   rmpath(genpath(ftExtDir));
  % end
  
  % remove fieldtrip's eeglab directory if it was added
  ftEeglabDir = dir(fullfile(ftDir,'external','eeglab*'));
  if ~isempty(ftEeglabDir)
    ftEeglabDir = fullfile(myMatlabDir,ftEeglabDir.name);
    rmpath(genpath(ftEeglabDir));
  end
  
  % functions that get in the way
  conflictFiles = {{'progress','MVPA'}}; % {{fxn, other package},{fxn, other package}}
  for i = 1:length(conflictFiles)
    if exist(fullfile(ftDir,'utilities','compat',[conflictFiles{i}{1},'.m']),'file')
      fprintf('Found an old/unnecessary FieldTrip function %s.m that conflicts with %s''s function. Moving FT''s to %s_old.m.\n',conflictFiles{i}{1},conflictFiles{i}{2},conflictFiles{i}{1});
      unix(sprintf('mv %s %s',fullfile(ftDir,'utilities','compat',[conflictFiles{i}{1},'.m']),fullfile(ftDir,'utilities','compat',[conflictFiles{i}{1},'_old.m'])));
    end
  end
  
end

%% setup DMLT
dmltDir = fullfile(myMatlabDir,'DMLT');
if exist(dmltDir,'dir')
  % add top folder and all subfolders
  addpath(genpath(dmltDir));
end

% if you get an error because libgfortran is not installed, follow the
% instructions from here <http://www.cs.ubc.ca/~hoffmanm/matlab.html> to
% download the libgfortran 4.2 tarball from here
% <http://r.research.att.com/tools/> listed under OS X 10.6 and copy
% libgfortran.* files to /usr/local/lib.

%% set up EP_Toolkit path
epDir = dir(fullfile(myMatlabDir,'EP_Toolkit*'));
if ~isempty(epDir)
  epDir = fullfile(myMatlabDir,epDir.name);
  % add top folder and all subfolders
  addpath(genpath(epDir));
end

%% add my other analysis scripts/files
addpath(genpath(fullfile(myMatlabDir,'recogmodel_mvm')));

%% add my experiment, fieldtrip, and RM ANOVA scripts
addpath(genpath(fullfile(myMatlabDir,'mat-mvm')));

%% put the ~/Documents/MATLAB folder at the top of the path
addpath(myMatlabDir);

%% remove CVS and .svn directories from path
entries = regexp(path, ['[^',pathsep,']*',pathsep], 'match');
for i = 1:length(entries)
  entry = char(entries{i});
  if ~isempty(strfind(entry, '.git'))
    rmpath(entry);
  end
  if ~isempty(strfind(entry, '.svn'))
    rmpath(entry);
  end
  if ~isempty(strfind(entry, 'CVS'))
    rmpath(entry);
  end
end

%% finish
%cd(myMatlabDir);

clear *Dir entries i entry conflictFiles

%clearvars
