% % if path is not getting set properly, run these commands:
% restoredefaultpath
% rehash toolboxcache
% savepath

%% initialize

myMatlabDir = fullfile(getenv('HOME'),'Documents','MATLAB');

%% set up Psychtoolbox path

ptbDir = fullfile(filesep,'Applications','Psychtoolbox');
if ~isempty(ptbDir)
  % add top folder and all subfolders
  addpath(genpath(ptbDir));
  
  obsoleteDir = fullfile(ptbDir,'PsychHardware','iViewXToolbox','Obsolete');
  if exist(obsoleteDir,'dir')
    %fprintf('Removing %s and its subdirectories from path.\n',obsoleteDir);
    rmpath(genpath(obsoleteDir));
  end
end

%% set up eeg_toolbox path

%eeg_toolboxDir = dir(fullfile(myMatlabDir,'eeg_toolbox'));
%if ~isempty(eeg_toolboxDir)

eeg_toolboxDir = fullfile(myMatlabDir,'eeg_toolbox');
if exist(eeg_toolboxDir,'dir')
  %eeg_toolboxDir = fullfile(myMatlabDir,eeg_toolboxDir.name);
  
  % add top folder and all subfolders
  addpath(genpath(eeg_toolboxDir));
end

%% set up eeglab path

eeglabDir = dir(fullfile(myMatlabDir,'eeglab*'));
if ~isempty(eeglabDir)
  if length(eeglabDir) == 1
    eeglabDir = fullfile(myMatlabDir,eeglabDir.name);
    % add top folder and all subfolders
    addpath(genpath(eeglabDir));
    
    % remove eeglab's external directory if it was added
    eeglabExtDir = fullfile(eeglabDir,'external');
    if exist(eeglabExtDir,'dir')
      %fprintf('Removing %s and its subdirectories from path.\n',eeglabExtDir);
      rmpath(genpath(eeglabExtDir));
    end
    
    % remove eeglab's octavefunc directory if it was added
    eeglabOctDir = fullfile(eeglabDir,'functions','octavefunc');
    if exist(eeglabOctDir,'dir')
      %fprintf('Removing %s and its subdirectories from path because firls.m interferes with MATLAB''s version.\n',eeglabOctDir);
      rmpath(genpath(eeglabOctDir));
    end
    
    % % remove eeglab's fieldtrip directory if it was added
    % eeglabFtDir = dir(fullfile(eeglabDir,'external','fieldtrip*'));
    % if ~isempty(eeglabFtDir)
    %   eeglabFtDir = fullfile(myMatlabDir,eeglabFtDir.name);
    %   rmpath(genpath(eeglabFtDir));
    % end
  else
    warning('More than one eeglab* directory found, not adding to path.')
  end
else
  warning('%s not found!',eeglabDir);
end

%% set up MVPA toolbox path

mvpaDir = fullfile(myMatlabDir,'mvpa');
if exist(mvpaDir,'dir')
  % add only the top folder
  addpath(mvpaDir);
  
  % % add the subdirectories that MVPA needs
  % mvpa_add_paths;
  
  % add the subdirectories for eeg_ana toolbox
  addpath(fullfile(mvpaDir,'core','learn'));
  addpath(fullfile(mvpaDir,'core','util'));
  
  %   % functions that get in the way; structure arguments as:
  %   %
  %   % conflictFiles = {{fxn, other package},{fxn, other package}};
  %   conflictFiles = {{'isrow','MATLAB'}};
  %   for i = 1:length(conflictFiles)
  %     if exist(fullfile(mvpaDir,'afni_matlab',[conflictFiles{i}{1},'.m']),'file')
  %       fprintf('Found a conflicting MVPA function %s.m that conflicts with %s''s function. Moving MVPA''s to %s_old.m.\n',conflictFiles{i}{1},conflictFiles{i}{2},conflictFiles{i}{1});
  %       unix(sprintf('mv %s %s',fullfile(mvpaDir,'afni_matlab',[conflictFiles{i}{1},'.m']),fullfile(mvpaDir,'afni_matlab',[conflictFiles{i}{1},'_old.m'])));
  %     end
  %   end
end

%% set up fieldtrip path

%ftDir = dir(fullfile(myMatlabDir,'fieldtrip-*'));
%if ~isempty(ftDir)

ftDir = fullfile(myMatlabDir,'fieldtrip');
if exist(ftDir,'dir')
  %ftDir = fullfile(myMatlabDir,ftDir.name);
  
  % add only the top folder
  addpath(ftDir);
  % add the subdirectories that FT needs
  ft_defaults;
  
  % add the peer directory
  %fprintf('Adding %s to path.\n',fullfile(ftDir,'peer'));
  addpath(fullfile(ftDir,'peer'));
  
  % add the multivariate directory (outdated; use DMLT)
  %addpath(genpath(fullfile(ftDir,'multivariate')));
  
  % add the SPM directory
  %fprintf('Adding %s to path.\n',fullfile(ftDir,'external','spm8'));
  %addpath(fullfile(ftDir,'external','spm8'));
  
  % % remove fieldtrip's external directory
  % ftExtDir = fullfile(ftDir,'external');
  % if ~isempty(ftExtDir)
  %   fprintf('Removing %s and its subdirectories from path.\n',ftExtDir);
  %   rmpath(genpath(ftExtDir));
  % end
  
  % remove fieldtrip's eeglab directory if it was added
  ftEeglabDir = dir(fullfile(ftDir,'external','eeglab*'));
  if ~isempty(ftEeglabDir)
    ftEeglabDir = fullfile(myMatlabDir,ftEeglabDir.name);
    %fprintf('Removing %s and its subdirectories from path.\n',ftEeglabDir);
    rmpath(genpath(ftEeglabDir));
  end
  
  % functions that get in the way; structure arguments as:
  %
  % conflictFiles = {{fxn, other package},{fxn, other package}};
  conflictFiles = {{'progress','MVPA'}};
  for i = 1:length(conflictFiles)
    if exist(fullfile(ftDir,'utilities','compat',[conflictFiles{i}{1},'.m']),'file')
      fprintf('Found an old/unnecessary FieldTrip function %s.m that conflicts with %s''s function. Moving FT''s to %s_old.m.\n',conflictFiles{i}{1},conflictFiles{i}{2},conflictFiles{i}{1});
      unix(sprintf('mv %s %s',fullfile(ftDir,'utilities','compat',[conflictFiles{i}{1},'.m']),fullfile(ftDir,'utilities','compat',[conflictFiles{i}{1},'_old.m'])));
    end
  end
else
  warning('%s not found!',ftDir);
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

%epDir = dir(fullfile(myMatlabDir,'EP_Toolkit*'));
%if ~isempty(epDir)

epDir = fullfile(myMatlabDir,'EP_Toolkit');
if exist(epDir,'dir')
  %epDir = fullfile(myMatlabDir,epDir.name);
  
  % add top folder and all subfolders
  addpath(genpath(epDir));
end

%% set up EEG Analysis Toolbox path

eatDir = dir(fullfile(myMatlabDir,'eeg_ana_*'));
if ~isempty(eatDir)
  if length(eatDir) == 1
    eatDir = fullfile(myMatlabDir,eatDir.name);
    
    % add top folder and all subfolders
    addpath(genpath(eatDir));
    
    % remove some toolboxes from eeg_ana's external directory
    eatExtDir = fullfile(eatDir,'external');
    if exist(eatExtDir,'dir')
      %fprintf('Removing %s and its subdirectories from path.\n',eatExtDir);
      rmpath(genpath(fullfile(eatExtDir,'eeglab')));
      rmpath(genpath(fullfile(eatExtDir,'mvpa')));
      rmpath(genpath(fullfile(eatExtDir,'eeg_toolbox')));
    end
  else
    warning('More than one eeg_ana_* directory found, not adding to path.')
  end
end

%% add my other analysis scripts/files

addpath(genpath(fullfile(myMatlabDir,'recogmodel_mvm')));

%% add my experiment, fieldtrip, and RM ANOVA scripts

addpath(genpath(fullfile(myMatlabDir,'mat-mvm')));

%% put the ~/Documents/MATLAB folder at the top of the path

addpath(myMatlabDir);

%% remove version control directories from path (.git, .svn, CVS)

entries = regexp(path, ['[^',pathsep,']*',pathsep], 'match');
% find the version control entries
vc_entries = cell2mat(cellfun(@(x) ~isempty(strfind(x,'.git')) | ~isempty(strfind(x,'.svn')) | ~isempty(strfind(x,'CVS')), entries, 'UniformOutput', false));
% remove them
rmpath(sprintf(repmat('%s',1,sum(vc_entries)),entries{vc_entries}));

%% set better default colors using linspecer

% http://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap

% if exist('linspecer','file')
%   set(0,'DefaultFigureColormap',linspecer);
% end

% % put it back to the default:
% set(0,'DefaultFigureColormap',jet);

% % Set the color order separately for each plot. Replace '8' with the
% % number of lines you have.
% set(0,'DefaultAxesColorOrder',linspecer(8))

% % Increase line thickness to make nicer plots and better perceived color
% % difference.

% set(0,'DefaultLineLineWidth',2);
% set(0,'DefaultLineLineWidth',1.2);
set(0,'DefaultLineLineWidth',1.0);

%% clean up

%cd(myMatlabDir);

clear *Dir *entries i entry conflictFiles

%clearvars
