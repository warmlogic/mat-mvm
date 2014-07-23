% % if path is not getting set properly, run these commands:
% restoredefaultpath
% rehash toolboxcache
% savepath

%% initialize

myMatlabDir = fullfile(getenv('HOME'),'Documents','MATLAB');

%% set up Psychtoolbox path

% https://psychtoolbox.org

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

%% set up MNE path

% http://martinos.org/mne/

mneDir = dir(fullfile(filesep,'Applications','MNE-*'));
if ~isempty(mneDir)
  setenv('MNE_ROOT',fullfile(filesep,'Applications',mneDir.name));
  
  mneDir = fullfile(filesep,'Applications',mneDir.name,'share','matlab');
  
  if exist(mneDir,'dir')
    % add top folder
    addpath(mneDir);
  end
end

% mnehome = getenv('MNE_ROOT');
% mnematlab = sprintf('%s/share/matlab',mnehome);
% if exist(mnematlab,'dir')
%   addpath(mnematlab);
% end
% clear mnehome mnematlab

%% set up eeg_toolbox path

% http://memory.psych.upenn.edu/Software

%eeg_toolboxDir = dir(fullfile(myMatlabDir,'eeg_toolbox'));
%if ~isempty(eeg_toolboxDir)

eeg_toolboxDir = fullfile(myMatlabDir,'eeg_toolbox');
if exist(eeg_toolboxDir,'dir')
  %eeg_toolboxDir = fullfile(myMatlabDir,eeg_toolboxDir.name);
  
  % add top folder and all subfolders
  addpath(genpath(eeg_toolboxDir));
end

%% set up eeglab path

% http://sccn.ucsd.edu/eeglab/

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

% https://code.google.com/p/princeton-mvpa-toolbox/

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

% http://fieldtrip.fcdonders.nl

% find the right FT directory
ftDir = fullfile(myMatlabDir,'fieldtrip');
if ~exist(ftDir,'dir')
  ftDir = dir(fullfile(myMatlabDir,'fieldtrip-*'));
  if ~isempty(ftDir)
    ftDir = fullfile(myMatlabDir,ftDir.name);
  else
    warning('No FieldTrip directory found!');
    ftDir = fullfile(myMatlabDir,'fieldtrip');
  end
end

if exist(ftDir,'dir')
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
  % if ~isempty(ftExtDir) && ~isempty(strfind(path,ftExtDir))
  %   fprintf('Removing %s and its subdirectories from path.\n',ftExtDir);
  %   rmpath(genpath(ftExtDir));
  % end
  
  % add java file for reading EGI MFF files
  %
  % see: http://fieldtrip.fcdonders.nl/getting_started/egi
  egi_mff_jar = dir(fullfile(ftDir,'external','egi_mff','java','MFF-*.jar'));
  if ~isempty(egi_mff_jar) && length(egi_mff_jar) == 1
    javaaddpath({fullfile(ftDir,'external','egi_mff','java',egi_mff_jar(1).name)});
    clear egi_mff_jar
  end
  
  % remove fieldtrip's external eeglab directory if it was added
  ftEeglabDir = dir(fullfile(ftDir,'external','eeglab*'));
  if ~isempty(ftEeglabDir)
    ftEeglabDir = fullfile(ftDir,'external',ftEeglabDir.name);
    %fprintf('Removing %s and its subdirectories from path.\n',ftEeglabDir);
    if ~isempty(strfind(path,ftEeglabDir))
      rmpath(genpath(ftEeglabDir));
    end
  end
  
  % remove fieldtrip's external MNE directory if it was added
  ftMneDir = fullfile(ftDir,'external','mne');
  if exist(ftMneDir,'dir') && ~isempty(strfind(path,ftMneDir))
    %fprintf('Removing %s and its subdirectories from path.\n',ftMneDir);
    rmpath(genpath(ftMneDir));
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

% https://github.com/distrep/DMLT

dmltDir = fullfile(myMatlabDir,'DMLT');
if exist(dmltDir,'dir')
  % add top folder and all subfolders
  addpath(genpath(dmltDir));
  
  % remove the murphy toolbox because pca conflicts with matlab's function
  dlmtExtDir = fullfile(dmltDir,'external');
  if exist(dlmtExtDir,'dir')
    %fprintf('Removing %s and its subdirectories from path.\n',fullfile(dlmtExtDir,'murphy'));
    rmpath(genpath(fullfile(dlmtExtDir,'murphy')));
  end
end

% if you get an error because libgfortran is not installed, follow the
% instructions from here <http://www.cs.ubc.ca/~hoffmanm/matlab.html> to
% download the libgfortran 4.2 tarball from here
% <http://r.research.att.com/tools/> listed under OS X 10.6 and copy
% libgfortran.* files to /usr/local/lib.

%% set up EP_Toolkit path

% http://sourceforge.net/projects/erppcatoolkit/

%epDir = dir(fullfile(myMatlabDir,'EP_Toolkit*'));
%if ~isempty(epDir)

epDir = fullfile(myMatlabDir,'EP_Toolkit');
if exist(epDir,'dir')
  %epDir = fullfile(myMatlabDir,epDir.name);
  
  % add top folder and all subfolders
  addpath(genpath(epDir));
end

%% set up EEG Analysis Toolbox path

% https://code.google.com/p/eeg-analysis-toolbox/

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

%% set up Phillip Gilley's NDTools path

% http://spot.colorado.edu/~gilley/ndtools.html

ndDir = fullfile(myMatlabDir,'ndtools');
if exist(ndDir,'dir')
  % add top folder and all subfolders
  addpath(genpath(ndDir));
end

%% add my other analysis scripts/files

rm_mvmDir = fullfile(myMatlabDir,'recogmodel_mvm');

if exist(rm_mvmDir,'dir')
  % add top folder and all subfolders
  addpath(genpath(rm_mvmDir));
end

%% add my experiment, fieldtrip, and RM ANOVA scripts

% https://github.com/warmlogic/mat-mvm/

mat_mvmDir = fullfile(myMatlabDir,'mat-mvm');

if exist(mat_mvmDir,'dir')
  % add top folder and all subfolders
  addpath(genpath(mat_mvmDir));
end

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
