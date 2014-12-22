% EBIRD classification

% to investigate: http://sccn.ucsd.edu/wiki/Minimalist_BCI

% common spatial pattern (DMLT; Eunho Noh)


%% load the aanalysis details

expName = 'EBIRD';

% subDir = '';
subDir = 'new_art500';
if ~isempty(subDir)
  warning('Loading data from subDir ''%s''...',subDir);
end
dataDir = fullfile(expName,'EEG','Sessions','ftpp',subDir);
% Possible locations of the data files (dataroot)
serverDir = fullfile(filesep,'Volumes','curranlab','Data');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dreamDir = fullfile(filesep,'data','projects','curranlab');
localDir = fullfile(getenv('HOME'),'data');

% pick the right dataroot
if exist('serverDir','var') && exist(serverDir,'dir')
  dataroot = serverDir;
  %runLocally = 1;
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
  %runLocally = 1;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
  %runLocally = 0;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

procDir = fullfile(dataroot,dataDir,'ft_data/data_art_nsClassic_ftAuto/tla');

subjects = {
  %'EBIRD049'; % Pilot. (due to short ses1 match, missing ses2 name)
  %'EBIRD002'; % Pilot. (due to short ses1 match, missing ses2 name)
  %'EBIRD003'; % Pilot. (due to missing ses7 name) - NB: LAST PILOT TO BE REPLACED
  %'EBIRD004'; % DNF. Dropout. Last session: 8.
  'EBIRD005';
%   %'EBIRD006'; % DNF. Dropout. Last session: 2.
%   'EBIRD007';
%   'EBIRD008';
%   'EBIRD009';
%   'EBIRD010';
%   'EBIRD011';
%   'EBIRD012';
%   %'EBIRD013'; % DNF. Dropout. Last session: 5. Lost session 6 in HD crash.
%   %'EBIRD014'; % DNF. Rejected. Last session: 1.
%   %'EBIRD015'; % DNF. Lost in HD crash.
%   %'EBIRD016'; % DNF. Lost in HD crash.
%   %'EBIRD017'; % DNF. Lost in HD crash.
%   'EBIRD018';
%   'EBIRD019';
%   'EBIRD020';
%   %'EBIRD021'; % Bad behavior
%   %'EBIRD022'; % DNF. Dropout. Last session: 8.
%   %'EBIRD023'; % DNF. Dropout. Last session: 1.
%   'EBIRD024';
%   'EBIRD025';
%   'EBIRD027';
%   'EBIRD029';
%   'EBIRD032';
%   'EBIRD034';
%   'EBIRD042';
  };

% only one cell, with all session names
% sesNames = {'session_1','session_8','session_9'};
sesNames = {'session_8'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot,[],subDir);

files.figPrintFormat = 'png';
files.saveFigs = true;

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% list the event values to analyze

% classifier asks: can we predict accurate comparison (discrimation
% between) posttest match trial stimuli? This would use data during and
% after stim2, which is when they would make the comparison.

% the classifier is trained on correctly and incorrectly named stimuli
% during training sessions. What about the fact that naming accuracy
% essentially reaches ceiling by the last training day?


% or we could ask: when during training do they start to show expert-like
% activity. This would be trained on in/accurate posttest trials and tested
% on training naming trials.



