% EBIRD: because of how the session number got passed to seg2ft.m (or
% really, did not get passed and what was passed was misinterpreted; fixed
% on april 16 2014), the wrong events were selected for sessions > 1, and
% therefore the wrong data is in trialinfo. This script uses the data from
% the right events session to reconstruct trialinfo, even after artifacts
% have been rejected (all due to exper.badEv for ever session).

expName = 'EBIRD';

subDir = '';
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
  %'EBIRD006'; % DNF. Dropout. Last session: 2.
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
%   'EBIRD021';
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
sesNames = {'session_1','session_8','session_9'};
% sesNames = {'session_1'};
% sesNames = {'session_8'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
% replaceDataroot = false;
replaceDataroot = true;

[full_exper,full_ana,full_dirs,full_files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%% constants

evSesNames = {'pretest','posttest','posttest_delay'};
% use the collapsed phasename
phaseName = 'match';
% evVal = {'match_stim'};

old_trl_order_match_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'stimNum', 'imgCond', 'isSubord', 'trained', 'sameSpecies', 'response', 'rt', 'acc'};
new_trl_order_match_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'stimNum', 'imgCond', 'isSubord', 'trained', 'sameTrained', 'sameSpecies', 'response', 'rt', 'acc'};

% index of sameTrained data
stInd = find(ismember(new_trl_order_match_stim,{'sameTrained'}));

eegFileNameRaw = 'data_raw_match_stim.mat';
eegFileNameProc = 'data_tla_match_stim.mat';

evTypes = {'MATCH_STIM1','MATCH_STIM2'};

image_conditions = sort({'color', 'g', 'g_hi8', 'g_lo8', 'normal'});
match_responses = {'same', 'diff'};

%% processing the data: put sameTrained into EEG data.trialinfo

for sub = 1:length(full_exper.subjects)
  fprintf('Processing subject %s...\n',full_exper.subjects{sub});
  
  eventsFile = fullfile(full_dirs.dataroot,full_dirs.behDir,full_exper.subjects{sub},'events','events.mat');
  fprintf('Loading events: %s...',eventsFile);
  load(eventsFile);
  fprintf('Done.\n');
  
  for ses = 1:length(full_exper.sessions)
    fprintf('Working on session %s...\n',evSesNames{ses});
    sesEv = events.(evSesNames{ses}).(phaseName).data;
    % select only the events that we segmented
    sesEv = sesEv(ismember({sesEv.type},evTypes));
    
    % only keep the good events
    sesEv = sesEv(~full_exper.badEv.(full_exper.sesStr{ses}).match_stim{1});
    
    sesDirRaw = fullfile(full_dirs.saveDirRaw,full_exper.subjects{sub},full_exper.sesStr{ses});
    sesDirProc = fullfile(full_dirs.saveDirProc,full_exper.subjects{sub},full_exper.sesStr{ses});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the Raw subject details file
    sdFileRaw = fullfile(sesDirRaw,'subjectDetails.mat');
    fprintf('Loading raw subject details: %s...',sdFileRaw);
    load(sdFileRaw);
    fprintf('Done.\n');
    
    if length(ana.trl_order.match_stim) ~= length(new_trl_order_match_stim)
      % overwrite trl_order for match_stim with new one
      ana.trl_order.match_stim = new_trl_order_match_stim;
      
      % save the new subject details file
      fprintf('Saving raw subject details: %s...',sdFileRaw);
      save(sdFileRaw,'exper','ana','dirs','files','cfg_pp','-v7');
      fprintf('Done.\n');
    end
    clear exper ana dirs files cfg_pp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the Proc subject details file
    sdFileProc = fullfile(sesDirProc,'subjectDetails.mat');
    fprintf('Loading processed subject details: %s...',sdFileProc);
    load(sdFileProc);
    fprintf('Done.\n');
    
    if length(ana.trl_order.match_stim) ~= length(new_trl_order_match_stim)
      % overwrite trl_order for match_stim with new one
      ana.trl_order.match_stim = new_trl_order_match_stim;
      
      % save the new subject details file
      fprintf('Saving processed subject details: %s...',sdFileProc);
      save(sdFileProc,'exper','ana','dirs','files','cfg_pp','cfg_proc','-v7');
      fprintf('Done.\n');
    end
    clear exper ana dirs files cfg_pp cfg_proc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the raw EEG file (variable: data)
    eegFileRaw = fullfile(sesDirRaw,eegFileNameRaw);
    fprintf('Loading raw EEG: %s...',eegFileRaw);
    load(eegFileRaw);
    fprintf('Done.\n');
    
    %     if size(data.trialinfo,2) ~= length(new_trl_order_match_stim)
    
    % add in sameTrained status
    fprintf('Fixing trialinfo in raw EEG trialinfo...');
    
    % initialize
    new_trialinfo = zeros(size(data.trialinfo,1),length(new_trl_order_match_stim));
    
    for i = 1:length(sesEv)
      
      % eventNumber match only had 1 event type so we can hardcode this
      eventNumber = 1;
      sesType = find(ismember(full_ana.sessionNames,full_ana.sessionNames{ses}));
      phaseType = find(ismember(full_ana.phaseNames{sesType},phaseName));
      phaseCount = sesEv(i).phaseCount;
      trial = sesEv(i).trial;
      familyNum = sesEv(i).familyNum;
      speciesNum = sesEv(i).speciesNum;
      exemplarNum = sesEv(i).exemplarNum;
      if strcmp(sesEv(i).type,'MATCH_STIM1')
        stimNum = 1;
      elseif strcmp(sesEv(i).type,'MATCH_STIM2')
        stimNum = 2;
      end
      imgCond = find(ismember(image_conditions,sesEv(i).imgCond));
      isSubord = sesEv(i).isSubord;
      trained = sesEv(i).trained;
      sameTrained = sesEv(i).sameTrained;
      sameSpecies = sesEv(i).sameSpecies;
      response = ismember(match_responses,sesEv(i).resp);
      if any(response)
        response = find(response);
      elseif ~any(response) && strcmp(sesEv(i).resp,'none')
        response = 0;
      else
        keyboard
      end
      rt = sesEv(i).rt;
      acc = sesEv(i).acc;
      
      for to = 1:length(new_trl_order_match_stim)
        thisInd = find(ismember(new_trl_order_match_stim,new_trl_order_match_stim{to}));
        if ~isempty(thisInd)
          if exist(new_trl_order_match_stim{to},'var')
            new_trialinfo(i,thisInd) = eval(new_trl_order_match_stim{to});
          else
            fprintf('variable %s does not exist!\n',new_trl_order_match_stim{to});
            keyboard
          end
        end
      end
      
    end
    
    % replace the old trialinfo with the new one
    data.trialinfo = new_trialinfo;
    fprintf('Done.\n');
    
    % save the updated raw EEG file
    fprintf('Saving raw EEG: %s...',eegFileRaw);
    save(eegFileRaw,'data','-v7');
    fprintf('Done.\n');
    
    %clear data new_trialinfo
    clear data
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the processed EEG file (variable: data)
    eegFileProc = fullfile(sesDirProc,eegFileNameProc);
    fprintf('Loading processed EEG: %s...',eegFileProc);
    load(eegFileProc);
    fprintf('Done.\n');
    
    % add in sameTrained status
    fprintf('Fixing trialinfo in processed EEG trialinfo...');
    
    % replace the old trialinfo with the new one
    timelock.trialinfo = new_trialinfo;
    fprintf('Done.\n');
    
    % save the updated processed EEG file
    fprintf('Saving processed EEG: %s...',eegFileProc);
    save(eegFileProc,'timelock','-v7');
    fprintf('Done.\n');
    clear timelock new_trialinfo
    %     else
    %       fprintf('already modified! moving on.\n');
    %       clear data
    %     end
  end
end
