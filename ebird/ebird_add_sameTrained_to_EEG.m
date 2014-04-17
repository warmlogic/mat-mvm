% EBIRD: add status of whether each stimulus presented during the match
% phase (MATCH_STIM1 and MATCH_STIM2, as well as MATCH_RESP) had the same
% training status. This allows us to figure out if it was a
% trained-trained, untrained-untrained, trained-untrained, or
% untrained-trained trial from only a single stimulus.

% Question: what happens to MATCH_RESP for the trained field? there are
% two entries, one for MATCH_STIM1 and one for MATCH_STIM2.
%
% Answer: only segmnted MATCH_STIM1 and MATCH_STIM2 (not MATCH_RESP).
% Therefore, this is not an issue.

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
  'EBIRD007';
  'EBIRD008';
  'EBIRD009';
  'EBIRD010';
  'EBIRD011';
  'EBIRD012';
  %'EBIRD013'; % DNF. Dropout. Last session: 5. Lost session 6 in HD crash.
  %'EBIRD014'; % DNF. Rejected. Last session: 1.
  %'EBIRD015'; % DNF. Lost in HD crash.
  %'EBIRD016'; % DNF. Lost in HD crash.
  %'EBIRD017'; % DNF. Lost in HD crash.
  'EBIRD018';
  'EBIRD019';
  'EBIRD020';
  'EBIRD021';
  %'EBIRD022'; % DNF. Dropout. Last session: 8.
  %'EBIRD023'; % DNF. Dropout. Last session: 1.
  'EBIRD024';
  'EBIRD025';
  'EBIRD027';
  'EBIRD029';
  'EBIRD032';
  'EBIRD034';
  'EBIRD042';
  };

% only one cell, with all session names
sesNames = {'session_1','session_8','session_9'};
% sesNames = {'session_1'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = false;
% replaceDataroot = true;

[full_exper,full_ana,full_dirs,full_files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%% constants

evSesNames = {'pretest','posttest','posttest_delay'};
% use the collapsed phasename
phaseName = 'match';

old_trl_order_match_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'stimNum', 'imgCond', 'isSubord', 'trained', 'sameSpecies', 'response', 'rt', 'acc'};
new_trl_order_match_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'stimNum', 'imgCond', 'isSubord', 'trained', 'sameTrained', 'sameSpecies', 'response', 'rt', 'acc'};

% index of sameTrained data
stInd = find(ismember(new_trl_order_match_stim,{'sameTrained'}));

eegFileNameRaw = 'data_raw_match_stim.mat';
eegFileNameProc = 'data_tla_match_stim.mat';

%% processing the data: put sameTrained into EEG data.trialinfo

for sub = 1:length(full_exper.subjects)
  fprintf('Processing subject %s...\n',full_exper.subjects{sub});
  
  eventsFile = fullfile(full_dirs.dataroot,full_dirs.behDir,full_exper.subjects{sub},'events','events.mat');
  fprintf('Loading events: %s...',eventsFile);
  load(eventsFile);
  fprintf('Done.\n');
  
  for ses = 1:length(full_exper.sessions)
    sesEv = events.(evSesNames{ses}).(phaseName).data;
    
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
    
    if size(data.trialinfo,2) ~= length(new_trl_order_match_stim)
      % add in sameTrained status
      fprintf('Adding sameTrained status to raw EEG trialinfo...');
      
      % put in another column
      new_trialinfo = cat(2,data.trialinfo(:,1:stInd-1),zeros(size(data.trialinfo,1),1),data.trialinfo(:,stInd:end));
      
      for i = 1:size(data.trialinfo,1)
        phaseCount = data.trialinfo(i,ismember(old_trl_order_match_stim,'phaseCount'));
        trial = data.trialinfo(i,ismember(old_trl_order_match_stim,'trial'));
        exemplarNum = data.trialinfo(i,ismember(old_trl_order_match_stim,'exemplarNum'));
        isSubord = data.trialinfo(i,ismember(old_trl_order_match_stim,'isSubord'));
        stimNum = data.trialinfo(i,ismember(old_trl_order_match_stim,'stimNum'));
        if stimNum == 1
          type = 'MATCH_STIM1';
        elseif stimNum == 2
          type = 'MATCH_STIM2';
        else
          fprintf('stimNum does not match either 1 or 2\n');
          keyboard
        end
        
        % find the corresponding event
        thisEv = sesEv([sesEv.phaseCount] == phaseCount & [sesEv.trial] == trial & [sesEv.exemplarNum] == exemplarNum & [sesEv.isSubord] == isSubord & ismember({sesEv.type},type));
        if length(thisEv) == 1
          % put sameTrained in the right place in new_trialinfo
          new_trialinfo(i,stInd) = thisEv.sameTrained;
        elseif length(thisEv) > 1
          fprintf('found too many events\n');
          keyboard
        elseif isempty(thisEv)
          fprintf('did not find any events\n');
          keyboard
        else
          fprintf('what happened?\n');
          keyboard
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
      fprintf('Adding sameTrained status to processed EEG trialinfo...');
      
      %     % put in another column
      %     new_trialinfo = cat(2,timelock.trialinfo(:,1:stInd-1),zeros(size(timelock.trialinfo,1),1),timelock.trialinfo(:,stInd:end));
      %
      %     for i = 1:size(timelock.trialinfo,1)
      %       phaseCount = timelock.trialinfo(i,ismember(old_trl_order_match_stim,'phaseCount'));
      %       trial = timelock.trialinfo(i,ismember(old_trl_order_match_stim,'trial'));
      %       exemplarNum = timelock.trialinfo(i,ismember(old_trl_order_match_stim,'exemplarNum'));
      %       isSubord = timelock.trialinfo(i,ismember(old_trl_order_match_stim,'isSubord'));
      %       stimNum = timelock.trialinfo(i,ismember(old_trl_order_match_stim,'stimNum'));
      %       if stimNum == 1
      %         type = 'MATCH_STIM1';
      %       elseif stimNum == 2
      %         type = 'MATCH_STIM2';
      %       else
      %         fprintf('stimNum does not match either 1 or 2\n');
      %         keyboard
      %       end
      %
      %       % find the corresponding event
      %       thisEv = sesEv([sesEv.phaseCount] == phaseCount & [sesEv.trial] == trial & [sesEv.exemplarNum] == exemplarNum & [sesEv.isSubord] == isSubord & ismember({sesEv.type},type));
      %       if length(thisEv) == 1
      %         % put sameTrained in the right place in new_trialinfo
      %         new_trialinfo(i,stInd) = thisEv.sameTrained;
      %       elseif length(thisEv) > 1
      %         fprintf('found too many events\n');
      %         keyboard
      %       elseif isempty(thisEv)
      %         fprintf('did not find any events\n');
      %         keyboard
      %       else
      %         fprintf('what happened?\n');
      %         keyboard
      %       end
      %     end
      
      % replace the old trialinfo with the new one
      timelock.trialinfo = new_trialinfo;
      fprintf('Done.\n');
      
      % save the updated processed EEG file
      fprintf('Saving processed EEG: %s...',eegFileProc);
      save(eegFileProc,'timelock','-v7');
      fprintf('Done.\n');
      clear timelock new_trialinfo
    else
      fprintf('already modified! moving on.\n');
      clear data
    end
  end
end
