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
    sesEv = events.(evSesNames{ses}).(phaseName).data;
    % select only the events that we segmented
    sesEv = sesEv(ismember({sesEv.type},evTypes));
    
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
    %clear exper ana dirs files cfg_pp
    
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
    %clear exper ana dirs files cfg_pp cfg_proc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the raw EEG file (variable: data)
    eegFileRaw = fullfile(sesDirRaw,eegFileNameRaw);
    fprintf('Loading raw EEG: %s...',eegFileRaw);
    load(eegFileRaw);
    fprintf('Done.\n');
    
    if size(data.trialinfo,2) ~= length(new_trl_order_match_stim)
      % add in sameTrained status
      fprintf('Fixing trialinfo in raw EEG trialinfo...');
      
      % initialize
      new_trialinfo = zeros(size(data.trialinfo,1),length(new_trl_order_match_stim));
      
      for i = 1:size(data.trialinfo,1)
        evCount = 0;
        if ~exper.badEv.(full_exper.sesStr{ses}).match_stim{1}(i)
          evCount = evCount + 1;
          
          % eventNumber match only had 1 event type so we can hardcode this
          eventNumber = 1;
          sesType = find(ismember(ana.sessionNames,ana.sessionNames{ses}));
          phaseType = find(ismember(ana.phaseNames{sesType},phaseName));
          phaseCount = events.(evSesNames{ses}).(phaseName).data(evCount).phaseCount;
          trial = events.(evSesNames{ses}).(phaseName).data(evCount).trial;
          familyNum = events.(evSesNames{ses}).(phaseName).data(evCount).familyNum;
          speciesNum = events.(evSesNames{ses}).(phaseName).data(evCount).speciesNum;
          exemplarNum = events.(evSesNames{ses}).(phaseName).data(evCount).exemplarNum;
          if strcmp(events.(evSesNames{ses}).(phaseName).data(evCount).type,'MATCH_STIM1')
            stimNum = 1;
          elseif strcmp(events.(evSesNames{ses}).(phaseName).data(evCount).type,'MATCH_STIM2')
            stimNum = 2;
          end
          imgCond = find(ismember(image_conditions,events.(evSesNames{ses}).(phaseName).data(evCount).imgCond));
          isSubord = events.(evSesNames{ses}).(phaseName).data(evCount).isSubord;
          trained = events.(evSesNames{ses}).(phaseName).data(evCount).trained;
          sameTrained = events.(evSesNames{ses}).(phaseName).data(evCount).sameTrained;
          sameSpecies = events.(evSesNames{ses}).(phaseName).data(evCount).sameSpecies;
          response = ismember(match_responses,events.(evSesNames{ses}).(phaseName).data(evCount).resp);
          if any(response)
            response = find(response);
          elseif ~any(response) && strcmp(this_event.resp,'none')
            response = 0;
          else
            keyboard
          end
          rt = events.(evSesNames{ses}).(phaseName).data(evCount).rt;
          acc = events.(evSesNames{ses}).(phaseName).data(evCount).acc;
          
          for to = 1:length(new_trl_order_match_stim)
            thisInd = find(ismember(new_trl_order_match_stim,new_trl_order_match_stim{to}));
            if ~isempty(thisInd)
              if exist(trl_order{to},'var')
                new_trialinfo(evCount,thisInd) = eval(new_trl_order_match_stim{to});
              else
                fprintf('variable %s does not exist!\n',new_trl_order_match_stim{to});
                keyboard
              end
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
