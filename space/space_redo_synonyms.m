% regrade recall responses, taking synonyms into account (given half
% credit, so they can be lumped in with correct or incorrect)

subjects = {
  'SPACE001';
  'SPACE002';
  'SPACE003';
  'SPACE004';
  'SPACE005';
  'SPACE006';
  'SPACE007';
  %'SPACE008'; % didn't perform task correctly, didn't perform well
  'SPACE009';
  'SPACE010';
  'SPACE011';
  'SPACE012';
  'SPACE013';
  'SPACE014';
  'SPACE015';
  'SPACE016';
  'SPACE017'; % old assessment: really noisy EEG, half of ICA components rejected
  'SPACE018';
  'SPACE019';
  'SPACE020';
  %   %'SPACE021'; % synonyms already scored from here
  %   %'SPACE022';
  %   %'SPACE027';
  %   %'SPACE029';
  };

procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';

sesNames = {'session_1'};
sessionName = 'oneDay';

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

%% backup the events file - rename old version

for sub = 1:length(subjects)
  fprintf('%s...',subjects{sub});
  eventsDir = fullfile(dirs.dataroot,dirs.behDir,subjects{sub},'events');
  
  if exist(fullfile(eventsDir,'events.mat'),'file')
    fprintf('Backing up events file...');
    [s] = unix(sprintf('mv %s %s',fullfile(eventsDir,'events.mat'),fullfile(eventsDir,'events_syn_backup.mat')));
    if s ~= 0
      warning('file was not moved correctly!');
    else
      fprintf('Done.\n');
    end
  else
    keyboard
  end
  
end

%% regrade subjects

space_prepData_events(subjects);

%% load EEG and transfer info

for sub = 1:length(subjects)
  fprintf('%s...',subjects{sub});
  
  fprintf('Loading events...');
  eventsDir = fullfile(dirs.dataroot,dirs.behDir,subjects{sub},'events');
  load(fullfile(eventsDir,'events.mat'));
  fprintf('Done.\n');
  
  for ses = 1:length(sesNames)
    %% load the exposure phase file - data_tla_expo_stim.mat
    
    phaseName = 'expo';
    crField = 'cr_recall_spellCorr';
    
    if any([events.(sessionName).(phaseName).data.(crField)] == 0.5)
      
      type = 'EXPO_IMAGE';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      %cols.cr_recog_acc = 13;
      %cols.cr_recall_resp = 14;
      cols.(crField) = 15;
      
      eegFile_expo = fullfile(procDir,subjects{sub},sesNames{ses},'data_tla_expo_stim.mat');
      load(eegFile_expo);
      
      %phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s...\n',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          if thisEv.(crField) == 0.5
            fprintf('\tsynonym\n');
            %keyboard
          end
          
          %timelock.trialinfo(i,cols.cr_recog_acc) = thisEv.cr_recog_acc;
          %timelock.trialinfo(i,cols.cr_recall_resp) = thisEv.cr_recall_resp;
          timelock.trialinfo(i,cols.(crField)) = thisEv.(crField);
        else
          fprintf('\tlength(thisEv) ~= 1\n');
          keyboard
        end
      end
      
      fprintf('\tresaving EEG...');
      save(eegFile_expo,'timelock');
      fprintf('Done.\n')
    else
      fprintf('\tNo synonyms found!\n');
    end
    
    %% load the multistudy phase file - data_tla_multistudy_image.mat
    
    phaseName = 'multistudy';
    crField = 'cr_recall_spellCorr';
    
    if any([events.(sessionName).(phaseName).data.(crField)] == 0.5)
      type = 'STUDY_IMAGE';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      %cols.cr_recog_acc = 14;
      %cols.cr_recall_resp = 15;
      cols.(crField) = 16;
      
      eegFile_msImg = fullfile(procDir,subjects{sub},sesNames{ses},'data_tla_multistudy_image.mat');
      load(eegFile_msImg);
      
      %phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s (image)...\n',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          if thisEv.(crField) == 0.5
            fprintf('\tsynonym\n');
            %keyboard
          end
          
          %timelock.trialinfo(i,cols.cr_recog_acc) = thisEv.cr_recog_acc;
          %timelock.trialinfo(i,cols.cr_recall_resp) = thisEv.cr_recall_resp;
          timelock.trialinfo(i,cols.(crField)) = thisEv.(crField);
        else
          fprintf('\tlength(thisEv) ~= 1\n');
          keyboard
        end
      end
      
      fprintf('\tresaving EEG...');
      save(eegFile_msImg,'timelock');
      fprintf('Done.\n')
    else
      fprintf('\tNo synonyms found!\n');
    end
    
    %% load the multistudy phase file - data_tla_multistudy_word.mat
    
    phaseName = 'multistudy';
    crField = 'cr_recall_spellCorr';
    
    if any([events.(sessionName).(phaseName).data.(crField)] == 0.5)
      type = 'STUDY_WORD';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      %cols.cr_recog_acc = 14;
      %cols.cr_recall_resp = 15;
      cols.(crField) = 16;
      
      eegFile_msWord = fullfile(procDir,subjects{sub},sesNames{ses},'data_tla_multistudy_word.mat');
      load(eegFile_msWord);
      
      %phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s (word)...\n',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          if thisEv.(crField) == 0.5
            fprintf('\tsynonym\n');
            %keyboard
          end
          
          %timelock.trialinfo(i,cols.cr_recog_acc) = thisEv.cr_recog_acc;
          %timelock.trialinfo(i,cols.cr_recall_resp) = thisEv.cr_recall_resp;
          timelock.trialinfo(i,cols.(crField)) = thisEv.(crField);
        else
          fprintf('\tlength(thisEv) ~= 1\n');
          keyboard
        end
      end
      
      fprintf('\tresaving EEG...');
      save(eegFile_msWord,'timelock');
      fprintf('Done.\n')
    else
      fprintf('\tNo synonyms found!\n');
    end
    
    %% load the cued recall phase file - data_tla_cued_recall_stim.mat
    
    phaseName = 'cued_recall';
    crField = 'recall_spellCorr';
    
    if any([events.(sessionName).(phaseName).data.(crField)] == 0.5)
      type = 'RECOGTEST_STIM';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      %cols.recog_acc = 13;
      %cols.recall_resp = 18;
      cols.(crField) = 19;
      
      eegFile_cr = fullfile(procDir,subjects{sub},sesNames{ses},'data_tla_cued_recall_stim.mat');
      load(eegFile_cr);
      
      %phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s...\n',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          if thisEv.(crField) == 0.5
            fprintf('\tsynonym\n');
            %keyboard
          end
          
          %if ~isempty(thisEv.recall_resp) && ~ismember(thisEv.recall_resp,{'NO_RESPONSE'})
          %  recall_resp = 1;
          %else
          %  recall_resp = 0;
          %end
          
          %timelock.trialinfo(i,cols.recog_acc) = thisEv.recog_acc;
          %timelock.trialinfo(i,cols.recall_resp) = recall_resp;
          timelock.trialinfo(i,cols.(crField)) = thisEv.(crField);
        else
          fprintf('\tlength(thisEv) ~= 1\n');
          keyboard
        end
      end
      
      fprintf('\tresaving EEG...');
      save(eegFile_cr,'timelock');
      fprintf('Done.\n')
    else
      fprintf('\tNo synonyms found!\n');
    end
  end
end

