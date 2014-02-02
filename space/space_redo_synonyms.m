% regrade synonyms

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
  %   'SPACE021';
  %   'SPACE022';
  %   'SPACE027';
  %   'SPACE029';
  };

procDir = '/Users/matt/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';

sesNames = {'session_1'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);



%% move the events file to the old version

for sub = 1:length(subjects)
  eventsDir = fullfile(dirs.dataroot,behDir,subjects{sub},'events');
  
  if exist(fullfile(eventsDir,'events.mat'),'file')
    unix(sprintf('mv %s %s',fullfile(eventsDir,'events.mat'),fullfile(eventsDir,'events_old.mat')));
  end
  
end

%% regrade subjects

space_prepData_events(subjects);

%% load EEG and transfer info

sessionName = 'oneDay';

% eegFilesToTransferInfo = {};

for sub = 1:length(subjects)
  eventsDir = fullfile(dirs.dataroot,dirs.behDir,subjects{sub},'events');
  
  load(fullfile(eventsDir,'events.mat'));
  
  if ~ismember(0.5,unique([events.oneDay.cued_recall.data.recall_spellCorr]))
    
    for ses = 1:length(sesNames)
      % load the exposure phase file - data_tla_expo_stim.mat
      type = 'EXPO_IMAGE';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      cols.cr_recog_acc = 13;
      cols.cr_recall_resp = 14;
      cols.cr_recall_spellCorr = 15;
      
      eegFile_expo = fullfile(procDir,subjects{sub},sesNames{ses},'data_tla_expo_stim.mat');
      load(eegFile_expo);
      
      phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s...',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          timelock.trialinfo(i,cols.cr_recog_acc) = thisEv.cr_recog_acc;
          timelock.trialinfo(i,cols.cr_recall_resp) = thisEv.cr_recall_resp;
          timelock.trialinfo(i,cols.cr_recall_spellCorr) = thisEv.cr_recall_spellCorr;
        end
      end
      
      fprintf('resaving EEG...');
      save(eegFile_expo,'timelock');
      fprintf('Done.\n')
      
      
      
      
      % load the multistudy phase file - data_tla_multistudy_image.mat
      type = 'STUDY_IMAGE';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      cols.cr_recog_acc = 14;
      cols.cr_recall_resp = 15;
      cols.cr_recall_spellCorr = 16;
      
      eegFile_msImg = fullfile(procDir,subjects{sub},'data_tla_multistudy_image.mat');
      load(eegFile_msImg);
      
      phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s...',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          timelock.trialinfo(i,cols.cr_recog_acc) = thisEv.cr_recog_acc;
          timelock.trialinfo(i,cols.cr_recall_resp) = thisEv.cr_recall_resp;
          timelock.trialinfo(i,cols.cr_recall_spellCorr) = thisEv.cr_recall_spellCorr;
        end
      end
      fprintf('resaving EEG...');
      save(eegFile_msImg,'timelock');
      fprintf('Done.\n')
      
      
      
      
      
      
      % load the multistudy phase file - data_tla_multistudy_word.mat
      eegFile_msWord = fullfile(procDir,subjects{sub},'data_tla_multistudy_word.mat');
      load(eegFile_msWord);
      
      type = 'STUDY_WORD';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      cols.cr_recog_acc = 14;
      cols.cr_recall_resp = 15;
      cols.cr_recall_spellCorr = 16;
      
      eegFile_msImg = fullfile(procDir,subjects{sub},'data_tla_multistudy_image.mat');
      load(eegFile_msImg);
      
      phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s...',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          timelock.trialinfo(i,cols.cr_recog_acc) = thisEv.cr_recog_acc;
          timelock.trialinfo(i,cols.cr_recall_resp) = thisEv.cr_recall_resp;
          timelock.trialinfo(i,cols.cr_recall_spellCorr) = thisEv.cr_recall_spellCorr;
        end
      end
      
      fprintf('resaving EEG...');
      save(eegFile_msWord,'timelock');
      fprintf('Done.\n')
      
      
      
      
      
      
      % load the cued recall phase file - data_tla_cued_recall_stim.mat
      eegFile_cr = fullfile(procDir,subjects{sub},'data_tla_cued_recall_stim.mat');
      load(eegFile_cr);
      
      type = 'RECOGTEST_STIM';
      
      cols = struct;
      cols.phaseName = 3;
      cols.phaseCount = 4;
      cols.trial = 5;
      cols.stimNum = 6;
      %cols.cr_recog_acc = 13;
      cols.cr_recall_resp = 18;
      cols.cr_recall_spellCorr = 19;
      
      eegFile_msImg = fullfile(procDir,subjects{sub},'data_tla_multistudy_image.mat');
      load(eegFile_msImg);
      
      phaseName = ana.phaseNames{ses}{timelock.trialinfo(1,cols.phaseName)};
      fprintf('\t%s...',phaseName);
      
      for i = 1:size(timelock.trialinfo,1)
        thisEv = events.(sessionName).(phaseName).data(...
          strcmp({events.(sessionName).(phaseName).data.phaseName},phaseName) & ...
          [events.(sessionName).(phaseName).data.phaseCount] == timelock.trialinfo(i,cols.phaseCount) & ...
          [events.(sessionName).(phaseName).data.trial] == timelock.trialinfo(i,cols.trial) & ...
          [events.(sessionName).(phaseName).data.stimNum] == timelock.trialinfo(i,cols.stimNum) & ...
          ismember({events.(sessionName).(phaseName).data.type},type) ...
          );
        
        if length(thisEv) == 1
          %timelock.trialinfo(i,cols.cr_recog_acc) = thisEv.cr_recog_acc;
          timelock.trialinfo(i,cols.cr_recall_resp) = thisEv.cr_recall_resp;
          timelock.trialinfo(i,cols.cr_recall_spellCorr) = thisEv.cr_recall_spellCorr;
        end
      end
      
      fprintf('resaving EEG...');
      save(eegFile_cr,'timelock');
      fprintf('Done.\n')
    end
  end
end



