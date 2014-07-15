% load and re-save RSA PCA pow results

%rsync -avzP matt@dreamio2.colorado.edu:/data/projects/curranlab/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow/RSA_PCA*cluster.mat ~/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow/

expName = 'SPACE';
saveDirProc = fullfile(filesep,'data','projects','curranlab',expName,'EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/pow');

subjects = {
  %'SPACE001'; % low trial counts
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
  %'SPACE017'; % old assessment: really noisy EEG, half of ICA components rejected
  'SPACE018';
  %'SPACE019';
  'SPACE020';
  'SPACE021';
  'SPACE022';
  'SPACE027';
  'SPACE029';
  'SPACE037';
  %'SPACE039'; % noisy EEG; original EEG analyses stopped here
  'SPACE023';
  'SPACE024';
  'SPACE025';
  'SPACE026';
  'SPACE028';
  %'SPACE030'; % low trial counts
  'SPACE032';
  'SPACE034';
  'SPACE047';
  'SPACE049';
  'SPACE036';
  };

% only one cell, with all session names
sesNames = {'session_1'};

analysisDate = '13-Jul-2014';

origDataType = 'pow';

avgovertime = 'yes';
avgoverfreq = 'yes';

% accurateClassifSelect = true;
accurateClassifSelect = false;
if accurateClassifSelect
  classif_str = 'classif';
else
  classif_str = 'noClassif';
end

dataTypes = {'img_RgH_rc_spac', 'img_RgH_rc_mass','img_RgH_fo_spac', 'img_RgH_fo_mass'};

allROIs = {{'LPI2','LPS','LT','RPI2','RPS','RT'},{'center109'},{'LPS','RPS'},{'LT','RT'},{'LPI2','RPI2'},{'LAS','FC','RAS'}};

allLats = {[0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
  0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
  0 0.3; 0.3 0.6; 0.6 0.9; ...
  0 0.5; 0.5 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0]};

% allFreqs = {[4 8; 8 12; 12 30; 30 50]};
% allFreqs = {[4 8] [8 12] [12 30] [30 50]};
% allFreqs = {[4 8; 8 12; 12 30; 30 50] [4 8] [8 12] [12 30] [30 50]};
allFreqs = {[3 8; 8 12; 12 30; 30 50] [3 8] [8 12] [12 30] [30 50]};

allSimMethod = {'cosine'};

% sim_method = 'cosine';
% sim_method = 'correlation';
% sim_method = 'spearman';
% sim_method = 'euclidean';

allEigCrit = {'kaiser'};
% allEigCrit = {'CV85','kaiser'};

% % keep components with eigenvalue >= 1
% eig_criterion = 'kaiser';

% % compute the percent explained variance expected from each component if
% % all events are uncorrelated with each other; keep it if above this level.
% % So, each component would explain 100/n, where n is the number of
% % events/components.
% eig_criterion = 'analytic';

% keep components that cumulatively explain at least 85% of the variance
% eig_criterion = 'CV85';

for r = 1:length(allROIs)
  thisROI = allROIs{r};
  if iscell(thisROI)
    roi_str = sprintf(repmat('%s',1,length(thisROI)),thisROI{:});
  elseif ischar(thisROI)
    roi_str = thisROI;
  end
  
  for l = 1:length(allLats)
    latencies = allLats{l};
    
    for f = 1:length(allFreqs)
      freqs = allFreqs{f};
      freq_str = sprintf('%dfreq%dto%d',size(freqs,1),freqs(1,1),freqs(end,end));
      
      for s = 1:length(allSimMethod)
        sim_method = allSimMethod{s};
        
        for e = 1:length(allEigCrit)
          eig_criterion = allEigCrit{e};
          
          % initialize
          similarity_all = cell(length(subjects),length(sesNames),length(dataTypes),size(latencies,1));
          similarity_ntrials = nan(length(subjects),length(sesNames),length(dataTypes),size(latencies,1));
          
          for sub = 1:length(subjects)
            for ses = 1:length(sesNames)
              savedFile = fullfile(saveDirProc,subjects{sub},sesNames{ses},sprintf('RSA_PCA_%s_%s_%s_%s_%s_%dlat_%s_%sAvgT_%sAvgF_%s.mat',origDataType,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,avgovertime,avgoverfreq,analysisDate));
              if exist(savedFile,'file')
                fprintf('Loading %s...\n',savedFile);
                subData = load(savedFile);
                fprintf('Done.\n');
              else
                error('Does not exist: %s',savedFile);
              end
              
              exper = subData.exper;
              cfg_sel = subData.cfg_sel;
              
              for d = 1:length(dataTypes)
                for lat = 1:size(latencies,1)
                  similarity_all{sub,ses,d,lat} = subData.similarity_all{1,1,d,lat};
                  similarity_ntrials(sub,ses,d,lat) = subData.similarity_ntrials(1,1,d,lat);
                end
              end
              
            end
          end
          
          exper.subjects = subjects;
          exper.sesNames = sesNames;
          
          saveFile = fullfile(saveDirProc,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%dlat_%s_%sAvgT_%sAvgF_%s_cluster.mat',origDataType,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,cfg_sel.avgovertime,cfg_sel.avgoverfreq,analysisDate));
          fprintf('Saving %s...\n',saveFile);
          save(saveFile,'exper','dataTypes','thisROI','cfg_sel','eig_criterion','latencies','freqs','similarity_all','similarity_ntrials');
          fprintf('Done.\n');
          
        end
      end
    end
  end
end
