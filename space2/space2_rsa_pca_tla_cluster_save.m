% load and re-save RSA PCA tla results.

% these results would be saved on dream using
% space2_rsa_pca_tla_classif_wrapper.m (which runs
% space2_rsa_pca_tla_classif_cluster.m).

% follow up by running space2_rsa_pca_tla_analysis.m

% can also save files processed using the hilbert method. I think the
% difference is that this tla script would not average over time, whereas
% hilbert (processed using space2_rsa_pca_tla_hilbert_classif_wrapper.m)
% would average over time.

%#rsync -avzP matt@dreamio2.colorado.edu:/data/projects/curranlab/SPACE2/EEG/Sessions/ftpp/ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_continuousICA_ftAuto/tla/RSA_PCA_tla*cluster.mat ~/data/SPACE/EEG/Sessions/ftpp/ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_continuousICA_ftAuto/tla/
%rsync -avzP matt@dreamio2.colorado.edu:/data/projects/curranlab/SPACE2/EEG/Sessions/ftpp/ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla/RSA_PCA_tla*cluster.mat ~/data/SPACE2/EEG/Sessions/ftpp/ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla/

expName = 'SPACE2';

saveDirProc = fullfile(filesep,'data','projects','curranlab',expName,'EEG/Sessions/ftpp/ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');
% saveDirProc = fullfile(filesep,'data','projects','curranlab',expName,'EEG/Sessions/ftpp/ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_continuousICA_ftAuto/tla');

subjects = {
  'SPACE2001';
  %'SPACE2002'; % really noisy EEG, finished both sessions, incl in beh
  %'SPACE2003'; % DNF session 2, exclude
  'SPACE2004';
  'SPACE2005';
  'SPACE2006';
  %'SPACE2007'; % bad performance, low trial counts, EXCLUDE
  'SPACE2008';
  %'SPACE2009'; % DNF session 2, exclude
  'SPACE2010';
  'SPACE2011';
  'SPACE2012';
  %'SPACE2013'; % didn't record EEG, stopped session 1 in middle, exclude
  'SPACE2014';
  'SPACE2015';
  %'SPACE2016'; % really noisy EEG, EXCLUDE
  'SPACE2017'; % not great performance, still including
  'SPACE2018';
  'SPACE2019';
  %'SPACE2020'; % DNF session 2, exclude
  'SPACE2021';
  'SPACE2022';
  %'SPACE2023'; % DNF session 2, exclude
  %'SPACE2024'; % DNF session 2, exclude
  %'SPACE2025'; % bad performance, low trial counts, EXCLUDE
  'SPACE2026';
  %'SPACE2027'; % really noisy EEG, DNF session 2, exclude
  'SPACE2028';
  'SPACE2029';
  'SPACE2030';
  'SPACE2031';
  'SPACE2032';
  'SPACE2033';
  'SPACE2034'; % low trial counts, EXCLUDE via threshold
  'SPACE2035'; % not great performance, still including
  'SPACE2036';
  'SPACE2037';
  'SPACE2038';
  %'SPACE2039'; % DNF session 2, exclude
  'SPACE2040';
  };

% only one cell, with all session names
sesNames = {'session_1_session_2'};

analysisDate = '22-Feb-2015';

origDataType = 'tla';
% origDataType = 'hilbert';

% avgovertime = 'yes';
avgovertime = 'no';

% accurateClassifSelect = true;
accurateClassifSelect = false;
if accurateClassifSelect
  classif_str = 'classif';
else
  classif_str = 'noClassif';
end

lpfilt = false;
if lpfilt
  lpfreq = 40;
%   lofiltord = 3;
%   lpfilttype = 'but';
  lpfilt_str = sprintf('lpfilt%d',lpfreq);
else
  lpfilt_str = 'lpfiltNo';
end

dataTypes = {'img_rc_mass','img_fo_mass', 'img_rc_spac2','img_fo_spac2', 'img_rc_spac12','img_fo_spac12', 'img_rc_spac32','img_fo_spac32'};

% dataTypes = {'word_rc_mass','word_fo_mass', 'word_rc_spac2','word_fo_spac2', 'word_rc_spac12','word_fo_spac12', 'word_rc_spac32','word_fo_spac32'};

% dataTypes = {'img_rc_mass','img_fo_mass', 'img_rc_spac2','img_fo_spac2', 'img_rc_spac12','img_fo_spac12', 'img_rc_spac32','img_fo_spac32', 'word_rc_mass','word_fo_mass', 'word_rc_spac2','word_fo_spac2', 'word_rc_spac12','word_fo_spac12', 'word_rc_spac32','word_fo_spac32'};

if all(ismember({'img_rc_mass','img_fo_mass', 'img_rc_spac2','img_fo_spac2', 'img_rc_spac12','img_fo_spac12', 'img_rc_spac32','img_fo_spac32', 'word_rc_mass','word_fo_mass', 'word_rc_spac2','word_fo_spac2', 'word_rc_spac12','word_fo_spac12', 'word_rc_spac32','word_fo_spac32'},dataTypes))
  data_str = 'img_word';
elseif all(ismember({'img_rc_mass','img_fo_mass', 'img_rc_spac2','img_fo_spac2', 'img_rc_spac12','img_fo_spac12', 'img_rc_spac32','img_fo_spac32'},dataTypes))
  data_str = 'img';
elseif all(ismember({'word_rc_mass','word_fo_mass', 'word_rc_spac2','word_fo_spac2', 'word_rc_spac12','word_fo_spac12', 'word_rc_spac32','word_fo_spac32'},dataTypes))
  data_str = 'word';
end

% allROIs = {{'LPI2','LPS','LT','RPI2','RPS','RT'},{'center109'},{'LPS','RPS'},{'LT','RT'},{'LPI2','RPI2'},{'LAS2','FS','RAS2'},{'LFP','FC','RFP'}};
allROIs = {{'FS'},{'C'},{'PS'},{'PI'},{'LT'},{'RT'},{'LPS'},{'RPS'}};

allLats = {[0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
  0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
  0 0.3; 0.3 0.6; 0.6 0.9; ...
  0 0.5; 0.5 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0]};

if strcmp(origDataType,'tla')
  allFreqs = {[]};
elseif strcmp(origDataType,'hilbert')
  % allFreqs = {[4 8; 8 12; 12 30; 30 50]};
  % allFreqs = {[4 8] [8 12] [12 30] [30 50]};
  % allFreqs = {[4 8; 8 12; 12 30; 30 50] [4 8] [8 12] [12 30] [30 50]};
  % allFreqs = {[3 7; 8 12; 13 20; 23 30; 32 47; 51 80] [3 7] [8 12] [13 20] [23 30] [32 47] [51 80]};
  freq = mm_freqSet('ndtools');
  allFreqs = {[ ...
    freq.theta; ...
    freq.alpha; ...
    freq.beta_lower; ...
    freq.beta_upper; ...
    freq.gamma_lower; ...
    freq.gamma_upper] ...
    freq.theta ...
    freq.alpha ...
    freq.beta_lower ...
    freq.beta_upper ...
    freq.gamma_lower ...
    freq.gamma_upper ...
    };
end

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
      if strcmp(origDataType,'tla')
        freq_str = 'noFreq';
      elseif strcmp(origDataType,'hilbert')
        freq_str = sprintf('%dfreq%dto%d',size(freqs,1),round(freqs(1,1)),round(freqs(end,end)));
      end
      
      for s = 1:length(allSimMethod)
        sim_method = allSimMethod{s};
        
        for e = 1:length(allEigCrit)
          eig_criterion = allEigCrit{e};
          
          % initialize
          similarity_all = cell(length(subjects),length(sesNames),length(dataTypes),size(latencies,1));
          similarity_ntrials = nan(length(subjects),length(sesNames),length(dataTypes),size(latencies,1));
          
          for sub = 1:length(subjects)
            for ses = 1:length(sesNames)
              if strcmp(origDataType,'tla')
                %savedFile = fullfile(saveDirProc,subjects{sub},sesNames{ses},sprintf('RSA_PCA_%s_%s_%s_%s_%s_%s_%dlat_%sAvgT_%s.mat',origDataType,data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),avgovertime,analysisDate));
                savedFile = fullfile(saveDirProc,subjects{sub},sesNames{ses},sprintf('RSA_PCA_%s_%s_%s_%s_%s_%s_%dlat_%sAvgT_%s_%s.mat',origDataType,data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),avgovertime,lpfilt_str,analysisDate));
              elseif strcmp(origDataType,'hilbert')
                savedFile = fullfile(saveDirProc,subjects{sub},sesNames{ses},sprintf('RSA_PCA_%s_%s_%s_%s_%s_%s_%dlat_%s_%sAvgT_%s.mat',origDataType,data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,avgovertime,analysisDate));
              end
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
          
          if strcmp(origDataType,'tla')
            saveFile = fullfile(saveDirProc,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%s_%dlat_%sAvgT_%s_%s_cluster.mat',origDataType,data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),cfg_sel.avgovertime,lpfilt_str,analysisDate));
          elseif strcmp(origDataType,'hilbert')
            saveFile = fullfile(saveDirProc,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%s_%dlat_%s_%sAvgT_%s_cluster.mat',origDataType,data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,cfg_sel.avgovertime,analysisDate));
          end
          
          fprintf('Saving %s...\n',saveFile);
          if strcmp(origDataType,'tla')
            save(saveFile,'exper','dataTypes','thisROI','cfg_sel','eig_criterion','latencies','similarity_all','similarity_ntrials');
          elseif strcmp(origDataType,'hilbert')
            save(saveFile,'exper','dataTypes','thisROI','cfg_sel','eig_criterion','latencies','freqs','similarity_all','similarity_ntrials');
          end
          fprintf('Done.\n');
          
        end
      end
    end
  end
end
