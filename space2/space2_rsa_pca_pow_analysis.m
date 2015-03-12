% analysis of files saved using space2_rsa_pca_pow_cluster_save.m

%% load

analysisDate = '28-Feb-2015';

data_str = 'img';

% data_str = 'img_word';

% word
% thisROI = {'center109'}; % spac, ~spac x lat
% thisROI = {'LPI2','LPS','LT','RPI2','RPS','RT'}; % spac, lat
% thisROI = {'LPS','RPS'}; % spac
% thisROI = {'LT','RT'}; % spac
% thisROI = {'LPI2','RPI2'}; % spac, ~lat
% thisROI = {'LAS2','FS','RAS2'}; % spac, mem, ~spac x lat
% thisROI = {'LFP','FC','RFP'}; % spac, spac x lat
% thisROI = {'FS'};
% thisROI = {'C'};
% thisROI = {'PS'};
% thisROI = {'PI'};
% thisROI = {'LT'};
% thisROI = {'RT'};
% thisROI = {'LPS'};
% thisROI = {'RPS'};

% img
% thisROI = {'center109'}; % spac, mem x lat, ~3-way
% thisROI = {'LPI2','LPS','LT','RPI2','RPS','RT'}; % spac, mem x lat, 3-way
% thisROI = {'LPS','RPS'}; % spac, mem x lat, 3-way
% thisROI = {'LT','RT'}; % spac, ~mem x lat
% thisROI = {'LPI2','RPI2'}; % spac, ~3-way
% thisROI = {'LAS2','FS','RAS2'}; % spac
% thisROI = {'LFP','FC','RFP'}; % spac
% thisROI = {'FS'}; % none
% thisROI = {'C'}; ~spac, 3-way (PS is stronger)
% thisROI = {'PS'}; % **
% thisROI = {'PI'}; % *
% thisROI = {'LT'}; % none
% thisROI = {'RT'}; % spac
% thisROI = {'LPS'}; % lat
% thisROI = {'RPS'}; % spac

% thisROI = {'LPI2','LPS','LT','RPI2','RPS','RT'}; % spac, mem x lat, 3-way
thisROI = {'LPS','RPS'}; % spac, mem x lat, 3-way
% thisROI = {'LT','RT'}; % spac, ~mem x lat
% thisROI = {'PS'}; % **

freq = mm_freqSet('ndtools');

freqs = [ ...
  freq.theta; ...
  freq.alpha; ...
  freq.beta_lower; ...
  freq.beta_upper; ...
  freq.gamma_lower; ...
  freq.gamma_upper];

% freqs = freq.theta;
% freqs = freq.alpha;
% freqs = freq.beta_lower;
% freqs = freq.beta_upper; % frontal mem?
% freqs = freq.gamma_lower;
% freqs = freq.gamma_upper;

freq_str = sprintf('%dfreq%dto%d',size(freqs,1),round(freqs(1,1)),round(freqs(end,end)));

if iscell(thisROI)
  roi_str = sprintf(repmat('%s',1,length(thisROI)),thisROI{:});
elseif ischar(thisROI)
  roi_str = thisROI;
end

% latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
%   0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
%   0 0.3; 0.3 0.6; 0.6 0.9; ...
%   0 0.5; 0.5 1.0; ...
%   0.3 0.8; ...
%   0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
%   0 0.8; 0.1 0.9; 0.2 1.0;
%   0 1.0];

latencies = [0.0 0.2; 0.22 0.4; 0.42 0.6; 0.62 0.8; 0.82 1.0; ...
  0.1 0.3; 0.32 0.5; 0.52 0.7; 0.72 0.9; ...
  0 0.3; 0.32 0.6; 0.62 0.9; ...
  0 0.5; 0.52 1.0; ...
  0.3 0.8; ...
  0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
  0 0.8; 0.1 0.9; 0.2 1.0;
  0 1.0];

origDataType = 'pow';

% avgovertime = 'yes';
avgovertime = 'no';
avgoverfreq = 'yes';

sim_method = 'cosine';
% sim_method = 'correlation';
% sim_method = 'spearman';

% accurateClassifSelect = true;
accurateClassifSelect = false;
if accurateClassifSelect
  classif_str = 'classif';
else
  classif_str = 'noClassif';
end

% eig_criterion = 'CV85';
eig_criterion = 'kaiser';
% eig_criterion = 'analytic';

expName = 'SPACE2';

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

procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_nsClassic_ftAuto/tla');
% procDir = fullfile(dataroot,dataDir,'ft_data/cued_recall_stim_multistudy_image_multistudy_word_art_continuousICA_ftAuto/tla');

if ~isempty(data_str)
  rsaFile = fullfile(procDir,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%s_%dlat_%s_%sAvgT_%sAvgF_%s_cluster.mat',origDataType,data_str,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,avgovertime,avgoverfreq,analysisDate));
else
  rsaFile = fullfile(procDir,sprintf('RSA_PCA_%s_%s_%s_%s_%s_%dlat_%s_%sAvgT_%sAvgF_%s_cluster.mat',origDataType,sim_method,classif_str,eig_criterion,roi_str,size(latencies,1),freq_str,avgovertime,avgoverfreq,analysisDate));
end
if exist(rsaFile,'file')
  fprintf('loading: %s...',rsaFile);
  load(rsaFile);
  fprintf('Done.\n');
else
  error('does not exist: %s',rsaFile);
end

subjects_all = exper.subjects;
sesNames_all = exper.sesNames;

%% stats

% % % word
% % nTrialThresh = 6; % 27
% % nTrialThresh = 7; % 25
% nTrialThresh = 8; % 22
% % nTrialThresh = 9; % 21
% % nTrialThresh = 10; % 19
% % nTrialThresh = 12; % 18
% % nTrialThresh = 14; % 15
% % nTrialThresh = 15; % 13

% % img
nTrialThresh = 6; % 28 **
% nTrialThresh = 7; % 26 *
% nTrialThresh = 8; % 25
% nTrialThresh = 9; % 24
% nTrialThresh = 10; % 22
% nTrialThresh = 11; % 21
% nTrialThresh = 12; % 20
% nTrialThresh = 14; % 19
% nTrialThresh = 15; % 15

plotit = false;

noNans = true(length(exper.subjects),length(exper.sesStr));

passTrlThresh = true(length(exper.subjects),length(exper.sesStr));

mean_similarity = struct;
for d = 1:length(dataTypes)
  mean_similarity.(dataTypes{d}) = nan(length(subjects_all),length(sesNames_all),size(latencies,1));
  for lat = 1:size(latencies,1)
    
    for sub = 1:length(subjects_all)
      for ses = 1:length(sesNames_all)
        %   for sub = 1:length(exper.subjects)
        %     for ses = 1:length(exper.sessions)
        
        % Average Pres1--Pres2 similarity
        dataTypePairSimilarity = diag(similarity_all{sub,ses,d,lat},size(similarity_all{sub,ses,d,lat},1) / 2);
        if ~isempty(dataTypePairSimilarity)
          if length(dataTypePairSimilarity) < nTrialThresh
            passTrlThresh(sub,ses) = false;
          end
          
          mean_similarity.(dataTypes{d})(sub,ses,lat) = mean(dataTypePairSimilarity);
          
          if plotit
            figure
            imagesc(similarity_all{sub,ses,d,lat});
            colorbar;
            axis square;
            title(sprintf('%s, %.2f to %.2f',strrep(dataTypes{d},'_','-'),latencies(lat,1),latencies(lat,2)));
          end
        else
          noNans(sub,ses) = false;
        end
      end
    end
    
  end
end

% disp(mean_similarity);

fprintf('Threshold: >= %d trials. Including %d subjects.\n',nTrialThresh,sum(noNans & passTrlThresh));

%% ANOVA: factors: spaced/massed, recalled/forgotten, latency

% latencies = [0.0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1.0; ...
%   0.1 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; ...
%   0 0.3; 0.3 0.6; 0.6 0.9; ...
%   0 0.5; 0.5 1.0; ...
%   0.3 0.8; ...
%   0 0.6; 0.1 0.7; 0.2 0.8; 0.3 0.9; 0.4 1.0; ...
%   0 0.8; 0.1 0.9; 0.2 1.0;
%   0 1.0];

% % 0 to 1, in 200 ms chunks
latInd = [1:5];

% % 0.1 to 0.9, in 200 ms chunks
% latInd = [6:9];

% % 0-0.3, 0.3-0.6, 0.6-0.9
% latInd = [10:12];

% % 0-0.5, 0.5-1 *****
% latInd = [13:14];

% % 0 to 1, in 600 ms chunks
% latInd = [16:20];

% % 0 to 1 in 800 ms chunks
% latInd = [21:23];

% % 0 to 1
% latInd = 24;

% =================================================

spacings = {'spac2','spac12','spac32','mass'};
memConds = {'rc', 'fo'};
latency = cell(1,length(latInd));
latencySec = cell(1,length(latInd));
for i = 1:length(latency)
  latency{i} = sprintf('%dto%d',latencies(latInd(1)+i-1,1)*1000,latencies(latInd(1)+i-1,2)*1000);
  latencySec{i} = sprintf('%.1f-%.1f',latencies(latInd(1)+i-1,1),latencies(latInd(1)+i-1,2));
end
latStr = sprintf(repmat('_%s',1,length(latency)),latency{:});
latStr = latStr(2:end);

factorNames = {'spacings', 'memConds', 'latency'};

nVariables = nan(size(factorNames));
keepTheseFactors = false(size(factorNames));
levelNames_teg = cell(size(factorNames)); % TEG
for c = 1:length(factorNames)
  nVariables(c) = length(eval(factorNames{c}));
  levelNames_teg{c} = eval(factorNames{c}); % TEG
  if length(eval(factorNames{c})) > 1
    keepTheseFactors(c) = true;
  end
end

variableNames = cell(1,prod(nVariables));
levelNames = cell(prod(nVariables),length(factorNames));

ses=1;
nSub = sum(noNans & passTrlThresh(:,ses));
anovaData = nan(nSub,prod(nVariables));
rmaov_data_teg = nan(nSub*prod(nVariables),length(factorNames) + 2);

fprintf('Collecting %s %s ANOVA data for %d subjects:\n\t',data_str,origDataType,nSub);
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(spacings)),spacings{:}),factorNames{1});
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(memConds)),memConds{:}),factorNames{2});
fprintf('%s (%s),',latStr,factorNames{3});
fprintf('\n\tROI: %s, Freq: %s, Eig: %s, Sim: %s...',roi_str,freq_str,eig_criterion,sim_method);

lnDone = false;
vnDone = false;
subCount = 0;
rmCount = 0;
for sub = 1:length(exper.subjects)
  if all(noNans(sub,:)) && all(passTrlThresh(sub,:))
    subCount = subCount + 1;
    for ses = 1:length(exper.sesStr)
      lnCount = 0;
      vnCount = 0;
      
      for sp = 1:length(spacings)
        for mc = 1:length(memConds)
          cond_str = sprintf('%s_%s_%s',data_str,memConds{mc},spacings{sp});
          
          for lat = 1:length(latInd)
            if ~lnDone
              lnCount = lnCount + 1;
              levelNames{lnCount,1} = spacings{sp};
              levelNames{lnCount,2} = memConds{mc};
              levelNames{lnCount,3} = latency{lat};
            end
            
            vnCount = vnCount + 1;
            if ~vnDone
              variableNames{vnCount} = sprintf('Y%d',vnCount);
            end
            
            anovaData(subCount,vnCount) = mean_similarity.(cond_str)(sub,ses,latInd(lat));
            
            rmCount = rmCount + 1;
            rmaov_data_teg(rmCount,:) = [mean_similarity.(cond_str)(sub,ses,latInd(lat)) sp mc lat sub];
          end
          
        end
      end
      
      lnDone = true;
      vnDone = true;
    end
  end
end

if any(~keepTheseFactors)
  factorNames = factorNames(keepTheseFactors);
  levelNames = levelNames(:,keepTheseFactors);
  nVariables = nVariables(keepTheseFactors);
  levelNames_teg = levelNames_teg(keepTheseFactors); % TEG
  
  rmaov_data_teg = rmaov_data_teg(:,[1 (find(keepTheseFactors) + 1) size(rmaov_data_teg,2)]); % TEG
  fprintf('\n\tOnly keeping factors:%s...',sprintf(repmat(' %s',1,length(factorNames)),factorNames{:}));
end
fprintf('Done.\n');

% TEG RM ANOVA

fprintf('=======================================\n');
fprintf('This ANOVA: %s %s, ROI: %s, Freq: %s, Latency:%s, %s, %s\n\n',data_str,origDataType,roi_str,freq_str,latStr,eig_criterion,sim_method);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: %s %s, ROI: %s, Freq: %s, Latency:%s, %s, %s\n',data_str,origDataType,roi_str,freq_str,latStr,eig_criterion,sim_method);
fprintf('=======================================\n');

%% Matlab RM ANOVA

t = array2table(anovaData,'VariableNames',variableNames);

within = cell2table(levelNames,'VariableNames',factorNames);

% rm = fitrm(t,'Y1-Y8~1','WithinDesign',within);
rm = fitrm(t,sprintf('%s-%s~1',variableNames{1},variableNames{end}),'WithinDesign',within);

fprintf('=======================================\n');
fprintf('This ANOVA: %s, ROI: %s, Freq: %s, Latency:%s, %s, %s\n\n',origDataType,roi_str,freq_str,latStr,eig_criterion,sim_method);

margmean(rm,factorNames)
% grpstats(rm,factorNames)

MODELSPEC = sprintf(repmat('*%s',1,length(factorNames)),factorNames{:});
MODELSPEC = MODELSPEC(2:end);

% Perform repeated measures analysis of variance.
if any(nVariables) > 2
  [ranovatbl,A,C,D] = ranova(rm, 'WithinModel',MODELSPEC)
  %Show epsilon values
  %I assume that HF epsilon values > 1 (in R) are truncated to 1 by epsilon.m
  [ranovatbl,A,C,D] = ranova(rm, 'WithinModel',MODELSPEC);
  for cn = 1:length(C)
    tbl = epsilon(rm, C(cn))
  end
else
  [ranovatbl] = ranova(rm, 'WithinModel',MODELSPEC)
end

fprintf('Prev ANOVA: %s, ROI: %s, Freq: %s, Latency:%s, %s, %s\n',origDataType,roi_str,freq_str,latStr,eig_criterion,sim_method);
fprintf('=======================================\n');

% multcompare(rm,'spacings','By','memConds')
% multcompare(rm,'spacings','By','oldnew')
% multcompare(rm,'memConds','By','oldnew')
% % multcompare(rm,'spacings','By','oldnew','By','memConds')

pairwiseComps = nchoosek(1:length(factorNames),2);
for i = 1:size(pairwiseComps,1)
  multcompare(rm,factorNames{pairwiseComps(i,1)},'By',factorNames{pairwiseComps(i,2)})
end

%% plot RSA spacing x subsequent memory interaction

theseSub = noNans & passTrlThresh;

ses=1;

rc_mass = squeeze(mean_similarity.(sprintf('%s_rc_mass',data_str))(theseSub,ses,latInd));
fo_mass = squeeze(mean_similarity.(sprintf('%s_fo_mass',data_str))(theseSub,ses,latInd));
rc_spac2 = squeeze(mean_similarity.(sprintf('%s_rc_spac2',data_str))(theseSub,ses,latInd));
fo_spac2 = squeeze(mean_similarity.(sprintf('%s_fo_spac2',data_str))(theseSub,ses,latInd));
rc_spac12 = squeeze(mean_similarity.(sprintf('%s_rc_spac12',data_str))(theseSub,ses,latInd));
fo_spac12 = squeeze(mean_similarity.(sprintf('%s_fo_spac12',data_str))(theseSub,ses,latInd));
rc_spac32 = squeeze(mean_similarity.(sprintf('%s_rc_spac32',data_str))(theseSub,ses,latInd));
fo_spac32 = squeeze(mean_similarity.(sprintf('%s_fo_spac32',data_str))(theseSub,ses,latInd));

% mean across time
rc_mass_sub = mean(rc_mass,2);
fo_mass_sub = mean(fo_mass,2);
rc_spac2_sub = mean(rc_spac2,2);
fo_spac2_sub = mean(fo_spac2,2);
rc_spac12_sub = mean(rc_spac12,2);
fo_spac12_sub = mean(fo_spac12,2);
rc_spac32_sub = mean(rc_spac32,2);
fo_spac32_sub = mean(fo_spac32,2);

plotMeanLine = false;
plotSub = true;

if plotMeanLine
  m_mark = 'bs--';
  s2_mark = 'ro-';
  s12_mark = 'rx-';
  s32_mark = 'r^-';
else
  m_mark = 'bs';
  s2_mark = 'ro';
  s12_mark = 'rx';
  s32_mark = 'r^';
end

if plotSub
  subSpacing = 0.1;
  m_mark_sub = 'bs';
  s2_mark_sub = 'ro';
  s12_mark_sub = 'rx';
  s32_mark_sub = 'r^';
end

meanSizeR = 20;
meanSizeF = 20;

figure
hold on

if plotSub
  % forgotten
  plot((1+subSpacing)*ones(sum(theseSub(:,ses)),1), fo_mass_sub,m_mark_sub,'LineWidth',1);
  plot((1-subSpacing)*ones(sum(theseSub(:,ses)),1), fo_spac2_sub,s2_mark_sub,'LineWidth',1);
  plot((1-subSpacing)*ones(sum(theseSub(:,ses)),1), fo_spac12_sub,s12_mark_sub,'LineWidth',1);
  plot((1-subSpacing)*ones(sum(theseSub(:,ses)),1), fo_spac32_sub,s32_mark_sub,'LineWidth',1);
  
  % recalled
  plot((2+subSpacing)*ones(sum(theseSub(:,ses)),1), rc_mass_sub,m_mark_sub,'LineWidth',1);
  plot((2-subSpacing)*ones(sum(theseSub(:,ses)),1), rc_spac2_sub,s2_mark_sub,'LineWidth',1);
  plot((2-subSpacing)*ones(sum(theseSub(:,ses)),1), rc_spac12_sub,s12_mark_sub,'LineWidth',1);
  plot((2-subSpacing)*ones(sum(theseSub(:,ses)),1), rc_spac32_sub,s32_mark_sub,'LineWidth',1);
end

if plotMeanLine
  hm = plot([mean(fo_mass_sub,1) mean(rc_mass_sub,1)],m_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  hs2 = plot([mean(fo_spac2_sub,1) mean(rc_spac2_sub,1)],s2_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  hs12 = plot([mean(fo_spac12_sub,1) mean(rc_spac12_sub,1)],s12_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  hs32 = plot([mean(fo_spac32_sub,1) mean(rc_spac32_sub,1)],s32_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  
%   % recalled
%   plot(2, ,s_mark,'LineWidth',3,'MarkerSize',meanSizeR);
%   plot(2, ,m_mark,'LineWidth',3,'MarkerSize',meanSizeR);
else
  % forgotten
  plot(1, mean(fo_mass_sub,1),m_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  plot(1, mean(fo_spac2_sub,1),s2_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  plot(1, mean(fo_spac12_sub,1),s12_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  plot(1, mean(fo_spac32_sub,1),s32_mark,'LineWidth',3,'MarkerSize',meanSizeF);
  
  % recalled
  hm = plot(2, mean(rc_mass_sub,1),m_mark,'LineWidth',3,'MarkerSize',meanSizeR);
  hs2 = plot(2, mean(rc_spac2_sub,1),s2_mark,'LineWidth',3,'MarkerSize',meanSizeR);
  hs12 = plot(2, mean(rc_spac12_sub,1),s12_mark,'LineWidth',3,'MarkerSize',meanSizeR);
  hs32 = plot(2, mean(rc_spac32_sub,1),s32_mark,'LineWidth',3,'MarkerSize',meanSizeR);
end

% horiz
plot([-3 3], [0 0],'k--','LineWidth',2);

hold off
axis square
% axis([0.75 2.25 75 110]);
xlim([0.75 2.25]);
ylim([-0.61 0.61]);

set(gca,'XTick', [1 2]);
set(gca,'XTickLabel',{'Forgot','Recalled'});

ylabel('Neural Similarity');

title(sprintf('Spacing \\times Subsequent Memory: %s',data_str));
legend([hm, hs2, hs12, hs32],{'Massed','Spaced 2','Spaced 12','Spaced 32'},'Location','North');

% ticFontSize = 20;
ticFontSize = 18;
publishfig(gcf,0,ticFontSize,[],[]);

% print(gcf,'-depsc2',sprintf('~/Desktop/similarity_spacXmem_%s_%s_%s_%s_%s_%s',data_str,origDataType,roi_str,latStr,eig_criterion,sim_method));
% % print(gcf,'-dpng',sprintf('~/Desktop/similarity_spacXmem_%s_%s_%s_%s_%s_%s',data_str,origDataType,roi_str,latStr,eig_criterion,sim_method));

%% plot RSA spacing x time interaction

theseSub = noNans & passTrlThresh;

ses=1;

rc_mass = squeeze(mean_similarity.(sprintf('%s_rc_mass',data_str))(theseSub,ses,latInd));
fo_mass = squeeze(mean_similarity.(sprintf('%s_fo_mass',data_str))(theseSub,ses,latInd));
rc_spac2 = squeeze(mean_similarity.(sprintf('%s_rc_spac2',data_str))(theseSub,ses,latInd));
fo_spac2 = squeeze(mean_similarity.(sprintf('%s_fo_spac2',data_str))(theseSub,ses,latInd));
rc_spac12 = squeeze(mean_similarity.(sprintf('%s_rc_spac12',data_str))(theseSub,ses,latInd));
fo_spac12 = squeeze(mean_similarity.(sprintf('%s_fo_spac12',data_str))(theseSub,ses,latInd));
rc_spac32 = squeeze(mean_similarity.(sprintf('%s_rc_spac32',data_str))(theseSub,ses,latInd));
fo_spac32 = squeeze(mean_similarity.(sprintf('%s_fo_spac32',data_str))(theseSub,ses,latInd));

% mean across rc/fo
mass = mean(cat(3,rc_mass,fo_mass),3);
spac2 = mean(cat(3,rc_spac2,fo_spac2),3);
spac12 = mean(cat(3,rc_spac12,fo_spac12),3);
spac32 = mean(cat(3,rc_spac32,fo_spac32),3);

plotMeanLines = true;
plotSub = true;
plotSubLines = false;

if plotMeanLines
  m_mark = 'bs--';
  s2_mark = 'ro-';
  s12_mark = 'rx-';
  s32_mark = 'r^-';
else
  m_mark = 'bs';
  s2_mark = 'ro';
  s12_mark = 'rX';
  s32_mark = 'r^';
end

meanSize = 20;

if plotSub
  subSpacing = 0.1;
  if plotSubLines
    m_mark_sub = 'bs--';
    s2_mark_sub = 'ro-';
    s12_mark_sub = 'rx-';
    s32_mark_sub = 'r^-';
  else
    m_mark_sub = 'bs';
    s2_mark_sub = 'ro';
    s12_mark_sub = 'rx';
    s32_mark_sub = 'r^';
  end
end

figure
hold on

if plotSub
  if plotSubLines
    for s = 1:sum(theseSub)
      plot(mass(s,:),m_mark_sub,'LineWidth',1);
      plot(spac2(s,:),s2_mark_sub,'LineWidth',1);
      plot(spac12(s,:),s12_mark_sub,'LineWidth',1);
      plot(spac32(s,:),s32_mark_sub,'LineWidth',1);
    end
  else
    for t = 1:length(latInd)
      if plotSub
        plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), mass(:,t),m_mark_sub,'LineWidth',1);
        plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), spac2(:,t),s2_mark_sub,'LineWidth',1);
        plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), spac12(:,t),s12_mark_sub,'LineWidth',1);
        plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), spac32(:,t),s32_mark_sub,'LineWidth',1);
      end
    end
  end
end
if plotMeanLines
  hm = plot(mean(mass,1),m_mark,'LineWidth',3,'MarkerSize',meanSize);
  hs2 = plot(mean(spac2,1),s2_mark,'LineWidth',3,'MarkerSize',meanSize);
  hs12 = plot(mean(spac12,1),s12_mark,'LineWidth',3,'MarkerSize',meanSize);
  hs32 = plot(mean(spac32,1),s32_mark,'LineWidth',3,'MarkerSize',meanSize);
else
  for t = 1:length(latInd)
    hm = plot(t, mean(mass(:,t),1),m_mark,'LineWidth',3,'MarkerSize',meanSize);
    hs2 = plot(t, mean(spac2(:,t),1),s2_mark,'LineWidth',3,'MarkerSize',meanSize);
    hs12 = plot(t, mean(spac12(:,t),1),s12_mark,'LineWidth',3,'MarkerSize',meanSize);
    hs32 = plot(t, mean(spac32(:,t),1),s32_mark,'LineWidth',3,'MarkerSize',meanSize);
  end
end

% horiz
plot([-length(latInd)-1, length(latInd)+1], [0 0],'k--','LineWidth',2);

hold off
axis square
xlim([0.75 length(latInd)+0.25]);
ylim([-0.65 0.65]);

set(gca,'XTick', 1:length(latInd));
set(gca,'XTickLabel',latencySec);
xlabel('Time (Sec)');

ylabel('Neural Similarity');

title(sprintf('Spacing \\times Time: %s',data_str));
legend([hm, hs2, hs12, hs32],{'Massed','Spaced 2','Spaced 12','Spaced 32'},'Location','North');

% ticFontSize = 20;
ticFontSize = 18;
publishfig(gcf,0,ticFontSize,[],[]);

% print(gcf,'-depsc2',sprintf('~/Desktop/similarity_spacXtime_%s_%s_%s_%s_%s_%s',data_str,origDataType,roi_str,latStr,eig_criterion,sim_method));
% % print(gcf,'-dpng',sprintf('~/Desktop/similarity_spacXtime_%s_%s_%s_%s_%s_%s',data_str,origDataType,roi_str,latStr,eig_criterion,sim_method));

%% plot RSA spacing x memory x time interaction

theseSub = noNans & passTrlThresh;

ses=1;

rc_mass = squeeze(mean_similarity.(sprintf('%s_rc_mass',data_str))(theseSub,ses,latInd));
fo_mass = squeeze(mean_similarity.(sprintf('%s_fo_mass',data_str))(theseSub,ses,latInd));
rc_spac2 = squeeze(mean_similarity.(sprintf('%s_rc_spac2',data_str))(theseSub,ses,latInd));
fo_spac2 = squeeze(mean_similarity.(sprintf('%s_fo_spac2',data_str))(theseSub,ses,latInd));
rc_spac12 = squeeze(mean_similarity.(sprintf('%s_rc_spac12',data_str))(theseSub,ses,latInd));
fo_spac12 = squeeze(mean_similarity.(sprintf('%s_fo_spac12',data_str))(theseSub,ses,latInd));
rc_spac32 = squeeze(mean_similarity.(sprintf('%s_rc_spac32',data_str))(theseSub,ses,latInd));
fo_spac32 = squeeze(mean_similarity.(sprintf('%s_fo_spac32',data_str))(theseSub,ses,latInd));

plotMeanLines = true;
plotSub = false;
plotSubLines = false;

if plotMeanLines
  m_rc_mark = 'bo-';
  m_fo_mark = 'cx--';
  s2_rc_mark = 'ro-';
  s2_fo_mark = 'mx--';
  s12_rc_mark = 'ro-';
  s12_fo_mark = 'mx--';
  s32_rc_mark = 'ro-';
  s32_fo_mark = 'mx--';
else
  m_rc_mark = 'bo';
  m_fo_mark = 'cx';
  s2_rc_mark = 'ro';
  s2_fo_mark = 'mx';
  s12_rc_mark = 'ro';
  s12_fo_mark = 'mx';
  s32_rc_mark = 'ro';
  s32_fo_mark = 'mx';
end

meanSizeS = 20;
meanSizeM = 20;

if plotSub
  subSpacing = 0.2;
  if plotSubLines
    m_rc_mark_sub = 'bo-';
    m_fo_mark_sub = 'cx--';
    s2_rc_mark_sub = 'ro-';
    s2_fo_mark_sub = 'mx--';
    s12_rc_mark_sub = 'ro-';
    s12_fo_mark_sub = 'mx--';
    s32_rc_mark_sub = 'ro-';
    s32_fo_mark_sub = 'mx--';
  else
    m_rc_mark_sub = 'bo';
    m_fo_mark_sub = 'cx';
    s2_rc_mark_sub = 'ro';
    s2_fo_mark_sub = 'mx';
    s12_rc_mark_sub = 'ro';
    s12_fo_mark_sub = 'mx';
    s32_rc_mark_sub = 'ro';
    s32_fo_mark_sub = 'mx';
  end
end

figure
hold on

if plotSub
  if plotSubLines
    plotSub = false;
    for s = 1:sum(theseSub)
      % massed
      plot(rc_mass(s,:),m_rc_mark_sub,'LineWidth',1);
      plot(fo_mass(s,:),m_fo_mark_sub,'LineWidth',1);
      % spaced 2
      plot(rc_spac2(s,:),s2_rc_mark_sub,'LineWidth',1);
      plot(fo_spac2(s,:),s2_fo_mark_sub,'LineWidth',1);
      % spaced 12
      plot(rc_spac12(s,:),s12_rc_mark_sub,'LineWidth',1);
      plot(fo_spac12(s,:),s12_fo_mark_sub,'LineWidth',1);
      % spaced 32
      plot(rc_spac32(s,:),s32_rc_mark_sub,'LineWidth',1);
      plot(fo_spac32(s,:),s32_fo_mark_sub,'LineWidth',1);
    end
  else
    for t = 1:length(latInd)
      % massed
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), rc_mass(:,t),m_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), fo_mass(:,t),m_fo_mark_sub,'LineWidth',1);
      % spaced 2
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), rc_spac2(:,t),s2_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), fo_spac2(:,t),s2_fo_mark_sub,'LineWidth',1);
      % spaced 12
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), rc_spac12(:,t),s2_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), fo_spac12(:,t),s2_fo_mark_sub,'LineWidth',1);
      % spaced 32
      plot((t-subSpacing)*ones(sum(theseSub(:,ses)),1), rc_spac32(:,t),s2_rc_mark_sub,'LineWidth',1);
      plot((t+subSpacing)*ones(sum(theseSub(:,ses)),1), fo_spac32(:,t),s2_fo_mark_sub,'LineWidth',1);
    end
  end
end
if plotMeanLines
  % massed
  hmr = plot(mean(rc_mass,1),m_rc_mark,'LineWidth',3,'MarkerSize',meanSizeM);
  hmf = plot(mean(fo_mass,1),m_fo_mark,'LineWidth',3,'MarkerSize',meanSizeM);
  % spaced
  hsr2 = plot(mean(rc_spac2,1),s2_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  hsf2 = plot(mean(fo_spac2,1),s2_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  hsr12 = plot(mean(rc_spac12,1),s12_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  hsf12 = plot(mean(fo_spac12,1),s12_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  hsr32 = plot(mean(rc_spac32,1),s32_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  hsf32 = plot(mean(fo_spac32,1),s32_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
else
  for t = 1:length(latInd)
    % massed
    hmr = plot(t,mean(rc_mass(:,t),1),m_rc_mark,'LineWidth',3,'MarkerSize',meanSizeM);
    hmf = plot(t,mean(fo_mass(:,t),1),m_fo_mark,'LineWidth',3,'MarkerSize',meanSizeM);
    % spaced
    hsr2 = plot(t,mean(rc_spac2(:,t),1),s2_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    hsf2 = plot(t,mean(fo_spac2(:,t),1),s2_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    hsr12 = plot(t,mean(rc_spac12(:,t),1),s12_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    hsf12 = plot(t,mean(fo_spac12(:,t),1),s12_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    hsr32 = plot(t,mean(rc_spac32(:,t),1),s32_rc_mark,'LineWidth',3,'MarkerSize',meanSizeS);
    hsf32 = plot(t,mean(fo_spac32(:,t),1),s32_fo_mark,'LineWidth',3,'MarkerSize',meanSizeS);
  end
end

% horiz
plot([-length(latInd)-1, length(latInd)+1], [0 0],'k--','LineWidth',2);

hold off
axis square
xlim([0.75 length(latInd)+0.25]);
ylim([-0.35 0.35]);

set(gca,'XTick', 1:length(latInd));
set(gca,'XTickLabel',latencySec);
xlabel('Time (Sec)');

ylabel('Neural Similarity');

title(sprintf('Spacing \\times Memory \\times Time: %s',data_str));
% legend([hmr, hmf, hsr, hsf],{'Mass Recalled','Mass Forgot','Space Recalled','Space Forgot'},'Location','North');
legend([hmr, hmf, hsr2, hsf2, hsr12, hsf12, hsr32, hsf32],{'Mass Recalled','Mass Forgot','Space 2 Recalled','Space 2 Forgot','Space 12 Recalled','Space 12 Forgot','Space 32 Recalled','Space 32 Forgot'},'Location','North');

% ticFontSize = 20;
ticFontSize = 18;
publishfig(gcf,0,ticFontSize,[],[]);

% print(gcf,'-depsc2',sprintf('~/Desktop/similarity_spacXmemXtime_%s_%s_%s_%s_%s_%s',data_str,origDataType,roi_str,latStr,eig_criterion,sim_method));
% % print(gcf,'-dpng',sprintf('~/Desktop/similarity_spacXmemXtime_%s_%s_%s_%s_%s_%s',data_str,origDataType,roi_str,latStr,eig_criterion,sim_method));
