% ebird behavioral analysis

%% load the analysis details

procDir = '/Users/matt/data/EBIRD/EEG/Sessions/ftpp/ft_data/cued_recall_stim_expo_stim_multistudy_image_multistudy_word_art_ftManual_ftICA/tla';

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

% % only one cell, with all session names
% sesNames = {'session_1','session_2','session_3','session_4','session_5','session_6','session_7','session_8','session_9'};
% 
% % replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
% replaceDataroot = false;
% 
% [exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot);

% %% decide who to kick out based on trial counts
% 
% % Subjects with bad behavior
% exper.badBehSub = {{}};
% 
% % exclude subjects with low event counts
% [exper,ana] = mm_threshSubs_multiSes(exper,ana,5,[],'vert');

%% collapsed or not

% collapsePhases = true;
collapsePhases = false;

%% load the behavioral data

if collapsePhases
  behfile = 'EBIRD_behav_results_collapsed.mat';
else
  behfile = 'EBIRD_behav_results.mat';
end

load(fullfile(getenv('HOME'),'data','EBIRD','Behavioral','Sessions',behfile));

exper.badSub = false(length(subjects),1);

%% ttest

alpha = 0.05;
tails = 'both';

phase = 'match_1';

trainCond = 'trained';
% trainCond = 'untrained';

imgCond = 'all';
% imgCond = 'color';
% imgCond = 'g';
% imgCond = 'g_hi8';
% imgCond = 'g_lo8';
% imgCond = 'normal';

% level = 'basic';
level = 'subord';

measure = 'dp';

ses = 'pretest';
% ses = 'posttest';

if strcmp(imgCond,'all')
  data1 = results.(ses).(phase).(trainCond).(level).(measure)(~exper.badSub);
else
  data1 = results.(ses).(phase).(trainCond).(imgCond).(level).(measure)(~exper.badSub);
end
data1_str = sprintf('%s %s %s %s img:%s %s',ses,phase,trainCond,level,imgCond,measure);

ses = 'posttest';
% ses = 'posttest_delay';

if strcmp(imgCond,'all')
  data2 = results.(ses).(phase).(trainCond).(level).(measure)(~exper.badSub);
else
  data2 = results.(ses).(phase).(trainCond).(imgCond).(level).(measure)(~exper.badSub);
end
data2_str = sprintf('%s %s %s %s img:%s %s',ses,phase,trainCond,level,imgCond,measure);

d = mm_effect_size('within',data1,data2);
[h, p, ci, stats] = ttest(data1,data2,alpha,tails);

fprintf('%s (M=%.2f; SEM=%.2f) vs\t%s (M=%.2f; SEM=%.2f):\n\tt(%d)=%.2f, d=%.2f, SD=%.2f, SEM=%.2f, p=%.5f\n', ...
  data1_str, ...
  mean(data1), ...
  std(data1) / sqrt(length(data1)), ...
  data2_str, ...
  mean(data2), ...
  std(data2) / sqrt(length(data2)), ...
  stats.df, ...
  stats.tstat, ...
  d, ...
  std(data1 - data2),...
  std(data1 - data2) / sqrt(length(data1)),...
  p);

%% RMANOVA - prepost

% cnds = {'img_RgH_rc_spac_word_RgH_rc_spac', 'img_RgH_fo_spac_word_RgH_fo_spac', 'img_RgH_rc_mass_word_RgH_rc_mass', 'img_RgH_fo_mass_word_RgH_fo_mass'};
% 
% nThresh = 1;
% 
% goodSub = ones(length(exper.subjects),1);
% 
% for ses = 1:length(exper.sesStr)
%   for cnd = 1:length(cnds)
%     goodSub = goodSub .* D.(cnds{cnd}).nTrial(:,ses) >= nThresh;
%   end
% end

sesNames = {'pretest','posttest','posttest_delay'};
% trainConds = {'trained','untrained'};
% trainConds = {'TT','UU','TU','UT'};
trainConds = {'TT','UU'};
imgConds = {'normal','color','g','g_hi8','g_lo8'};
famLevel = {'basic','subord'};

phaseName = 'match_1';

measure = 'dp';

anovaData = [];

for sub = 1:length(subjects)
  
  if ~exper.badSub(sub)
    
    theseData = [];
    
    for ses = 1:length(sesNames)
      for trn = 1:length(trainConds)
        for img = 1:length(imgConds)
          for fam = 1:length(famLevel)
            
            theseData = cat(2,theseData,results.(sesNames{ses}).(phaseName).(trainConds{trn}).(imgConds{img}).(famLevel{fam}).(measure)(sub));
            
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'session','training','imgConds', 'basicSubord'};
O = teg_repeated_measures_ANOVA(anovaData, [length(sesNames) length(trainConds) length(imgConds) length(famLevel)], varnames);

%% write to file, with header

anovaFile = fullfile(getenv('HOME'),'data','EBIRD','Behavioral','Sessions','ANOVA',sprintf('EBIRD_ANOVA_prepost_%s.txt',measure));
fid = fopen(anovaFile,'w+');

thisHeader = [];
for i = 1:length(sesNames)
  thisHeader = cat(2,thisHeader,sprintf('%s%s',sesNames{i},repmat(sprintf('\t'),1,length(trainConds) * length(imgConds) * length(famLevel))));
end
fprintf(fid,'\t%s\n',thisHeader);

thisHeader = [];
for i = 1:length(trainConds)
  thisHeader = cat(2,thisHeader,sprintf('%s%s',trainConds{i},repmat(sprintf('\t'),1,length(imgConds) * length(famLevel))));
end
thisHeader = repmat(thisHeader,1,length(sesNames));
fprintf(fid,'\t%s\n',thisHeader);

thisHeader = [];
for i = 1:length(imgConds)
  thisHeader = cat(2,thisHeader,sprintf('%s%s',imgConds{i},repmat(sprintf('\t'),1,length(famLevel))));
end
thisHeader = repmat(thisHeader,1,length(sesNames) * length(trainConds));
fprintf(fid,'\t%s\n',thisHeader);

thisHeader = [];
for i = 1:length(famLevel)
  thisHeader = cat(2,thisHeader,sprintf('%s%s',famLevel{i},repmat(sprintf('\t'),1,length(famLevel) - 1)));
end
thisHeader = repmat(thisHeader,1,length(sesNames) * length(trainConds) * length(imgConds));
fprintf(fid,'\t%s\n',thisHeader);

for i = 1:size(anovaData,1)
  tData = sprintf('%s',subjects{i});
  %tData = sprintf('%.4f',anovaData(i,1));
  tData = sprintf('%s%s\n',tData,sprintf(repmat('\t%.4f',1,length(anovaData(i,:))),anovaData(i,:)));
  
  fprintf(fid,'%s',tData);
end

fclose(fid);

%% RMANOVA - train - need to collapse phases

% cnds = {'img_RgH_rc_spac_word_RgH_rc_spac', 'img_RgH_fo_spac_word_RgH_fo_spac', 'img_RgH_rc_mass_word_RgH_rc_mass', 'img_RgH_fo_mass_word_RgH_fo_mass'};
% 
% nThresh = 1;
% 
% goodSub = ones(length(exper.subjects),1);
% 
% for ses = 1:length(exper.sesStr)
%   for cnd = 1:length(cnds)
%     goodSub = goodSub .* D.(cnds{cnd}).nTrial(:,ses) >= nThresh;
%   end
% end

sesNames = {'train1','train2','train3','train4','train5','train6'};
% trainConds = {'trained','untrained'};
% imgConds = {'color','g','g_hi8','g_lo8','normal'};
famLevel = {'basic','subord'};

% phaseNames = {'name_1','name_2','name_3','name_4'};
phaseName = 'name';

measure = 'hr';
% measure = 'rt_hit';

anovaData = [];

for sub = 1:length(subjects)
  
  if ~exper.badSub(sub)
    
    theseData = [];
    
    for ses = 1:length(sesNames)
      %for trn = 1:length(trainConds)
      %for img = 1:length(imgConds)
      for fam = 1:length(famLevel)
        
        %theseData = cat(2,theseData,results.(sesNames{ses}).(phaseName).(trainConds{trn}).(imgConds{img}).(famLevel{fam}).(measure)(sub));
        theseData = cat(2,theseData,results.(sesNames{ses}).(phaseName).(famLevel{fam}).(measure)(sub));
        
      end
      %end
      %end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'session', 'basicSubord'};
O = teg_repeated_measures_ANOVA(anovaData, [length(sesNames) length(famLevel)], varnames);

%% write to file, with header

anovaFile = fullfile(getenv('HOME'),'data','EBIRD','Behavioral','Sessions','ANOVA',sprintf('EBIRD_ANOVA_train_%s.txt',measure));
fid = fopen(anovaFile,'w+');

thisHeader = [];
for i = 1:length(sesNames)
  thisHeader = cat(2,thisHeader,sprintf('%s%s',sesNames{i},repmat(sprintf('\t'),1,length(famLevel))));
end
fprintf(fid,'\t%s\n',thisHeader);

% thisHeader = [];
% for i = 1:length(trainConds)
%   thisHeader = cat(2,thisHeader,sprintf('%s%s',trainConds{i},repmat(sprintf('\t'),1,length(imgConds) * length(famLevel))));
% end
% thisHeader = repmat(thisHeader,1,length(sesNames));
% fprintf(fid,'\t%s\n',thisHeader);
% 
% thisHeader = [];
% for i = 1:length(imgConds)
%   thisHeader = cat(2,thisHeader,sprintf('%s%s',imgConds{i},repmat(sprintf('\t'),1,length(famLevel))));
% end
% thisHeader = repmat(thisHeader,1,length(sesNames) * length(trainConds));
% fprintf(fid,'\t%s\n',thisHeader);

thisHeader = [];
for i = 1:length(famLevel)
  thisHeader = cat(2,thisHeader,sprintf('%s%s',famLevel{i},repmat(sprintf('\t'),1,length(famLevel) - 1)));
end
thisHeader = repmat(thisHeader,1,length(sesNames));
fprintf(fid,'\t%s\n',thisHeader);

for i = 1:size(anovaData,1)
  tData = sprintf('%s',subjects{i});
  %tData = sprintf('%.4f',anovaData(i,1));
  tData = sprintf('%s%s\n',tData,sprintf(repmat('\t%.4f',1,length(anovaData(i,:))),anovaData(i,:)));
  
  fprintf(fid,'%s',tData);
end

fclose(fid);



%% RMANOVA - pre-test vs post-test x basic/subord x training, for each image manipulation group

% cnds = {'img_RgH_rc_spac_word_RgH_rc_spac', 'img_RgH_fo_spac_word_RgH_fo_spac', 'img_RgH_rc_mass_word_RgH_rc_mass', 'img_RgH_fo_mass_word_RgH_fo_mass'};
% 
% nThresh = 1;
% 
% goodSub = ones(length(exper.subjects),1);
% 
% for ses = 1:length(exper.sesStr)
%   for cnd = 1:length(cnds)
%     goodSub = goodSub .* D.(cnds{cnd}).nTrial(:,ses) >= nThresh;
%   end
% end

% trainConds = {'TT','UU'};
% imgConds = {'normal','color','g'};

% trainConds = {'TT','UU','TU','UT'};
trainConds = {'TT','UU'};
imgConds = {'g','g_hi8','g_lo8'};

% sesNames = {'pretest','posttest','posttest_delay'};

% sesNames = {'posttest', 'pretest'};
% sesNames = {'posttest_delay', 'pretest'};

sesDiff = {{'posttest','pretest'},{'posttest_delay','pretest'}};

famLevel = {'basic','subord'};

phaseName = 'match_1';

measure = 'dp';

anovaData = [];

for sub = 1:length(subjects)
  
  if ~exper.badSub(sub)
    
    theseData = [];
    
    %for ses = 1:length(sesNames)
    for ses = 1:length(sesDiff)
      for trn = 1:length(trainConds)
        for img = 1:length(imgConds)
          for fam = 1:length(famLevel)
            
            %dp_diff = results.(sesNames{1}).(phaseName).(trainConds{trn}).(imgConds{img}).(famLevel{fam}).(measure)(sub) - results.(sesNames{2}).(phaseName).(trainConds{trn}).(imgConds{img}).(famLevel{fam}).(measure)(sub);
            dp_diff = results.(sesDiff{ses}{1}).(phaseName).(trainConds{trn}).(imgConds{img}).(famLevel{fam}).(measure)(sub) - results.(sesDiff{ses}{2}).(phaseName).(trainConds{trn}).(imgConds{img}).(famLevel{fam}).(measure)(sub);
            
            theseData = cat(2,theseData,dp_diff);
            
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'sesDiff', 'training', 'imgConds', 'basicSubord'};
O = teg_repeated_measures_ANOVA(anovaData, [length(sesDiff) length(trainConds) length(imgConds) length(famLevel)], varnames);


