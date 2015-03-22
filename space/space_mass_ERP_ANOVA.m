
%% initialize

allPeakInfo = struct;

cfg = [];

%% gather data

nPoints_sub_latency = 10;

% spacings = {'mass', 'spac', 'onePres'};
% % oldnew = {'p1'};
% oldnew = {'p2'};
% memConds = {'all'};

% spacings = {'mass', 'spac'};
% oldnew = {'p1', 'p2'};
% memConds = {'all'};

% didn't test new words, so can't assess memory, but can use p1
spacings = {'mass', 'spac'};
oldnew = {'p1', 'p2'};
memConds = {'rc','fo'};

erpComponents = {'LPC','N400','N1'};
% erpComponents = {'N1'};

% roi = { ...
%   {{'LPS2','RPS2'}} ...
%   {{'C'}} ...
%   {{'E50','E51','E57','E58','E59','E64','E65','E90','E91','E95','E96','E97','E100','E101'}} ...
%   };
% roi = { ...
%   {{'RPS2'}} ...
%   {{'C'}} ...
%   {{'E50','E51','E57','E58','E59','E64','E65'}} ... % E58
%   };

roi = { ...
  {{'E62','E72','E76','E77','E78','E84','E85'}} ... % E77
  {{'C'}} ...
  {{'E50','E51','E57','E58','E59','E64','E65'}} ... % E58
  };
% roi = { ...
%   {{'E62','E72','E76','E77','E78','E84','E85'}} ... % E77
%   };

cfg.nPoints_sub = nPoints_sub_latency;

% make sure roi has the same length as erpComponents
if length(erpComponents) ~= length(roi)
  error('roi must have the same length as erpComponents');
end

lpcPeak = 0.596; % centered on E77
% lpcPeak = 0.576; % LPS2+RPS2
% lpcPeak = 0.500; % RPS2
% lpcPeak = 0.600; % general

% n400Peak = 0.360; % FS2
n400Peak = 0.372; % C

n1Peak = 0.172;

for sp = 1:length(spacings)
  
  for on = 1:length(oldnew)
    
    for mc = 1:length(memConds)
      
      cond_str = [];
      
      if strcmp(spacings{sp},'onePres')
        % % single presentation or first presentation
        if strcmp(memConds{mc},'all');
          cfg.conditions = {'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1'};
          %cfg.conditions = {'word_onePres'};
          
          cond_str = sprintf('%s_%s',spacings{sp},memConds{mc});
        end
      elseif strcmp(spacings{sp},'mass') || strcmp(spacings{sp},'spac')
        
        if strcmp(memConds{mc},'all');
          cfg.conditions = {sprintf('word_RgH_rc_%s_%s',spacings{sp},oldnew{on}),sprintf('word_RgH_fo_%s_%s',spacings{sp},oldnew{on})};
        else
          cfg.conditions = {sprintf('word_RgH_%s_%s_%s',memConds{mc},spacings{sp},oldnew{on})};
        end
        
        cond_str = sprintf('%s_%s_%s',spacings{sp},oldnew{on},memConds{mc});
      end
      
      fprintf('================================================\n');
      fprintf('Condition: %s\n',cond_str);
      disp(cfg.conditions);
      fprintf('================================================\n');
      
      for er = 1:length(erpComponents)
        for r = 1:length(roi{er})
          
          cfg.roi = roi{er}{r};
          roi_str = sprintf('%s%s',erpComponents{er},sprintf(repmat('_%s',1,length(cfg.roi)),cfg.roi{:}));
          
          if strcmp(erpComponents{er},'LPC')
            % LPC
            cfg.order = 'descend'; % descend = positive peaks first
            %cfg.latency = [lpcPeak-0.05 lpcPeak+0.05]; % LPC - around GA peak (space+mass) +/- 50
            cfg.latency = [lpcPeak-0.1 lpcPeak+0.1]; % LPC - around GA peak (space+mass) +/- 100
            %cfg.latency = [lpcPeak-0.15 lpcPeak+0.15]; % LPC - around GA peak (space+mass) +/- 150
          elseif strcmp(erpComponents{er},'N400')
            % N400
            cfg.order = 'ascend'; % ascend = negative peaks first
            cfg.latency = [n400Peak-0.05 n400Peak+0.05]; % N400 - around GA peak (space+mass) +/- 50
            %cfg.latency = [n400Peak-0.1 n400Peak+0.1]; % N400 - around GA peak (space+mass) +/- 100
          elseif strcmp(erpComponents{er},'N1')
            cfg.order = 'ascend'; % ascend = negative peaks first
            cfg.latency = [n1Peak-0.05 n1Peak+0.05]; % around GA peak (space+mass) +/- 50
            %cfg.latency = [n1Peak-0.1 n1Peak+0.1]; % around GA peak (space+mass) +/- 100
          end
          
          % % average across electrodes
          cfg.datadim = 'time';
          % and time points
          cfg.avgovertime = true;
          
          cfg.is_ga = false;
          cfg.outputSubjects = true;
          % cfg.is_ga = true;
          cfg.sesNum = 1;
          
          cfg.plotit = false;
          cfg.voltlim = [-3 3]; % LPC
          % cfg.voltlim = [-2 2]; % N400
          % cfg.voltlim = [-1 5];
          
          % peakInfo = mm_findPeak(cfg,ana,exper,ga_tla);
          peakInfo = mm_findPeak(cfg,ana,exper,data_tla);
          
          allPeakInfo.(cond_str).(roi_str) = peakInfo;
        end
      end
    end
  end
end

%% ANOVA: factors: spaced/massed, recalled/forgotten, old/new, ROI

% spacings = {'mass', 'spac', 'onePres'};
% oldnew = {'p2'};
% memConds = {'all'};

% spacings = {'mass', 'spac'};
% oldnew = {'p1', 'p2'};
% memConds = {'all'};

% didn't test new words, so can't assess memory, but can use p1
spacings = {'mass', 'spac'};
oldnew = {'p1', 'p2'};
% oldnew = {'p2'};
memConds = {'rc','fo'};

measure = 'latency';
% measure = 'voltage';

if strcmp(measure,'latency')
  %nPoints = cfg.nPoints_sub;
  nPoints = 10;
else
  nPoints = 1;
end

% erpComp = 'N1';
% roi = {'E50_E51_E57_E58_E59_E64_E65'}; % centered on E58/T5

% erpComp = 'N400';
% roi = {'C'}; % centered on Cz
% % roi = {'FS2'}; % centered on E6

erpComp = 'LPC';
% roi = {'RPS2'}; % centered on E85
roi = {'E62_E72_E76_E77_E78_E84_E85'}; % centered on E77

factorNames = {'spacings', 'oldnew', 'memConds', 'roi'};
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

anovaData = nan(sum(~exper.badSub),prod(nVariables));
rmaov_data_teg = [];

fprintf('Getting data for %s%s %s ANOVA...',erpComp,sprintf(repmat(' %s',1,length(roi)),roi{:}),measure);

lnDone = false;
vnDone = false;
for sub = 1:sum(~exper.badSub)
  lnCount = 0;
  vnCount = 0;
  %theseData = [];
  
  for sp = 1:length(spacings)
    for on = 1:length(oldnew)
      for mc = 1:length(memConds)
        cond_str = [];
        if strcmp(spacings{sp},'onePres')
          % single presentation or first presentation
          if strcmp(memConds{mc},'all');
            cond_str = sprintf('%s_%s',spacings{sp},memConds{mc});
          end
        elseif strcmp(spacings{sp},'mass') || strcmp(spacings{sp},'spac')
          cond_str = sprintf('%s_%s_%s',spacings{sp},oldnew{on},memConds{mc});
        end
        
        for r = 1:length(roi)
          roi_str = sprintf('%s_%s',erpComp,roi{r});
          
          if ~lnDone
            lnCount = lnCount + 1;
            levelNames{lnCount,1} = spacings{sp};
            levelNames{lnCount,2} = oldnew{on};
            levelNames{lnCount,3} = memConds{mc};
            levelNames{lnCount,4} = roi{r};
          end
          
          vnCount = vnCount + 1;
          if ~vnDone
            variableNames{vnCount} = sprintf('Y%d',vnCount);
          end
          
          anovaData(sub,vnCount) = mean(allPeakInfo.(cond_str).(roi_str).subjects.(measure)(sub,nPoints));
          
          rmaov_data_teg = cat(1,rmaov_data_teg,[allPeakInfo.(cond_str).(roi_str).subjects.(measure)(sub,1) sp on mc r sub]);
        end
      end
    end
  end
  lnDone = true;
  vnDone = true;
end

if any(~keepTheseFactors)
  factorNames = factorNames(keepTheseFactors);
  levelNames = levelNames(:,keepTheseFactors);
  nVariables = nVariables(keepTheseFactors);
  levelNames_teg = levelNames_teg(keepTheseFactors); % TEG
  
  rmaov_data_teg = rmaov_data_teg(:,[1 (find(keepTheseFactors) + 1) size(rmaov_data_teg,2)]);
end

fprintf('Done.\n');

% TEG ANOVA

fprintf('================================================================\n');
fprintf('This ANOVA: %s:%s, %s\n',erpComp,sprintf(repmat(' %s',1,length(roi)),roi{:}),measure);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: %s:%s, %s\n',erpComp,sprintf(repmat(' %s',1,length(roi)),roi{:}),measure);
fprintf('================================================================\n');

%% Matlab ANOVA

t = array2table(anovaData,'VariableNames',variableNames);

within = cell2table(levelNames,'VariableNames',factorNames);

% rm = fitrm(t,'Y1-Y8~1','WithinDesign',within);
rm = fitrm(t,sprintf('%s-%s~1',variableNames{1},variableNames{end}),'WithinDesign',within);

fprintf('================================================================\n');
fprintf('This ANOVA: %s:%s, %s\n',erpComp,sprintf(repmat(' %s',1,length(roi)),roi{:}),measure);

margmean(rm,factorNames)
% grpstats(rm,factorNames)

% Perform repeated measures analysis of variance.
if any(nVariables) > 2
  [ranovatbl,A,C,D] = ranova(rm, 'WithinModel','spacings*oldnew*memConds')
  %Show epsilon values
  %I assume that HF epsilon values > 1 (in R) are truncated to 1 by epsilon.m
  [ranovatbl,A,C,D] = ranova(rm, 'WithinModel','spacings*oldnew*memConds');
  for cn = 1:length(C)
    tbl = epsilon(rm, C(cn))
  end
else
  [ranovatbl] = ranova(rm, 'WithinModel','spacings*oldnew*memConds')
end

fprintf('Prev ANOVA: %s:%s, %s\n',erpComp,sprintf(repmat(' %s',1,length(roi)),roi{:}),measure);
fprintf('================================================================\n');

% multcompare(rm,'spacings','By','memConds')
% multcompare(rm,'spacings','By','oldnew')
% multcompare(rm,'memConds','By','oldnew')
% % multcompare(rm,'spacings','By','oldnew','By','memConds')

pairwiseComps = nchoosek(1:length(factorNames),2);
for i = 1:size(pairwiseComps,1)
  multcompare(rm,factorNames{pairwiseComps(i,1)},'By',factorNames{pairwiseComps(i,2)})
end

%% gather data for pairwise t-tests

erpComp = 'LPC_E62_E72_E76_E77_E78_E84_E85';
% erpComp = 'N400_C';
% erpComp = 'N1_E50_E51_E57_E58_E59_E64_E65';

spac_str = 'spac_p2_rc';
mass_str = 'mass_p2_rc';

fprintf('%s\n',erpComp);

space_peak = allPeakInfo.(spac_str).(erpComp).subjects.latency(:,1);
space_volt = allPeakInfo.(spac_str).(erpComp).subjects.voltage(:,1);

space_mem_volt = allPeakInfo.spac_p2_rc.(erpComp).subjects.voltage(:,1) - allPeakInfo.spac_p2_fo.(erpComp).subjects.voltage(:,1);

mass_peak = allPeakInfo.(mass_str).(erpComp).subjects.latency(:,1);
mass_volt = allPeakInfo.(mass_str).(erpComp).subjects.voltage(:,1);

mass_mem_volt = allPeakInfo.mass_p2_rc.(erpComp).subjects.voltage(:,1) - allPeakInfo.mass_p2_fo.(erpComp).subjects.voltage(:,1);

spacEff_rc_volt = allPeakInfo.spac_p2_rc.(erpComp).subjects.voltage(:,1) - allPeakInfo.mass_p2_rc.(erpComp).subjects.voltage(:,1);

spac_repEff_rc_volt = allPeakInfo.spac_p2_all.(erpComp).subjects.voltage(:,1) - allPeakInfo.onePres_all.(erpComp).subjects.voltage(:,1);
mass_repEff_rc_volt = allPeakInfo.mass_p2_all.(erpComp).subjects.voltage(:,1) - allPeakInfo.onePres_all.(erpComp).subjects.voltage(:,1);

% one_peak = allPeakInfo.once_all.(erpComp).subjects.latency(:,1);
% one_volt = allPeakInfo.once_all.(erpComp).subjects.voltage(:,1);

%% ttest - latency

fprintf('space: %.4f sec\n',mean(space_peak));
fprintf('mass: %.4f sec\n',mean(mass_peak));
% fprintf('one: %.4f sec\n',mean(one_peak));

[h,p,ci,stats] = ttest(space_peak,mass_peak,'alpha',0.05,'tail','both');
fprintf('Space vs mass: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

% [h,p,ci,stats] = ttest(space_peak,one_peak,'alpha',0.05,'tail','both');
% fprintf('Space vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

% [h,p,ci,stats] = ttest(mass_peak,one_peak,'alpha',0.05,'tail','both');
% fprintf('Mass vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

%% ttest - voltage

fprintf('space: %.4f uV\n',mean(space_volt));
fprintf('mass: %.4f uV\n',mean(mass_volt));
% fprintf('one: %.4f uV\n',mean(one_volt));

[h,p,ci,stats] = ttest(space_volt,mass_volt,'alpha',0.05,'tail','both');
fprintf('Space vs mass: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

% [h,p,ci,stats] = ttest(space_volt,one_volt,'alpha',0.05,'tail','both');
% fprintf('Space vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

% [h,p,ci,stats] = ttest(mass_volt,one_volt,'alpha',0.05,'tail','both');
% fprintf('Mass vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);


%% voltage interaction

fprintf('space - one: %.4f uV\n',mean(space_volt - one_volt));
fprintf('mass - one: %.4f uV\n',mean(mass_volt - one_volt));

[h,p,ci,stats] = ttest([space_volt - one_volt],[mass_volt - one_volt],'alpha',0.05,'tail','both');
fprintf('Space/one vs mass/one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);


%% correlation - loading

subDir = '';
behDir = fullfile(exper.name,'Behavioral','Sessions',subDir);
behDir = fullfile(dirs.dataroot,behDir);

collapsePhases = true;
if collapsePhases
  collapseStr = '_collapsed';
else
  collapseStr = '';
end

% split into quantile divisions?
nDivisions = 1;
% nDivisions = 2;
% nDivisions = 3;
% nDivisions = 4;

if nDivisions > 1
  quantStr = sprintf('_%dquantileDiv',nDivisions);
else
  quantStr = '';
end

% load the behavioral data
%resultsFile = fullfile(dataroot,dirs.behDir,sprintf('%s_behav_results%s%s.mat',expName,quantStr,collapseStr));
resultsFile = fullfile(behDir,sprintf('%s_behav_results%s%s.mat',exper.name,quantStr,collapseStr));

fprintf('Loading %s...',resultsFile);
load(resultsFile);
fprintf('Done.\n');

%% correlation

ses = 'oneDay';
phase = 'cued_recall';

test = 'recall';
% measure = 'recall_nHit';
measure = 'recall_hr';

% test = 'recog';
% % measure = 'recall_nHit';
% measure = 'recog_dp';

manip = 'spaced';
space_measure = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);
manip = 'massed';
mass_measure = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);

% [rho,p] = corr(data1,space_volt);
% fprintf('rho=%.4f, p=%.5f\n',rho,p);

[rho,p] = corr(space_measure,space_mem_volt);
fprintf('space mem: rho=%.4f, p=%.5f\n',rho,p);

[rho,p] = corr(mass_measure,mass_mem_volt);
fprintf('mass mem: rho=%.4f, p=%.5f\n',rho,p);

% [rho,p] = corr(data1,space_volt - one_volt);

% [rho,p] = corr(data1,mass_volt - one_volt);

spacEff_rc_measure = results.(ses).(phase).spaced.(test).(measure)(~exper.badSub) - results.(ses).(phase).massed.(test).(measure)(~exper.badSub);
[rho,p] = corr(spacEff_rc_measure,spacEff_rc_volt);
fprintf('spacing effect rc: rho=%.4f, p=%.5f\n',rho,p);

spacEff_rc_measure = results.(ses).(phase).spaced.(test).(measure)(~exper.badSub) - ((results.(ses).(phase).spaced.(test).(measure)(~exper.badSub) + results.(ses).(phase).massed.(test).(measure)(~exper.badSub)) ./ 2);
[rho,p] = corr(spacEff_rc_measure,spac_repEff_rc_volt);
fprintf('spac: rep effect corr w spac minus avg recall: rho=%.4f, p=%.5f\n',rho,p);

massEff_rc_measure = results.(ses).(phase).massed.(test).(measure)(~exper.badSub) - ((results.(ses).(phase).spaced.(test).(measure)(~exper.badSub) + results.(ses).(phase).massed.(test).(measure)(~exper.badSub)) ./ 2);
[rho,p] = corr(massEff_rc_measure,mass_repEff_rc_volt);
fprintf('mass: rep effect corr w mass minus avg recall: rho=%.4f, p=%.5f\n',rho,p);
