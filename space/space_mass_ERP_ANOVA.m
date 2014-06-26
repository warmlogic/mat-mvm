
%% initialize

allPeakInfo = struct;

cfg = [];

%% gather data

% spacings = {'mass', 'spac', 'onePres'};
% oldnew = {'p2'};
% memConds = {'all'};

% I didn't test new words, so they can't be recalled/forgotten, but we can
% use p1
spacings = {'mass', 'spac'};
oldnew = {'p1', 'p2'};
memConds = {'rc','fo'};

erpComponents = {'LPC','N400'};

for sp = 1:length(spacings)
  
  for on = 1:length(oldnew)
  
  for mc = 1:length(memConds)
    
    if strcmp(spacings{sp},'onePres')
      % % single presentation or first presentation
      if strcmp(memConds{mc},'all');
        cfg.conditions = {'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1'};
        %cfg.conditions = {'word_onePres'};
      end
    elseif strcmp(spacings{sp},'mass') || strcmp(spacings{sp},'spac')
      
      if strcmp(memConds{mc},'all');
        cfg.conditions = {sprintf('word_RgH_rc_%s_%s',spacings{sp},oldnew{on}),sprintf('word_RgH_fo_%s_%s',spacings{sp},oldnew{on})};
      else
        cfg.conditions = {sprintf('word_RgH_%s_%s_%s',memConds{mc},spacings{sp},oldnew{on})};
        
      end
      
    end
    
    disp(cfg.conditions);
      
%     if strcmp(spacings{sp},'spaced')
%       
%       spac_str = 'spac';
%       
%       % spaced
%       if strcmp(memConds{mc},'all');
%         cfg.conditions = {'word_RgH_rc_spac_p2','word_RgH_fo_spac_p2'};
%       elseif strcmp(memConds{mc},'recalled');
%         cfg.conditions = {'word_RgH_rc_spac_p2'};
%       elseif strcmp(memConds{mc},'forgot');
%         cfg.conditions = {'word_RgH_fo_spac_p2'};
%       end
%       
%     elseif strcmp(spacings{sp},'massed')
%       % % massed
%       if strcmp(memConds{mc},'all');
%         cfg.conditions = {'word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};
%       elseif strcmp(memConds{mc},'recalled');
%         cfg.conditions = {'word_RgH_rc_mass_p2'};
%       elseif strcmp(memConds{mc},'forgot');
%         cfg.conditions = {'word_RgH_fo_mass_p2'};
%       end
%       
%     elseif strcmp(spacings{sp},'once')
%       % % single presentation or first presentation
%       if strcmp(memConds{mc},'all');
%         cfg.conditions = {'word_onePres','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1'};
%         %cfg.conditions = {'word_onePres'};
%       end
%     end
    
    for er = 1:length(erpComponents)
      if strcmp(erpComponents{er},'LPC')
        % LPC
        cfg.order = 'descend'; % descend = positive peaks first
        cfg.roi = {'Pz'};
        % cfg.latency = [0.4 0.8];
        lpcPeak = 0.592;
        % cfg.latency = [lpcPeak-0.05 lpcPeak+0.05]; % LPC - around GA peak (space+mass) +/- 50
        cfg.latency = [lpcPeak-0.1 lpcPeak+0.1]; % LPC - around GA peak (space+mass) +/- 100
      elseif strcmp(erpComponents{er},'N400')
        % N400
        cfg.order = 'ascend'; % ascend = negative peaks first
        cfg.roi = {'Cz'};
        % cfg.latency = [0.2 0.6];
        n400Peak = 0.364;
        % % cfg.latency = [n400Peak-0.05 n400Peak+0.05]; % N400 - around GA peak (space+mass) +/- 50
        cfg.latency = [n400Peak-0.1 n400Peak+0.1]; % N400 - around GA peak (space+mass) +/- 100
      end
      
      % % average across time window
      %cfg.datadim = 'elec';
      % cfg.roi = {'center101'};
      % % cfg.roi = {'PS2'};
      % % cfg.roi = {'LPI3','RPI3'};
      % cfg.latency = [0.4 0.8]; % LPC
      % % % cfg.latency = [0.3 0.5]; % N400
      % % % cfg.latency = [0.35 0.45]; % N400
      % % cfg.latency = [0.314 0.414]; % N400
      
      % % average across electrodes
      cfg.datadim = 'time';
      % and time points
      cfg.avgovertime = true;
      
      % % cfg.roi = {'Cz'};
      % % cfg.roi = {'LPI3','RPI3'};
      % % cfg.roi = {'Pz'};
      % % cfg.roi = {'PS2'};
      % % cfg.roi = {'RPI3'};
      % % cfg.roi = {'E84'}; % center of RPI3
      % % cfg.roi = {'RPS2'};
      % % cfg.roi = {'E85'}; % center of RPS2
      % % cfg.roi = {'LPS2'};
      % % cfg.latency = [0 1.0];
      % % cfg.latency = [0.2 0.9];
      
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
      
      allPeakInfo.(sprintf('%s_%s_%s',spacings{sp},oldnew{on},memConds{mc})).(erpComponents{er}) = peakInfo;
      
    end
  end
  end
end

%% ANOVA: factors: spaced/massed, recalled/forgotten, old/new???

spacings = {'massed', 'spaced', 'once'};
memConds = {'all'};

% I didn't test new words, so they can't be recalled/forgotten
% spacings = {'massed', 'spaced'};
% memConds = {'recalled','forgot'};

measure = 'latency';
% measure = 'voltage';

erpComp = 'LPC';
% erpComp = 'N400';

anovaData = [];

for sub = 1:sum(~exper.badSub)
  theseData = [];
  
  for sp = 1:length(spacings)
    for mc = 1:length(memConds)
      theseData = cat(2,theseData,allPeakInfo.(sprintf('%s_%s',spacings{sp},memConds{mc})).(erpComp).subjects.(measure)(sub,1));
    end
  end
  anovaData = cat(1,anovaData,theseData);
end

fprintf('================================================================\n');
fprintf('This test: %s %s\n',erpComp,measure);

if length(memConds) > 1
  levelnames = {spacings memConds};
  varnames = {'spacing', 'memory'};
  O = teg_repeated_measures_ANOVA(anovaData, [length(spacings) length(memConds)], varnames,[],[],[],[],[],[],levelnames);
else
  levelnames = {spacings};
  varnames = {'spacing'};
  O = teg_repeated_measures_ANOVA(anovaData, [length(spacings)], varnames,[],[],[],[],[],[],levelnames);
end

fprintf('Prev test: %s %s\n',erpComp,measure);
fprintf('================================================================\n');

%% gather data for pairwise t-tests

% erpComp = 'LPC';
% erpComp = 'N400';

% erpComp = 'lpc';
erpComp = 'n400';

fprintf('%s\n',erpComp);

space_peak = allPeakInfo.spaced_all.(erpComp).subjects.latency(:,1);
space_volt = allPeakInfo.spaced_all.(erpComp).subjects.voltage(:,1);

mass_peak = allPeakInfo.massed_all.(erpComp).subjects.latency(:,1);
mass_volt = allPeakInfo.massed_all.(erpComp).subjects.voltage(:,1);

one_peak = allPeakInfo.once_all.(erpComp).subjects.latency(:,1);
one_volt = allPeakInfo.once_all.(erpComp).subjects.voltage(:,1);

%% ttest - latency

fprintf('space: %.4f sec\n',mean(space_peak));
fprintf('mass: %.4f sec\n',mean(mass_peak));
fprintf('one: %.4f sec\n',mean(one_peak));

[h,p,ci,stats] = ttest(space_peak,mass_peak,'alpha',0.05,'tail','both');
fprintf('Space vs mass: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(space_peak,one_peak,'alpha',0.05,'tail','both');
fprintf('Space vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(mass_peak,one_peak,'alpha',0.05,'tail','both');
fprintf('Mass vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

%% ttest - voltage

fprintf('space: %.4f uV\n',mean(space_volt));
fprintf('mass: %.4f uV\n',mean(mass_volt));
fprintf('one: %.4f uV\n',mean(one_volt));

[h,p,ci,stats] = ttest(space_volt,mass_volt,'alpha',0.05,'tail','both');
fprintf('Space vs mass: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(space_volt,one_volt,'alpha',0.05,'tail','both');
fprintf('Space vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(mass_volt,one_volt,'alpha',0.05,'tail','both');
fprintf('Mass vs one: t(%d)=%.4f, p=%.8f\n',stats.df,stats.tstat,p);


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
measure = 'recall_nHit';
% measure = 'recall_hr';

manip = 'massed';
data1 = results.(ses).(phase).(manip).(test).(measure)(~exper.badSub);

[rho,p] = corr(data1,space_volt - one_volt);

% [rho,p] = corr(data1,mass_volt - one_volt);
