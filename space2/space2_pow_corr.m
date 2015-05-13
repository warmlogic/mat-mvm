

cfg_ana = [];
cfg_ana.latencies = [-0.18:0.1:0.92; -0.1:0.1:1.0]';

%% RM ANOVA: run through everything we need from cluster statistics in one go

% NB: need to run "cluster statistics" cell (without executing script)

%%%%%%%%%%%%%%%%
% GENERAL
%%%%%%%%%%%%%%%%

makePlots = false;
% makePlots = true;

% theseFreqs = cfg_ana.frequencies;
% theseFreqs = [ ...
%   ana.freq.theta; ...
%   %ana.freq.alpha; ...
%   ana.freq.alpha_lower; ...
%   ana.freq.alpha_upper; ...
%   ana.freq.beta_lower; ...
%   %ana.freq.beta_upper; ...
%   %ana.freq.gamma_lower; ...
%   %ana.freq.gamma_upper; ...
%   ];
theseFreqs = [ana.freq.theta];
% theseFreqs = [ana.freq.alpha_lower];
% theseFreqs = [ana.freq.alpha_upper];
% theseFreqs = [ana.freq.beta_lower];

% stimTypes = {'word_', 'img_'};
stimTypes = {'word_'};
% stimTypes = {'img_'};

badSub = {'SPACE2002', 'SPACE2007', 'SPACE2025', 'SPACE2016', 'SPACE2034'};

% behCats = {'once', 'massed', 'lag2', 'lag12', 'lag32'};
% eegCats = {'onePres', 'mass', 'spac2', 'spac12', 'spac32'};

behCats = {'massed', 'lag2', 'lag12', 'lag32'};
eegCats = {'mass', 'spac2', 'spac12', 'spac32'};

behSess = {'day1', 'day2'};
eegSess = 'session_1_session_2';

performMeasure = 'recall_hr';
eegMeasure = 'powspctrm';

% NB: set subjects_beh to the list of subjects from space2_behavioral (len=32)      

%% run it

for st = 1:length(stimTypes)
  stimType = stimTypes{st};
  
%   files.figPrintFormat = 'png';
  
  %%%%%%%%%%%%%%%%
  % FIND SIGNIFICANT ELECTRODES
  %%%%%%%%%%%%%%%%
  
  if strcmp(stimType,'word_')
    sigElecConditions = {{'word_rc_mass_p2','word_fo_mass_p2','word_rc_spac2_p2','word_fo_spac2_p2','word_rc_spac12_p2','word_fo_spac12_p2','word_rc_spac32_p2','word_fo_spac32_p2'}};
  elseif strcmp(stimType,'img_')
    sigElecConditions = {{'img_rc_mass_p2','img_fo_mass_p2','img_rc_spac2_p2','img_fo_spac2_p2','img_rc_spac12_p2','img_fo_spac12_p2','img_rc_spac32_p2','img_fo_spac32_p2'}};
  end
  
  % significantly different in half of comparisons
  nSigComparisons = floor(size(nchoosek(sigElecConditions{1},2),1) / 2);
  % % significantly different in one fourth of comparisons
  % nSigComparisons = floor(size(sigElecsAcrossComparisons.pos{1},2) / 4);
  % % significantly different in 5 of 28 of comparisons
  % nSigComparisons = floor(size(sigElecsAcrossComparisons.pos{1},2) / 5);
  
  %%%%%%%%%%%%%%%%
  % CONTRAST TOPO
  %%%%%%%%%%%%%%%%
  
  % TODO: more condition comparisons, see cluster stats above
  
%   if strcmp(stimType,'word_')
%     contrastConditions = {...
%       {'word_rc_mass_p2', 'word_rc_onePres'} ... % P2 Rc Spacing
%       {'word_rc_mass_p2', 'word_rc_spac2_p2'} ... % P2 Rc Spacing
%       {'word_rc_mass_p2', 'word_rc_spac12_p2'} ... % P2 Rc Spacing
%       {'word_rc_mass_p2', 'word_rc_spac32_p2'} ... % P2 Rc Spacing
%       {'word_fo_mass_p2', 'word_fo_onePres'} ... % P2 Fo Spacing
%       {'word_fo_mass_p2', 'word_fo_spac2_p2'} ... % P2 Fo Spacing
%       {'word_fo_mass_p2', 'word_fo_spac12_p2'} ... % P2 Fo Spacing
%       {'word_fo_mass_p2', 'word_fo_spac32_p2'} ... % P2 Fo Spacing
%       };
%   elseif strcmp(stimType,'img_')
%     contrastConditions = {...
%       {'img_rc_mass_p2', 'img_rc_onePres'} ... % P2 Rc Spacing
%       {'img_rc_mass_p2', 'img_rc_spac2_p2'} ... % P2 Rc Spacing
%       {'img_rc_mass_p2', 'img_rc_spac12_p2'} ... % P2 Rc Spacing
%       {'img_rc_mass_p2', 'img_rc_spac32_p2'} ... % P2 Rc Spacing
%       {'img_fo_mass_p2', 'img_fo_onePres'} ... % P2 Fo Spacing
%       {'img_fo_mass_p2', 'img_fo_spac2_p2'} ... % P2 Fo Spacing
%       {'img_fo_mass_p2', 'img_fo_spac12_p2'} ... % P2 Fo Spacing
%       {'img_fo_mass_p2', 'img_fo_spac32_p2'} ... % P2 Fo Spacing
%       };
%   end
  
  % topoLatencies = [0 0.5; 0.52 1.0]; % time
  
  %%%%%%%%%%%%%%%%
  % SMOOTH LINE PLOTS
  %%%%%%%%%%%%%%%%
  
%   linePlotLatencies = [-0.1:0.02:0.98; -0.08:0.02:1.0]';
%   
%   if strcmp(stimType,'word_')
%     linePlotConditions = {{{'word_rc_mass_p2', 'word_fo_mass_p2', 'word_rc_spac2_p2', 'word_fo_spac2_p2'}},{{'word_rc_spac12_p2', 'word_fo_spac12_p2', 'word_rc_spac32_p2', 'word_fo_spac32_p2'}}};
%   elseif strcmp(stimType,'img_')
%     linePlotConditions = {{{'img_rc_mass_p2', 'img_fo_mass_p2', 'img_rc_spac2_p2', 'img_fo_spac2_p2'}},{{'img_rc_spac12_p2', 'img_fo_spac12_p2', 'img_rc_spac32_p2', 'img_fo_spac32_p2'}}};
%   end
%   linePlotConditions_rename = {{{'Mass P2 Recall','Mass P2 Forgot','Space2 P2 Recalled','Space2 P2 Forgot'}},{{'Space12 P2 Recalled','Space12 P2 Forgot','Space32 P2 Recalled','Space32 P2 Forgot'}}};
  
%   linePlotLinestyle = {{'-','--','-','--'},{'-','--','-','--'}};
%   if exist('linspecer','file')
%     thisColor = linspecer(length(cellflat(linePlotConditions)));
%     linePlotColors = cell(1,2);
%     linePlotColors{1} = thisColor(1:4,:);
%     linePlotColors{2} = thisColor(5:8,:);
%   else
%     linePlotColors = {'bcrm','bcrm'};
%   end
  
  % onePres_colors = [[0 0 0]; [0 1 0]];
  
  %%%%%%%%%%%%%%%%
  % AVERAGE PLOTS
  %%%%%%%%%%%%%%%%
  
%   avgPlotConds = sigElecConditions{1};
%   avgPlotConds_rename = {'Mass P2 Recall','Mass P2 Forgot','Space2 P2 Recalled','Space2 P2 Forgot','Space12 P2 Recalled','Space12 P2 Forgot','Space32 P2 Recalled','Space32 P2 Forgot'};
  
  %%%%%%%%%%%%%%%%
  % RUN IT
  %%%%%%%%%%%%%%%%
  
  for f = 1:size(theseFreqs,1)
    % find the significant electodes
    files.saveFigs = false;
    sigElecs = space2_pow_sigElecs(sigElecConditions,theseFreqs(f,:),cfg_ana.latencies,cfg_ana.latencies,nSigComparisons,ana,exper,files,dirs,ga_pow);
    close all
    
%     if makePlots
%       % make the smooth line plots
%       files.saveFigs = true;
%       for lpc = 1:size(linePlotConditions,2)
%         space2_pow_linePlot(linePlotConditions{lpc},linePlotConditions_rename{lpc},theseFreqs(f,:),linePlotLatencies,sigElecs,linePlotColors{lpc},linePlotLinestyle{lpc},ana,exper,files,dirs,ga_pow);
%       end
%       close all
%     end
    
%     % set latencies for average plots, topoplots, and ANOVA
%     if (theseFreqs(f,1) == 11 && theseFreqs(f,2) == 12) || (theseFreqs(f,1) == 13.1 && theseFreqs(f,2) == 20.5)
%       latencies = [0.02:0.32:0.92; 0.32:0.32:1.0]'; % 300 no overlap
% %       avgPlotLatencies = [0 0.333; 0.333 0.666; 0.666 1.0];
%     else
%       latencies = [0.02:0.5:0.92; 0.5:0.5:1.0]'; % 500 no overlap
% %       avgPlotLatencies = [0 0.5; 0.5 1.0];
%     end
    
    latencies = [0.52, 1.0];
    lat = 1;
    
%     if makePlots
%       % make the average plots
%       files.saveFigs = true;
%       space2_pow_avgPlot(avgPlotConds,avgPlotConds_rename,theseFreqs(f,:),avgPlotLatencies,sigElecs,ana,exper,files,dirs,data_pow);
%       close all
%     end
    
%     % make the contrast topoplots
%     if makePlots
%       for t = 1:size(latencies,1)
%         files.saveFigs = true;
%         space2_pow_contrastTopo(contrastConditions,theseFreqs(f,:),latencies(t,:),sigElecs,ana,exper,files,dirs,ga_pow);
%         close all
%       end
%     end
    
    % run the RM ANOVA
    %space2_pow_rmanova(theseFreqs(f,:),latencies,sigElecs,stimType,ana,exper,data_pow);
    
    %roi = {sigElecs};
    chanIdx = ismember(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,sigElecs)})));

    freqIdx = (nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,theseFreqs(f,1)):nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,theseFreqs(f,2)));
    
    tbeg = nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.time,latencies(lat,1));
    tend = nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.time,latencies(lat,2));
    timeIdx = tbeg:tend;
    
    behData = nan(sum(~exper.badSub),length(behCats));
    eegData = nan(sum(~exper.badSub),length(eegCats));
    
    subCount = 0;
    
    for sub = 1:length(exper.subjects)
      if ~ismember(exper.subjects{sub},badSub)
        subCount = subCount + 1;
        behInd = find(ismember(subjects_beh,exper.subjects{sub}));
        for c = 1:length(behCats)
          %behData(subCount,c) = results.(behSess{1}).cued_recall_only.(behCats{c}).recall.(performMeasure)(behInd);
          behData(subCount,c) = (results.(behSess{1}).cued_recall_only.(behCats{c}).recall.(performMeasure)(behInd) +...
            results.(behSess{2}).cued_recall_only.(behCats{c}).recall.(performMeasure)(behInd)) / 2;
          
          if ~strcmp(eegCats{c},'onePres')
            eegCat = sprintf('%src_%s_p2',stimType,eegCats{c});
          else
            eegCat = sprintf('%src_%s',stimType,eegCats{c});
          end
          eegData(subCount,c) = mean(mean(mean(data_pow.(exper.sesStr{1}).(eegCat).sub(sub).data.(eegMeasure)(chanIdx,freqIdx,timeIdx),3),2),1);
          
        end
      end
    end
    
    % correlation - not sure this is how to run it
    [rho,p] = corr(behData,eegData);

  end

end



%%

roi = {'sigElecs'};

spacings = {'mass', 'spac2', 'spac12', 'spac32'};
% oldnew = {'p1', 'p2'};
oldnew = {'p2'};
% memConds = {'rc','fo'};
memConds = {'rc'};

% measure = 'latency';
measure = 'voltage';

if strcmp(measure,'latency')
  %nPoints = cfg.nPoints_sub;
  nPoints = 10;
%   nPoints = 1;
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
% roi = {'E62_E72_E76_E77_E78_E84_E85'}; % centered on E77

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

%%

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
        
        cond_str = sprintf('%s_%s_%s',spacings{sp},oldnew{on},memConds{mc});
        
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
          
          %subData = mean(allPeakInfo.(cond_str).(roi_str).subjects.(measure)(sub,1:nPoints));
          %anovaData(sub,vnCount) = subData;
          
          %rmaov_data_teg = cat(1,rmaov_data_teg,[subData, sp, on, mc, r, sub]);
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
  %levelNames_teg = levelNames_teg(keepTheseFactors); % TEG
  
  %rmaov_data_teg = rmaov_data_teg(:,[1 (find(keepTheseFactors) + 1) size(rmaov_data_teg,2)]);
end

fprintf('Done.\n');


%% output

fid = fopen('~/Desktop/SPACE2_eeg_recall.csv','w');

% for sub = 1:size(behData,1),
%   [rho, p] = corr(behData(sub,:)', eegData(sub,:)');
%   fprintf('%.4f, %.4f\n',rho,p);
%   allrho=cat(1,allrho,rho);allp=cat(1,allp,p);
% end

fprintf(fid,'sub,space,theta,recall\n');

for sub = 1:size(behData,1)
  for sp = 1:length(spacings)
    fprintf(fid,'%d',sub);
    fprintf(fid,',%s',spacings{sp});
    fprintf(fid,',%.6f,%.6f\n',eegData(sub,sp),behData(sub,sp));
  end
end

fclose(fid);

%% input

% nSub = 27;
% nCond = 4;
% 
% eegData2 = [];
% behData2 = [];

fid = fopen('~/Desktop/SPACE2_eeg_recall.csv','r');
C = textscan(fid,'%d%s%f%f','Delimiter',',','Headerlines',1);
fclose(fid);

%%

% t = cell2table(C,'VariableNames',{'subjects','spacings','y1','y2'});
t = table(C{1},C{2},C{3},C{4},'VariableNames',{'subject','spacings','y1','y2'});

within = table([1 2]','VariableNames',{'Measurements'});

% within = cell2table({'mass', 'spac2', 'spac12', 'spac32'}','VariableNames',{'spacings'});
% within = cell2table({'mass', 'spac2', 'spac12', 'spac32'}','VariableNames',{'spacings'});

% Fit a repeated measures model where the measurements are theta and recall and spacings is the predictor variable.

rm = fitrm(t,'y1-y2~spacings','WithinDesign',within);

manova(rm)

manova(rm,'By','spacings')

