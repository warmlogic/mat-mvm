function space_pow_rmanova(freqs,latencies,sigElecs,stimType,ana,exper,data_pow)

% stimType = 'word_';
% stimType = 'img_';
memType = 'RgH_';

% spacings = {'mass', 'spac', 'onePres'};
% % oldnew = {'p1'};
% oldnew = {'p2'};
% memConds = {'all'};

spacings = {'mass', 'spac'};
% spacings = {'spac'};
% oldnew = {'p1', 'p2'};
% oldnew = {'p1'};
oldnew = {'p2'};
memConds = {'rc','fo'};
% memConds = {'rc'};

measure = 'powspctrm';

% latencies = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
% latencies = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap

% latencies = [0.02:0.1:0.92; 0.1:0.1:1.0]'; % 100 no overlap
% latencies = [0.02:0.2:0.92; 0.2:0.2:1.0]'; % 200 no overlap
% latencies = [0.02:0.25:0.92; 0.25:0.25:1.0]'; % 250 no overlap
% latencies = [0.02:0.32:0.92; 0.32:0.32:1.0]'; % 300 no overlap
% latencies = [0.02:0.5:0.92; 0.5:0.5:1.0]'; % 500 no overlap

roi = {sigElecs};
% freqs = ana.freq.theta;
% % freqs = ana.freq.alpha;
% freqs = ana.freq.alpha_lower;
% freqs = ana.freq.alpha_upper;
% freqs = ana.freq.beta_lower;


% % % theta
% freqs = ana.freq.theta;
% roi = {sigElecs};
% 
% % latencies = [0.6 1.0]; % word
% % latencies = [0.1 0.4]; % img
% % latencies = [0.02 0.5; 0.52 1.0];
% 
% % img
% % roi = {'LAI'}; % yes, neg **
% % roi = {'LFP'};
% % roi = {'FC'}; % pos
% % roi = {'RFP'};
% % roi = {'RAI'}; % something awry word mass p1 forgot, values are too high
% 
% % roi = {{'C','FS'}}; % 
% % roi = {{'LAS2','C','FS','LPS'}}; % 
% % roi = {'anterior_noPeriph'};
% % roi = {{'LAS2','FS'}}; % 
% % roi = {'LAS2'}; % yes, pos ***
% % roi = {'LAS'}; % yes, pos **
% % roi = {'FS'}; % yes, pos **
% % roi = {{'FS','LAS2'},{'LPS','PS'}}; % yes, pos **
% % roi = {{'FC','LFP'},{'LPS','PS'}}; % yes, pos **
% % roi = {'C'}; % yes, pos
% % roi = {'RAS'}; % yes, pos *
% % roi = {'RAS2'}; % yes
% 
% % roi = {{'LT','FS'}}; % 
% % roi = {'LT'}; % 
% % roi = {'LPS2'}; % 
% % roi = {'LPS'}; % 
% % roi = {'PS'}; % 
% % roi = {'RPS'}; % 
% % roi = {'RPS2'}; % 
% % roi = {'RT'}; % 
% 
% % roi = {'LPI2'}; % yes
% % roi = {'PI'}; % no
% % roi = {'RPI2'}; % yes
% 
% % % word
% % roi = {'LAI'}; % yes, neg **
% % roi = {'LFP'};
% % roi = {'FC'};
% % roi = {'RFP'};
% % roi = {'RAI'}; % something awry word mass p1 forgot, values are too high
% 
% % roi = {{'LAS2','FS'}}; % 
% % roi = {'LAS2'}; % yes, pos ***
% % roi = {'LAS'}; % yes, pos **
% % roi = {'FS'}; % yes, pos **
% % roi = {'C'}; % yes, pos
% % roi = {'RAS'}; % yes, pos *
% % roi = {'RAS2'}; % yes
% 
% % roi = {{'LT','FS'}}; % 
% % roi = {'LT'}; % 
% % roi = {'LPS2'}; % 
% % roi = {'LPS'}; % 
% % roi = {'PS'}; % 
% % roi = {'RPS'}; % 
% % roi = {'RPS2'}; % 
% % roi = {'RT'}; % 
% 
% % roi = {'LPI2'}; % yes
% % roi = {'PI'}; % no
% % roi = {'RPI2'}; % yes
% 
% % roi = {{'E23','Fz','E3'}}; % AF3 Fz AF4


% % % alpha
% % freqs = ana.freq.alpha;
% freqs = ana.freq.alpha_lower;
% % freqs = ana.freq.alpha_upper;
% % freqs = ana.freq.beta_lower;
% roi = {sigElecs};

% % freqs = [11 20.5];
% % latencies = [0.6 1.0]; % img
% % roi = {'LAS2'};
% 
% % roi = {{'E67','Pz','E77'}}; % PO3 Pz P04
% 
% % latencies = [0.1 0.3]; % word, LT, early alpha effect Spac x Mem
% % latencies = [0.4 0.7]; % 
% % latencies = [0.5 0.9]; % 
% % latencies = [0.4 1.0]; % 
% % roi = {'LT'}; % word **
% 
% % latencies = [0.3 0.7]; % word
% % latencies = [0.5 0.7]; % word
% % latencies = [0.8 1.0]; % word
% % roi = {'LAS2'}; %
% % roi = {{'LAS2','LT','LPS','C','PS','RPS','RT','RAS2'}}; % 
% % roi = {{'LT','LPS2','PS','PI','RPS2','RT'}}; % 
% % roi = {{'LT','PS','RT'}}; % 
% % roi = {{'LPI2','PI','RPI2'}}; % 
% % roi = {{'LPS','C','PS','RPS'}}; % 
% % roi = {'PS'}; % word **
% % roi = {'LPS'}; % 
% % roi = {'PI'}; % word ** 
% % roi = {'posterior_noPeriph'};
% 
% % latencies = [0.4 0.6];
% % latencies = [0.4 0.7];
% % roi = {'FC'}; % no
% % roi = {'FS'}; % maybe
% % roi = {'FS2'}; % yes
% % roi = {'PS'}; %
% % roi = {'LPI2'}; %

latency = cell(1,size(latencies,1));
for i = 1:length(latency)
  latency{i} = sprintf('%dto%d',round(latencies(i,1)*1000),round(latencies(i,2)*1000));
end
latStr = sprintf(repmat('_%s',1,length(latency)),latency{:});
latStr = latStr(2:end);

freq_str = sprintf('%dfreq%dto%d',size(freqs,1),round(freqs(1,1)),round(freqs(end,end)));
freqIdx = (nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,freqs(1,1)):nearest(data_pow.(exper.sesStr{1}).(ana.eventValues{1}{1}{1}).sub(1).data.freq,freqs(1,2)));

factorNames = {'spacings', 'oldnew', 'memConds', 'roi', 'latency'};

nVariables = nan(size(factorNames));
keepTheseFactors = false(size(factorNames));
levelNames_teg = cell(size(factorNames)); % TEG
for c = 1:length(factorNames)
  % need to have a variable set to this exact factor name
  thisFac = eval(factorNames{c});
  nVariables(c) = length(thisFac);
  if nVariables(c) > 1
    keepTheseFactors(c) = true;
  end
  
  % go through the levels within this factor and save strings
  thisFac_str = cell(1,nVariables(c));
  for l = 1:length(thisFac)
    if iscell(thisFac{l})
      tf = thisFac{l};
      thisFac_str{l} = sprintf(repmat('_%s',1,length(tf)),tf{:});
      thisFac_str{l} = thisFac_str{l}(2:end);
    elseif ischar(thisFac{l})
      thisFac_str{l} = thisFac{l};
    end
  end
  levelNames_teg{c} = thisFac_str; % TEG
end

variableNames = cell(1,prod(nVariables));
levelNames = cell(prod(nVariables),length(factorNames));

% ses=1;
nSub = sum(~exper.badSub);
anovaData = nan(nSub,prod(nVariables));
rmaov_data_teg = nan(nSub*prod(nVariables),length(factorNames) + 2);

fprintf('Collecting ANOVA data for %d subjects:\n\t',nSub);
fprintf('%s\n',stimType);
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(spacings)),spacings{:}),factorNames{1});
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(oldnew)),oldnew{:}),factorNames{2});
fprintf('%s (%s),',sprintf(repmat(' %s',1,length(memConds)),memConds{:}),factorNames{3});
if iscell(roi{1})
  fprintf('%d ROIs (%s),',length(roi),factorNames{4});
elseif ischar(roi{1})
  fprintf('%s (%s),',sprintf(repmat(' %s',1,length(roi)),roi{:}),factorNames{4});
end
fprintf('%s (%s),',latStr,factorNames{5});
fprintf('\n\tFreq: %s...',freq_str);

lnDone = false;
vnDone = false;
subCount = 0;
rmCount = 0;
for sub = 1:length(exper.subjects)
  if ~exper.badSub(sub)
    subCount = subCount + 1;
  else
    continue
  end
  for ses = 1:length(exper.sesStr)
    lnCount = 0;
    vnCount = 0;
    
    for sp = 1:length(spacings)
      for on = 1:length(oldnew)
        for mc = 1:length(memConds)
          cond_str = [];
          cond_str_tmp = [];
          if strcmp(spacings{sp},'onePres')
            % single presentation or first presentation
            %cond_str = sprintf('%s%s%s_%s',stimType,memType,memConds{mc},spacings{sp});
            if strcmp(memConds{mc},'all')
              cond_str = sprintf('%s%s',stimType,spacings{sp});
            %  cond_str_tmp = {sprintf('%s%s%s_%s',stimType,memType,'rc',spacings{sp}), sprintf('%s%s%s_%s',stimType,memType,'fo',spacings{sp})};
            end
          elseif strcmp(spacings{sp},'mass') || strcmp(spacings{sp}(1:4),'spac')
            cond_str = sprintf('%s%s%s_%s_%s',stimType,memType,memConds{mc},spacings{sp},oldnew{on});
            if strcmp(memConds{mc},'all')
              % manually input rc and fo
              cond_str_tmp = {sprintf('%s%s%s_%s_%s',stimType,memType,'rc',spacings{sp},oldnew{on}), sprintf('%s%s%s_%s_%s',stimType,memType,'fo',spacings{sp},oldnew{on})};
            end
          end
          
          for r = 1:length(roi)
            if iscell(roi{r})
              roi_str = sprintf(repmat('%s',1,length(roi{r})),roi{r}{:});
            elseif ischar(roi{r})
              roi_str = roi{r};
            end
            if ~isempty(cond_str_tmp)
              chanIdx = ismember(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,roi{r})})));
            else
              chanIdx = ismember(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.label,unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,roi{r})})));
            end
            
            for lat = 1:length(latency)
              if ~isempty(cond_str_tmp)
                tbeg = nearest(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.time,latencies(lat,1));
                tend = nearest(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.time,latencies(lat,2));
              else
                tbeg = nearest(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,1));
                tend = nearest(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.time,latencies(lat,2));
              end
              latIdx = tbeg:tend;
              
              if ~lnDone
                lnCount = lnCount + 1;
                levelNames{lnCount,1} = spacings{sp};
                levelNames{lnCount,2} = oldnew{on};
                levelNames{lnCount,3} = memConds{mc};
                levelNames{lnCount,4} = roi{r};
                levelNames{lnCount,5} = latency{lat};
              end
              
              vnCount = vnCount + 1;
              if ~vnDone
                variableNames{vnCount} = sprintf('Y%d',vnCount);
              end
              
              if ~isempty(cond_str_tmp)
                thisData = (mean(mean(mean(data_pow.(exper.sesStr{ses}).(cond_str_tmp{1}).sub(sub).data.(measure)(chanIdx,freqIdx,latIdx),3),2),1) + ...
                  mean(mean(mean(data_pow.(exper.sesStr{ses}).(cond_str_tmp{2}).sub(sub).data.(measure)(chanIdx,freqIdx,latIdx),3),2),1)) ./ 2;
              else
                thisData = mean(mean(mean(data_pow.(exper.sesStr{ses}).(cond_str).sub(sub).data.(measure)(chanIdx,freqIdx,latIdx),3),2),1);
              end
              anovaData(subCount,vnCount) = thisData;

              rmCount = rmCount + 1;
              rmaov_data_teg(rmCount,:) = [anovaData(subCount,vnCount) sp on mc r lat sub];
            end
          end
        end
      end
    end
    lnDone = true;
    vnDone = true;
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
fprintf('This ANOVA: %s, Freq: %s\n\n',stimType,freq_str);

O = teg_repeated_measures_ANOVA(anovaData, nVariables, factorNames,[],[],[],[],[],[],levelNames_teg,rmaov_data_teg);

fprintf('Prev ANOVA: %s, Freq: %s\n',stimType,freq_str);
fprintf('=======================================\n');

end
