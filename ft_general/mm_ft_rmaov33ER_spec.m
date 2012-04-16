function mm_ft_rmaov33ER_spec(cfg_ana,exper,ana,data)
%MM_FT_RMAOV33ER_SPEC 3-way RMANOVA; 2- & 1-way RMANOVAs collapsing across IVs
%
% Can specify the setup
%

if ~isfield(cfg_ana,'parameter')
  error('Must specify cfg_ana.parameter, denoting the data to test (e.g., ''avg'' or ''individual'')');
end

if ~isfield(cfg_ana,'alpha')
  cfg_ana.alpha = 0.05;
end
if ~isfield(cfg_ana,'showtable')
  cfg_ana.showtable = 1;
end
if ~isfield(cfg_ana,'calcGGHF')
  cfg_ana.calcGGHF = 1;
end
if ~isfield(cfg_ana,'printTable_tex')
  cfg_ana.printTable_tex = 0;
end

% get the label info for this data struct
%if isfield(data.(exper.eventValues{1}).sub(1).ses(1).data,'label');
if isfield(ana.elec,'label');
  %lab = data.(exper.eventValues{1}).sub(1).ses(1).data.label;
  lab = ana.elec.label;
else
  error('label information not found in ana struct');
end

% set the channel information
if ~isfield(cfg_ana,'roi')
  error('Must specify either ROI names or channel names in cfg_ana.roi');
elseif isfield(cfg_ana,'roi')
  if ~iscell(cfg_ana.roi)
    cfg_ana.roi = {cfg_ana.roi};
  end
  cfg_ana.chan_str = sprintf(repmat(' %s',1,length(cfg_ana.roi)),cfg_ana.roi{:});
  
  if strcmp(cfg_ana.roi,'all')
    cfg_ana.roi = ft_channelselection(cfg_ana.roi,lab);
  end
end

% exclude the bad subjects from the subject count
cfg_ana.numSub = length(exper.subjects) - sum(exper.badSub);

% Initialize the matrices for storing the ANOVA data
%cfg_ana.anovamat = [];
cfg_ana.anovamat_iv1_iv2_iv3 = [];
cfg_ana.anovamat_iv1 = [];
cfg_ana.anovamat_iv3 = [];
cfg_ana.anovamat_iv2 = [];
cfg_ana.anovamat_iv1_iv3 = [];
cfg_ana.anovamat_iv1_iv2 = [];
cfg_ana.anovamat_iv2_iv3 = [];

% make sure cfg_ana.conditions is set correctly
if ~isfield(cfg_ana,'condMethod') || isempty(cfg_ana.condMethod)
  if ~iscell(cfg_ana.conditions) && (strcmp(cfg_ana.conditions,'all') || strcmp(cfg_ana.conditions,'all_across_types') || strcmp(cfg_ana.conditions,'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  elseif iscell(cfg_ana.conditions) && ~iscell(cfg_ana.conditions{1}) && length(cfg_ana.conditions) == 1 && (strcmp(cfg_ana.conditions{1},'all') || strcmp(cfg_ana.conditions{1},'all_across_types') || strcmp(cfg_ana.conditions{1},'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  elseif iscell(cfg_ana.conditions) && iscell(cfg_ana.conditions{1}) && length(cfg_ana.conditions{1}) == 1 && (strcmp(cfg_ana.conditions{1},'all') || strcmp(cfg_ana.conditions{1},'all_across_types') || strcmp(cfg_ana.conditions{1},'all_within_types'))
    cfg_ana.condMethod = 'pairwise';
  else
    cfg_ana.condMethod = [];
  end
end
cfg_ana.conditions = mm_ft_checkConditions(cfg_ana.conditions,ana,cfg_ana.condMethod);

% IV1
for iv1 = 1:length(cfg_ana.IV1.cond)
  % IV2
  for iv2 = 1:length(cfg_ana.IV2.cond)
%     cfg_ana.cond = cfg_ana.IV2.cond{iv2};
%     
%     if length(cfg_ana.cond) ~= length(cfg_ana.condCommon)
%       error('cfg_ana.condCommon and the number of events in this type of cfg_ana.conditions do not match');
%     end
    
    % IV3
    for iv3 = 1:length(cfg_ana.IV3.cond)
      
      evVal = sprintf('%s%s%s',cfg_ana.IV1.cond{iv1},cfg_ana.IV2.cond{iv2},cfg_ana.IV3.cond{iv3});
      if ~isfield(data,evVal)
        error('Condition name for %s was entered incorrectly!',evVal);
      end
      
      cfg_ana.values.(evVal) = nan(cfg_ana.numSub,length(exper.sessions));
      
      goodSubInd = 0;
      for sub = 1:length(exper.subjects)
        ses = 1;
        %for ses = 1:length(exper.sessions)
        if exper.badSub(sub,ses)
          %fprintf('Skipping bad subject: %s\n',exper.subjects{sub});
          continue
        else
          goodSubInd = goodSubInd + 1;
          
          % get the right channels (on an individual subject basis)
          if ismember(cfg_ana.roi,ana.elecGroupsStr)
            cfg_ana.channel = unique(cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi)}));
            
            % all good channels for this subject
            cfg_ana.chansel = ismember(data.(evVal).sub(sub).ses(ses).data.label,cfg_ana.channel);
            
            % common good channels
            cfg_ana.chansel = ismember(data.(evVal).sub(sub).ses(ses).data.label,cfg_ana.channel);
          else
            % find the channel indices for averaging
            cfg_ana.chansel = ismember(data.(evVal).sub(sub).ses(ses).data.label,cfg_ana.roi);
          end
          
          cfg_ana.timesel.(evVal) = find(data.(evVal).sub(sub).ses(ses).data.time >= cfg_ana.latency(1) & data.(evVal).sub(sub).ses(ses).data.time <= cfg_ana.latency(2));
          cfg_ana.values.(evVal)(goodSubInd,ses) = mean(mean(data.(evVal).sub(sub).ses(ses).data.(cfg_ana.parameter)(cfg_ana.chansel,cfg_ana.timesel.(evVal)),1),2);
        end
        
        % format for RMAOV3_mod: [data iv1 iv2 iv3 subNum]
        %cfg_ana.anovamat = cat(1,cfg_ana.anovamat,[cfg_ana.values.(cfg_ana.IV2.cond{typ}).(cfg_ana.cond{c})(goodSubInd,ses) r typ c goodSubInd]);
        
        cfg_ana.anovamat_iv1_iv2_iv3 = cat(1,cfg_ana.anovamat_iv1_iv2_iv3,[cfg_ana.values.(evVal)(goodSubInd,ses) iv1 iv2 iv3 goodSubInd]);
        cfg_ana.anovamat_iv1 = cat(1,cfg_ana.anovamat_iv1,[cfg_ana.values.(evVal)(goodSubInd,ses) iv1 goodSubInd]);
        cfg_ana.anovamat_iv2 = cat(1,cfg_ana.anovamat_iv2,[cfg_ana.values.(evVal)(goodSubInd,ses) iv2 goodSubInd]);
        cfg_ana.anovamat_iv3 = cat(1,cfg_ana.anovamat_iv3,[cfg_ana.values.(evVal)(goodSubInd,ses) iv3 goodSubInd]);
        cfg_ana.anovamat_iv1_iv3 = cat(1,cfg_ana.anovamat_iv1_iv3,[cfg_ana.values.(evVal)(goodSubInd,ses) iv1 iv3 goodSubInd]);
        cfg_ana.anovamat_iv1_iv2 = cat(1,cfg_ana.anovamat_iv1_iv2,[cfg_ana.values.(evVal)(goodSubInd,ses) iv1 iv2 goodSubInd]);
        cfg_ana.anovamat_iv2_iv3 = cat(1,cfg_ana.anovamat_iv2_iv3,[cfg_ana.values.(evVal)(goodSubInd,ses) iv2 iv3 goodSubInd]);
        %end % ses
      end % sub
    end % iv3
  end % iv2
end % iv1

% need to do Greenhouse-Geisser or Huynh-Feldt correction?

if length(cfg_ana.IV1.cond) > 2 ||  length(cfg_ana.IV2.cond) > 2 ||  length(cfg_ana.IV3.cond) > 2
  cfg_ana.calcGGHF = 1;
else
  cfg_ana.calcGGHF = 0;
end

% if length(cfg_ana.conditions) > 1
%   if length(cfg_ana.IV1.cond) > 2 || length(cfg_ana.IV2.cond) > 2 || sign(sum(cellfun('length',cfg_ana.conditions) > 2))
%     cfg_ana.calcGGHF = 1;
%   else
%     cfg_ana.calcGGHF = 0;
%   end
% else
%   if length(cfg_ana.IV1.cond) > 2 || length(cfg_ana.IV2.cond) > 2 || length(cfg_ana.conditions) > 2
%     cfg_ana.calcGGHF = 1;
%   else
%     cfg_ana.calcGGHF = 0;
%   end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('====================== RMAOV33_mod =========================\n');

%if ~strcmp(cfg_ana.IV1.type,'latency') &&  ~strcmp(cfg_ana.IV1.type,'latency') &&  ~strcmp(cfg_ana.IV1.type,'latency')
%end

fprintf('%.1f--%.1f s,%s: IV1: %s (%d;%s), IV2: %s (%d;%s), IV3: %s (%d;%s)\n',...
  cfg_ana.latency(1),cfg_ana.latency(2),cfg_ana.chan_str,...
  cfg_ana.IV1.name,length(cfg_ana.IV1.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV1.cond)),cfg_ana.IV1.cond{:}),...
  cfg_ana.IV2.name,length(cfg_ana.IV2.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV2.cond)),cfg_ana.IV2.cond{:}),...
  cfg_ana.IV3.name,length(cfg_ana.IV3.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV3.cond)),cfg_ana.IV3.cond{:}));
[cfg_ana.p] = RMAOV33_mod(cfg_ana.anovamat_iv1_iv2_iv3,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,cfg_ana.printTable_tex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do a bunch of tests collapsing over one IV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IV1xIV2
fprintf('====================== RMAOV2_mod =========================\n');
fprintf('%.1f--%.1f s,%s: IV1: %s (%d;%s), IV2: %s (%d;%s), collapsing over %s (%d;%s)\n',...
  cfg_ana.latency(1),cfg_ana.latency(2),cfg_ana.chan_str,...
  cfg_ana.IV1.name,length(cfg_ana.IV1.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV1.cond)),cfg_ana.IV1.cond{:}),...
  cfg_ana.IV2.name,length(cfg_ana.IV2.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV2.cond)),cfg_ana.IV2.cond{:}),...
  cfg_ana.IV3.name,length(cfg_ana.IV3.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV3.cond)),cfg_ana.IV3.cond{:}));
[cfg_ana.p,cfg_ana.X] = RMAOV2_mod(cfg_ana.anovamat_iv1_iv2,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,0,1);
% at each level of IV1, do pairwise comparisons of IV2's levels
levels1 = unique(cfg_ana.X(:,2));
tt_combos1 = nchoosek(unique(cfg_ana.X(:,3)),2);
for i = 1:length(levels1)
  for j = 1:size(tt_combos1,1)
    cond1 = cfg_ana.X(cfg_ana.X(:,2) == levels1(i) & cfg_ana.X(:,3) == tt_combos1(j,1),1);
    cond2 = cfg_ana.X(cfg_ana.X(:,2) == levels1(i) & cfg_ana.X(:,3) == tt_combos1(j,2),1);
    [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
    cohens_d = mm_effect_size('within',cond1,cond2);
    fprintf('ttest: %s: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
      cfg_ana.IV1.cond{levels1(i,1)},...
      cfg_ana.IV2.cond{tt_combos1(j,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
      cfg_ana.IV2.cond{tt_combos1(j,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
      cfg_ana.tstats.df,...
      cfg_ana.tstats.tstat,...
      cohens_d,...
      cfg_ana.tstats.sd,...
      (cfg_ana.tstats.sd/sqrt(length(cond1))),...
      cfg_ana.tp);
    if cfg_ana.tp < cfg_ana.alpha
      fprintf(' *');
    end
    fprintf('\n');
  end
end
fprintf('\n');
% at each level of IV2, do pairwise comparisons of IV1's levels
levels1 = unique(cfg_ana.X(:,3));
tt_combos1 = nchoosek(unique(cfg_ana.X(:,2)),2);
for i = 1:length(levels1)
  for j = 1:size(tt_combos1,1)
    cond1 = cfg_ana.X(cfg_ana.X(:,3) == levels1(i) & cfg_ana.X(:,2) == tt_combos1(j,1),1);
    cond2 = cfg_ana.X(cfg_ana.X(:,3) == levels1(i) & cfg_ana.X(:,2) == tt_combos1(j,2),1);
    [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
    cohens_d = mm_effect_size('within',cond1,cond2);
    fprintf('ttest: %s: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
      cfg_ana.IV2.cond{levels1(i,1)},...
      cfg_ana.IV1.cond{tt_combos1(j,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
      cfg_ana.IV1.cond{tt_combos1(j,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
      cfg_ana.tstats.df,...
      cfg_ana.tstats.tstat,...
      cohens_d,...
      cfg_ana.tstats.sd,...
      (cfg_ana.tstats.sd/sqrt(length(cond1))),...
      cfg_ana.tp);
    if cfg_ana.tp < cfg_ana.alpha
      fprintf(' *');
    end
    fprintf('\n');
  end
end
fprintf('\n');

% IV1xIV3
fprintf('====================== RMAOV2_mod =========================\n');
fprintf('%.1f--%.1f s,%s: IV1: %s (%d;%s), IV2: %s (%d;%s), collapsing over %s (%d;%s)\n',...
  cfg_ana.latency(1),cfg_ana.latency(2),cfg_ana.chan_str,...
  cfg_ana.IV1.name,length(cfg_ana.IV1.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV1.cond)),cfg_ana.IV1.cond{:}),...
  cfg_ana.IV3.name,length(cfg_ana.IV3.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV3.cond)),cfg_ana.IV3.cond{:}),...
  cfg_ana.IV2.name,length(cfg_ana.IV2.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV2.cond)),cfg_ana.IV2.cond{:}));
[cfg_ana.p,cfg_ana.X] = RMAOV2_mod(cfg_ana.anovamat_iv1_iv3,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,0,1);
% at each level of IV1, do pairwise comparisons of IV3's levels
levels1 = unique(cfg_ana.X(:,2));
tt_combos1 = nchoosek(unique(cfg_ana.X(:,3)),2);
for i = 1:length(levels1)
  for j = 1:size(tt_combos1,1)
    cond1 = cfg_ana.X(cfg_ana.X(:,2) == levels1(i) & cfg_ana.X(:,3) == tt_combos1(j,1),1);
    cond2 = cfg_ana.X(cfg_ana.X(:,2) == levels1(i) & cfg_ana.X(:,3) == tt_combos1(j,2),1);
    [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
    cohens_d = mm_effect_size('within',cond1,cond2);
    fprintf('ttest: %s: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
      cfg_ana.IV1.cond{levels1(i,1)},...
      cfg_ana.IV3.cond{tt_combos1(j,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
      cfg_ana.IV3.cond{tt_combos1(j,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
      cfg_ana.tstats.df,...
      cfg_ana.tstats.tstat,...
      cohens_d,...
      cfg_ana.tstats.sd,...
      (cfg_ana.tstats.sd/sqrt(length(cond1))),...
      cfg_ana.tp);
    if cfg_ana.tp < cfg_ana.alpha
      fprintf(' *');
    end
    fprintf('\n');
  end
end
fprintf('\n');
% at each level of IV3, do pairwise comparisons of IV1's levels
levels1 = unique(cfg_ana.X(:,3));
tt_combos1 = nchoosek(unique(cfg_ana.X(:,2)),2);
for i = 1:length(levels1)
  for j = 1:size(tt_combos1,1)
    cond1 = cfg_ana.X(cfg_ana.X(:,3) == levels1(i) & cfg_ana.X(:,2) == tt_combos1(j,1),1);
    cond2 = cfg_ana.X(cfg_ana.X(:,3) == levels1(i) & cfg_ana.X(:,2) == tt_combos1(j,2),1);
    [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
    cohens_d = mm_effect_size('within',cond1,cond2);
    fprintf('ttest: %s: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
      cfg_ana.IV3.cond{levels1(i,1)},...
      cfg_ana.IV1.cond{tt_combos1(j,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
      cfg_ana.IV1.cond{tt_combos1(j,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
      cfg_ana.tstats.df,...
      cfg_ana.tstats.tstat,...
      cohens_d,...
      cfg_ana.tstats.sd,...
      (cfg_ana.tstats.sd/sqrt(length(cond1))),...
      cfg_ana.tp);
    if cfg_ana.tp < cfg_ana.alpha
      fprintf(' *');
    end
    fprintf('\n');
  end
end
fprintf('\n');

% IV2xIV3
fprintf('====================== RMAOV2_mod =========================\n');
fprintf('%.1f--%.1f s,%s: IV1: %s (%d;%s), IV2: %s (%d;%s), collapsing over %s (%d;%s)\n',...
  cfg_ana.latency(1),cfg_ana.latency(2),cfg_ana.chan_str,...
  cfg_ana.IV2.name,length(cfg_ana.IV2.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV2.cond)),cfg_ana.IV2.cond{:}),...
  cfg_ana.IV3.name,length(cfg_ana.IV3.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV3.cond)),cfg_ana.IV3.cond{:}),...
  cfg_ana.IV1.name,length(cfg_ana.IV1.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV1.cond)),cfg_ana.IV1.cond{:}));
[cfg_ana.p,cfg_ana.X] = RMAOV2_mod(cfg_ana.anovamat_iv2_iv3,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,0,1);
fprintf('\n');
% at each level of IV2, do pairwise comparisons of IV3's levels
levels1 = unique(cfg_ana.X(:,2));
tt_combos1 = nchoosek(unique(cfg_ana.X(:,3)),2);
for i = 1:length(levels1)
  for j = 1:size(tt_combos1,1)
    cond1 = cfg_ana.X(cfg_ana.X(:,2) == levels1(i) & cfg_ana.X(:,3) == tt_combos1(j,1),1);
    cond2 = cfg_ana.X(cfg_ana.X(:,2) == levels1(i) & cfg_ana.X(:,3) == tt_combos1(j,2),1);
    [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
    cohens_d = mm_effect_size('within',cond1,cond2);
    fprintf('ttest: %s: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
      cfg_ana.IV2.cond{levels1(i,1)},...
      cfg_ana.IV3.cond{tt_combos1(j,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
      cfg_ana.IV3.cond{tt_combos1(j,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
      cfg_ana.tstats.df,...
      cfg_ana.tstats.tstat,...
      cohens_d,...
      cfg_ana.tstats.sd,...
      (cfg_ana.tstats.sd/sqrt(length(cond1))),...
      cfg_ana.tp);
    if cfg_ana.tp < cfg_ana.alpha
      fprintf(' *');
    end
    fprintf('\n');
  end
end
fprintf('\n');
% at each level of IV3, do pairwise comparisons of IV2's levels
levels1 = unique(cfg_ana.X(:,3));
tt_combos1 = nchoosek(unique(cfg_ana.X(:,2)),2);
for i = 1:length(levels1)
  for j = 1:size(tt_combos1,1)
    cond1 = cfg_ana.X(cfg_ana.X(:,3) == levels1(i) & cfg_ana.X(:,2) == tt_combos1(j,1),1);
    cond2 = cfg_ana.X(cfg_ana.X(:,3) == levels1(i) & cfg_ana.X(:,2) == tt_combos1(j,2),1);
    [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
    cohens_d = mm_effect_size('within',cond1,cond2);
    fprintf('ttest: %s: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
      cfg_ana.IV3.cond{levels1(i,1)},...
      cfg_ana.IV2.cond{tt_combos1(j,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
      cfg_ana.IV2.cond{tt_combos1(j,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
      cfg_ana.tstats.df,...
      cfg_ana.tstats.tstat,...
      cohens_d,...
      cfg_ana.tstats.sd,...
      (cfg_ana.tstats.sd/sqrt(length(cond1))),...
      cfg_ana.tp);
    if cfg_ana.tp < cfg_ana.alpha
      fprintf(' *');
    end
    fprintf('\n');
  end
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do a bunch of tests collapsing over two IVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IV1
fprintf('====================== RMAOV1_mod =========================\n');
fprintf('%.1f--%.1f s,%s: IV1: %s (%d;%s), collapsing over %s (%d;%s) and %s (%d;%s)\n',...
  cfg_ana.latency(1),cfg_ana.latency(2),cfg_ana.chan_str,...
  cfg_ana.IV1.name,length(cfg_ana.IV1.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV1.cond)),cfg_ana.IV1.cond{:}),...
  cfg_ana.IV2.name,length(cfg_ana.IV2.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV2.cond)),cfg_ana.IV2.cond{:}),...
  cfg_ana.IV3.name,length(cfg_ana.IV3.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV3.cond)),cfg_ana.IV3.cond{:}));
[cfg_ana.p,cfg_ana.X] = RMAOV1_mod(cfg_ana.anovamat_iv1,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,0,1);
fprintf('\n');
tt_combos = nchoosek(unique(cfg_ana.X(:,2)),2);
for i = 1:size(tt_combos,1)
  cond1 = cfg_ana.X(cfg_ana.X(:,2) == tt_combos(i,1),1);
  cond2 = cfg_ana.X(cfg_ana.X(:,2) == tt_combos(i,2),1);
  [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
  cohens_d = mm_effect_size('within',cond1,cond2);
  fprintf('ttest: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
    cfg_ana.IV1.cond{tt_combos(i,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
    cfg_ana.IV1.cond{tt_combos(i,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
    cfg_ana.tstats.df,...
    cfg_ana.tstats.tstat,...
    cohens_d,...
    cfg_ana.tstats.sd,...
    (cfg_ana.tstats.sd/sqrt(length(cond1))),...
    cfg_ana.tp);
  if cfg_ana.tp < cfg_ana.alpha
    fprintf(' *');
  end
  fprintf('\n');
end
fprintf('\n');

% IV2
fprintf('====================== RMAOV1_mod =========================\n');
fprintf('%.1f--%.1f s,%s: IV1: %s (%d;%s), collapsing over %s (%d;%s) and %s (%d;%s)\n',...
  cfg_ana.latency(1),cfg_ana.latency(2),cfg_ana.chan_str,...
  cfg_ana.IV2.name,length(cfg_ana.IV2.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV2.cond)),cfg_ana.IV2.cond{:}),...
  cfg_ana.IV1.name,length(cfg_ana.IV1.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV1.cond)),cfg_ana.IV1.cond{:}),...
  cfg_ana.IV3.name,length(cfg_ana.IV3.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV3.cond)),cfg_ana.IV3.cond{:}));
[cfg_ana.p,cfg_ana.X] = RMAOV1_mod(cfg_ana.anovamat_iv2,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,0,1);
fprintf('\n');
tt_combos = nchoosek(unique(cfg_ana.X(:,2)),2);
for i = 1:size(tt_combos,1)
  cond1 = cfg_ana.X(cfg_ana.X(:,2) == tt_combos(i,1),1);
  cond2 = cfg_ana.X(cfg_ana.X(:,2) == tt_combos(i,2),1);
  [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
  cohens_d = mm_effect_size('within',cond1,cond2);
  fprintf('ttest: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
    cfg_ana.IV2.cond{tt_combos(i,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
    cfg_ana.IV2.cond{tt_combos(i,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
    cfg_ana.tstats.df,...
    cfg_ana.tstats.tstat,...
    cohens_d,...
    cfg_ana.tstats.sd,...
    (cfg_ana.tstats.sd/sqrt(length(cond1))),...
    cfg_ana.tp);
  if cfg_ana.tp < cfg_ana.alpha
    fprintf(' *');
  end
  fprintf('\n');
end
fprintf('\n');

% IV3
fprintf('====================== RMAOV1_mod =========================\n');
fprintf('%.1f--%.1f s,%s: IV1: %s (%d;%s), collapsing over %s (%d;%s) and %s (%d;%s)\n',...
  cfg_ana.latency(1),cfg_ana.latency(2),cfg_ana.chan_str,...
  cfg_ana.IV3.name,length(cfg_ana.IV3.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV3.cond)),cfg_ana.IV3.cond{:}),...
  cfg_ana.IV1.name,length(cfg_ana.IV1.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV1.cond)),cfg_ana.IV1.cond{:}),...
  cfg_ana.IV2.name,length(cfg_ana.IV2.cond),sprintf(repmat(' %s',1,length(cfg_ana.IV2.cond)),cfg_ana.IV2.cond{:}));
[cfg_ana.p,cfg_ana.X] = RMAOV1_mod(cfg_ana.anovamat_iv3,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,0,1);
fprintf('\n');
tt_combos = nchoosek(unique(cfg_ana.X(:,2)),2);
for i = 1:size(tt_combos,1)
  cond1 = cfg_ana.X(cfg_ana.X(:,2) == tt_combos(i,1),1);
  cond2 = cfg_ana.X(cfg_ana.X(:,2) == tt_combos(i,2),1);
  [cfg_ana.th,cfg_ana.tp,cfg_ana.tci,cfg_ana.tstats] = ttest(cond1,cond2,cfg_ana.alpha,'both');
  cohens_d = mm_effect_size('within',cond1,cond2);
  fprintf('ttest: %s (M=%.3f, SEM=%.3f) vs\t%s (M=%.3f, SEM=%.3f):\tt(%d)=%.4f, d=%.3f, SD=%.2f, SEM=%.2f, p=%.10f',...
    cfg_ana.IV3.cond{tt_combos(i,1)},mean(cond1,1),std(cond1)/sqrt(length(cond1)),...
    cfg_ana.IV3.cond{tt_combos(i,2)},mean(cond2,1),std(cond2)/sqrt(length(cond2)),...
    cfg_ana.tstats.df,...
    cfg_ana.tstats.tstat,...
    cohens_d,...
    cfg_ana.tstats.sd,...
    (cfg_ana.tstats.sd/sqrt(length(cond1))),...
    cfg_ana.tp);
  if cfg_ana.tp < cfg_ana.alpha
    fprintf(' *');
  end
  fprintf('\n');
end
fprintf('\n');

end
