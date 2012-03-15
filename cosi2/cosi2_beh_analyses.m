
%% initialize randomness

s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);

%% data

exper.cosi2.subjects = {...
  'COSI2008';
  'COSI2009';
  'COSI2010';
  'COSI2012';
  'COSI2013';
  'COSI2015';
  'COSI2016';
  'COSI2017';
  'COSI2018';
  'COSI2019';
  'COSI2020';
  'COSI2021';
  'COSI2022';
  'COSI2023';
  'COSI2024';
  'COSI2025';
  'COSI2026';
  'COSI2027';
  'COSI2028';
  'COSI2029';
  'COSI2030';
  'COSI2032';
  'COSI2033';
  'COSI2034';
  'COSI2035';
  'COSI2036';
  'COSI2037';
  'COSI2038';
  'COSI2039';
  'COSI2040';
  'COSI2042';
  'COSI2043';
  'COSI2044';
  'COSI2045';
  };


exper.cosi2.color.goodSub = {
    'COSI2012';
    'COSI2013';
    'COSI2017';
    'COSI2018';
    'COSI2024';
    'COSI2028';
    'COSI2032';
    'COSI2034';
    'COSI2036';
    'COSI2037';
    'COSI2042';
    'COSI2044';
    'COSI2045';
  };

exper.cosi2.side.goodSub = {
    'COSI2010';
    'COSI2015';
    'COSI2019';
    'COSI2021';
    'COSI2022';
    'COSI2023';
    'COSI2026';
    'COSI2027';
    'COSI2030';
    'COSI2033';
    'COSI2035';
    'COSI2040';
    'COSI2043';
  };

% % exclude for a number of reasons
% exper.cosi2.badSub = {'COSI2008','COSI2009','COSI2020','COSI2025', 'COSI2011', 'COSI2014', 'COSI2031', 'COSI2041', 'COSI2016', 'COSI2029', 'COSI2039', 'COSI2038'};

% % 8, 9, 20, 25: no F responses in one color/side SC/SI bin

% % 11, 14, 31, 41: no session_1

% % 16, 29 have fewer than 15 trials for Side-SI

% % 39 has fewer than 15 trials for Color-SI

% % 38: potentially bad session_1 (puker)

% exclude based on automatic equating of accuracy
exper.cosi2.color.badSub = exper.cosi2.subjects(~ismember(exper.cosi2.subjects,exper.cosi2.color.goodSub));
exper.cosi2.side.badSub = exper.cosi2.subjects(~ismember(exper.cosi2.subjects,exper.cosi2.side.goodSub));


COSI2_dp_color = [0.7234, 5.5809, 0.596, 1.0686, 1.386, 0.8479, 2.6736, 1.0338, 1.0819, 0.1147, 1.9433, 1.1353, 0.1524, 0.2168, 0.6119, 2.1662, 0.238, 0.8939, 0.604, 1.35, 0.549, 1.7306, 0.6466, 0.8577, 0.7876, 0.994, 0.3608, 1.044, 3.0561, 0.4539, 0.7793, 0.1625, 0.4866, 0.5906];
COSI2_dp_side = [1.7697, 3.2063, 1.4826, 2.5392, 2.5688, 1.4101, 3.6724, 1.5205, 0.5383, 0.5133, 3.1705, 1.2578, 0.6185, 1.2063, 1.5154, 2.7215, 1.1305, 1.8154, 1.6792, 2.8938, 1.5142, 2.6603, 1.4934, 1.2753, 1.0709, 2.5971, 0.5229, 1.2748, 3.1693, 0.6425, 2.3447, 0.4569, 0.9257, 1.4539];

COSI2_F_color = [0.4667, 0.9988, 0.5545, 0.5143, 0.5333, 0.5, 0.6441, 0.5102, 0.5781, 0.449, 0.0013, 0.4231, 0.5333, 0.4375, 0.5263, 0.9988, 0.4198, 0.5278, 0.5769, 0.5352, 0.52, 0.7143, 0.4643, 0.5882, 0.4844, 0.5625, 0.5427, 0.6047, 0.75, 0.5556, 0.5392, 0.5079, 0.5, 0.5088];
COSI2_F_side = [0.9988, 0.8235, 0.5333, 0.6765, 0.5294, 0.4921, 0.7895, 0.5914, 0.5, 0.5051, 0.5, 0.5135, 0.4528, 0.5714, 0.6154, 0.9988, 0.5319, 0.5319, 0.5983, 0.62, 0.5063, 0.8, 0.5606, 0.5584, 0.4833, 0.7143, 0.5354, 0.6111, 0.625, 0.5063, 0.5873, 0.4797, 0.5405, 0.5714];


nTries = 10000;
p_thresh = 0.1;

%% try to equate

%foundIt = false;

for i = 1:nTries
  fprintf('.');
  subNums = randperm(length(exper.cosi2.subjects));
  subNums = subNums(~ismember(exper.cosi2.subjects(subNums),exper.cosi2.badSub));
  
  colorHalf = subNums(1:length(subNums)/2);
  sideHalf = subNums((length(subNums)/2) + 1:end);
  
  [h,p,ci,stats] = ttest2(COSI2_dp_color(colorHalf),COSI2_dp_side(sideHalf),0.05,'both');
  if p > p_thresh
    [h_cs,p_cs,ci_cs,stats_cs] = ttest2(COSI2_F_color(colorHalf),COSI2_F_side(sideHalf),0.05,'both');
    if p_cs < .05
      fprintf('Found p > %.2f and p_cs < .05 after %d tries.\n',p_thresh,i);
      
      fprintf('ttest2: Color half source d''=%.2f vs Side half source d''=%.2f: t(%d)=%.2f, p=%.3f',mean(COSI2_dp_color(colorHalf)),mean(COSI2_dp_side(sideHalf)),stats.df,stats.tstat,p);
      
      fprintf('\nColor half:');
      sort(exper.cosi2.subjects(colorHalf))
      fprintf('\nSide half:');
      sort(exper.cosi2.subjects(sideHalf))
      
      [h,p,ci,stats] = ttest2(COSI2_F_color(colorHalf),COSI2_F_side(sideHalf),0.05,'both');
      fprintf('ttest2: Color half F=%.2f vs Side half F=%.2f: t(%d)=%.2f, p=%.3f\n',mean(COSI2_F_color(colorHalf)),mean(COSI2_F_side(sideHalf)),stats.df,stats.tstat,p);
      
      [h,p,ci,stats] = ttest(COSI2_F_color(colorHalf),repmat(0.5,size(colorHalf)),0.05,'both');
      fprintf('ttest: Color half F=%.2f vs chance (%.1f): t(%d)=%.2f, p=%.3f\n',mean(COSI2_F_color(colorHalf)),0.5,stats.df,stats.tstat,p);
      
      [h,p,ci,stats] = ttest(COSI2_F_side(sideHalf),repmat(0.5,size(sideHalf)),0.05,'both');
      fprintf('ttest: Side half F=%.2f vs chance (%.1f): t(%d)=%.2f, p=%.3f\n',mean(COSI2_F_side(sideHalf)),0.5,stats.df,stats.tstat,p);
      
      %foundIt = true;
      break
    end
  end
  if i == nTries
    fprintf('Hit the %d-try limit.\n',i);
  end
end

%% ERP voltage

% 300--500 ms

volt.cosi2.RCR.FS = [];
volt.cosi2.RHSC.FS = [];
volt.cosi2.RHSI.FS = [];

volt.cosi2.RCR.LAS = [];
volt.cosi2.RHSC.LAS = [];
volt.cosi2.RHSI.LAS = [];

volt.cosi2.RCR.RAS = [];
volt.cosi2.RHSC.RAS = [];
volt.cosi2.RHSI.RAS = [];

% 500--800 ms

volt.cosi2.RCR.LPS = [];
volt.cosi2.RHSC.LPS = [];
volt.cosi2.RHSI.LPS = [];

volt.cosi2.RCR.RPS = [];
volt.cosi2.RHSC.RPS = [];
volt.cosi2.RHSI.RPS = [];



%% Set up the rest of the ANOVA

% IMPORTANT: This requires an equal number of subjects in each experiment

% IV1 (between-subjects factor): Experiment
exper.name = {'sosi','soco'};
% IV2 (within-subjects RM factor): ROI
exper.roi = {'LAS','RAS'};
% IV3 (within-subjects RM factor): Accuracy Condition
exper.cond = {'RCR','RHSC','RHSI'};
%exper.cond = {'RHSC','RHSI'};

% anovamat - data matrix (Size of matrix must be n-by-5;dependent variable=column 1;
%            independent variable 1=column 2;independent variable 2 (within-subjects)=column 3;
%            independent variable 3 (within-subjects)=column 4; subject=column 5 [Be careful
%            how the code for the subjects are entered.]) 

anovamat = [];

for e = 1:length(exper.name)
  goodSubInd = 0;
  for s = 1:length(exper.(exper.name{e}).subjects)
    if ~ismember(exper.(exper.name{e}).subjects(s),exper.(exper.name{e}).badSub)
      goodSubInd = goodSubInd + 1;
      for r = 1:length(exper.roi)
        for c = 1:length(exper.cond)
          anovamat = [anovamat; volt.(exper.name{e}).(exper.cond{c}).(exper.roi{r})(s) e r c goodSubInd];
        end
      end
    end
  end
end

alpha = 0.05;
%showtable = 1;
%calcGGHF = 0;

%[P] = RMAOV32_mod(anovamat,alpha,showtable,calcGGHF);
[P] = RMAOV32(anovamat,alpha);


%% dependent-samples t-test (Within an experiment)

alpha = 0.05;

rois = {'LAS','RAS'};
expNames = {'soco','sosi'};
condPairs = {{'RHSC','RHSI'},{'RHSC','RCR'},{'RHSI','RCR'}};

for e = 1:length(expNames)
  expName = expNames{e};
  for r = 1:length(rois)
    roi = rois{r};
    for c = 1:length(condPairs)
      conds = condPairs{c};
      [h,p,ci,stats] = ttest(...
        volt.(expName).(conds{1}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub)),...
        volt.(expName).(conds{2}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub)),...
        alpha,'both');
      fprintf('ttest: %s: %s: %s (M=%.2f) vs %s (M=%.2f): t(%d)=%.4f, p=%.10f',roi,expName,...
        conds{1},mean(volt.(expName).(conds{1}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub))),...
        conds{2},mean(volt.(expName).(conds{2}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub))),...
        stats.df,stats.tstat,p);
      if p < alpha
        fprintf(' *');
      end
      fprintf('\n');
    end
  end
end

