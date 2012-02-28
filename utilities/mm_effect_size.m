function [d] = mm_effect_size(method,varargin)
%MM_EFFECT_SIZE - Calculate effect size (Cohen's d)
%
%  [d] = mm_effect_size(method,varargin)
%
% A number of calculation methods are possible, and you need to carefully
% choose the method depending on the data. I HOPE these are all correct,
% but I do see a few small discrepancies (e.g., 1a vs 1c/1d; 2a vs 2b/2c).
%
% NB: If you're on a Mac, I recommend using G*Power 3 to double-check your
%     answers:
%     http://www.psycho.uni-duesseldorf.de/abteilungen/aap/gpower3/
%
% NB: There was a thread on the FieldTrip list about effect sizes:
%     http://mailman.science.ru.nl/pipermail/fieldtrip/2011-September/004224.html
%     They cited this paper, which seems useful: Rosnow & Rosenthal, Effect
%     Sizes for Experimenting Psychologists, Canadian Journal of
%     Experimental Psychology, 2003, 57:3, 221-237.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   method   = String describing the method. See below.
%   varargin = Depends on the method. See below.
%
% Output:
%   d        = Cohen's d
%
% Cohen's d interpretation (Cohen, 1988):
%   0.2 to 0.3 = "small" effect
%   around 0.5 = "medium" effect
%   0.8 to infinity = a "large" effect
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1: BETWEEN-SUBJECTS DATA (independent groups/samples)
%
% NB: There is a small difference between 1a and 1c/1d. I don't understand
%     why this happens. 1b agrees with G*Power 3 and is a slightly
%     different method than 1a.
%
% 1a. Divide the difference of the means by the pooled standard
%     deviation (independent-samples, equal or unequal sample sizes). Uses
%     N (population) as denominator for pooled SD (Cohen's d does this).
%     This does not agree with G*Power 3 (see 1b).
%
%   method = 'between_d'
%   data1  = data for group 1
%   data2  = data for group 2
%
% 1b. Divide the difference of the means by the pooled standard
%     deviation (independent-samples, equal or unequal sample sizes). Uses
%     N-1 (within-sample, unbiased) as denominator for pooled SD (Hedge's g
%     does this; Hedges and Olkin, 1985). This agrees with G*Power 3.
%
%   method = 'between_g'
%   data1  = data for group 1
%   data2  = data for group 2
%
% 1c. Using t-test results (we do not run TTEST2.M for you)
%
%   method = 'between_t'
%   t      = t-value
%   df     = degrees of freedom
%
% 1d. Using t-test results; we run TTEST2.M for you. Assumes equal
%     variances. This is the same as 1c, unless if you want unequal
%     variance, in which case you should use method 1c with the 'unequal'
%     option.
%
%   method = 'between_t_run'
%   data1  = data for group 1
%   data2  = data for group 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2: WITHIN-SUBJECTS DATA (dependent groups / paired samples)
%
% NB: 2b and 2c do not correct for any correlation between the two samples,
%     but this is supposedly OK because a paired-sample t-test takes the
%     correlated design into account. However, this page
%     <http://www.uccs.edu/~faculty/lbecker/es.htm> says that it is not OK,
%     so you probably shouldn't use 2b or 2c. I'm not sure, but it sounds
%     like the page also says to use the between-subjects method (e.g.,
%     1a/1b).
%
% 2a. Divide the difference of the means by the pooled standard deviation,
%     corrected for correlation between the two samples (using CORR.M).
%     This agrees with G*Power 3.
%
%   method = 'within'
%   data1  = data for group 1
%   data2  = data for group 2
%
% 2b. Using t-test results (we do not run TTEST.M for you).  This does not
%     agree with G*Power 3, but it was mentioned on the FieldTrip email
%     list linked to at the top of this help.
%
%   method = 'within_t'
%   t      = t-value
%   df     = degrees of freedom
%
% 2c. Using t-test results; we run TTEST.M for you. This is the same as 2b,
%
%   method = 'within_t_run'
%   data1  = data for group 1
%   data2  = data for group 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3: SINGLE-SAMPLE ONLY, comparison to a constant
%
% For control condition comparison only (e.g., vs chance; only 1 standard
% deviation is available)
%
%   method   = 'single'
%   data     = data for group
%   constant = single constant number to which to compare data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also: TTEST, TTEST2, CORR, STD
%

% TODO: report correlation coefficient



% The root mean squared standard deviation method is commented out because
% it is the same as using the pooled standard deviation, but doesn't handle
% unequal sample sizes

% if strcmp(method,'between_equal')
%   % 1a. Using the root mean squared standard deviation
%   % (independent-samples, equal sample sizes)
%   
%   data1 = varargin{1};
%   data2 = varargin{2};
%   
%   d = abs(mean(data1) - mean(data2)) / ...
%     sqrt((std(data1)^2 + std(data2)^2) / 2);


if strcmp(method,'between_d')
  % 1a. Using the pooled standard deviation (independent-samples, equal or
  % unequal sample sizes)
  
  data1 = varargin{1};
  data2 = varargin{2};
  
  data1Len = length(data1);
  data2Len = length(data2);
  
  s_pooled = sqrt((((data1Len - 1) * std(data1)^2) + ((data2Len - 1) * std(data2)^2)) / (data1Len + data2Len));
  
  d = abs(mean(data1) - mean(data2)) / s_pooled;
  
elseif strcmp(method,'between_g')
  % 1b. Using the pooled standard deviation (independent-samples, equal or
  % unequal sample sizes)
  
  data1 = varargin{1};
  data2 = varargin{2};
  
  data1Len = length(data1) - 1;
  data2Len = length(data2) - 1;
  
  s_pooled = sqrt(((data1Len * std(data1)^2) + (data2Len * std(data2)^2)) / (data1Len + data2Len));
  
  d = abs(mean(data1) - mean(data2)) / s_pooled;
  
elseif strcmp(method,'between_t')
  % 1c. Using t-test results and pooled standard deviation
  % (independent-samples, unequal sample sizes)
  
  t = varargin{1};
  df = varargin{2};
  
  d = 2 * abs(t) / sqrt(df);
  
elseif strcmp(method,'between_t_run')
  % 1d. run the ttest2
  
  data1 = varargin{1};
  data2 = varargin{2};
  
  [h,p,ci,stats] = ttest2(data1,data2,0.05,'both');
  d = 2 * abs(stats.tstat) / sqrt(stats.df);
  
elseif strcmp(method,'within')
  % 2: WITHIN-SUBJECTS DATA (paired-sample t-test)
  
  % 2a. Divide the difference of the means by the pooled standard
  % deviation, corrected for correlation
  
  data1 = varargin{1};
  data2 = varargin{2};
  
  if size(data1,1) == 1
    [rho] = corr(data1',data2');
  elseif size(data1,2) == 1
    [rho] = corr(data1,data2);
  else
    error('Not sure how to use correlation with these data.');
  end
  
  s_pooled = sqrt((std(data1)^2 + std(data2)^2) - (2 * rho * std(data1) * std(data2)));
  
  d = abs(mean(data1) - mean(data2)) / s_pooled;
  
elseif strcmp(method,'within_t')
  % 2b. Using t-test results and pooled standard deviation
  %     (dependent-samples)
  
  t = varargin{1};
  df = varargin{2};
  
  d = abs(t) / sqrt(df);
  
elseif strcmp(method,'within_t_run')
  % 2c. run the ttest
  
  data1 = varargin{1};
  data2 = varargin{2};
  
  [h,p,ci,stats] = ttest(data1,data2,0.05,'both');
  d = abs(stats.tstat) / sqrt(stats.df);
  
elseif strcmp(method,'single')
  % 3: SINGLE-SAMPLE ONLY, comparison to a constant
  
  data = varargin{1};
  constant = varargin{2};
  
  d = abs(mean(data) - constant) / std(data);
  
end

end

