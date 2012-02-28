function [d] = mm_effect_size(method,varargin)
%MM_EFFECT_SIZE - Calculate effect size (Cohen's d)
%
%  [d] = mm_effect_size(method,varargin)
%
% A number of calculation methods are possible, and you need to carefully
% choose the method depending on the data. I HOPE these are all correct,
% but I do see a small discrepancy betwen 1a and 1c/1d.
%
% NB: if you're on a Mac I recommend using G*Power 3 to double-check your
%     answers:
%     http://www.psycho.uni-duesseldorf.de/abteilungen/aap/gpower3/
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
% Cohen's d interpretation (Cohen, 1969):
%   0.2 to 0.3 = "small" effect
%   around 0.5 = "medium" effect
%   0.8 to infinity = a "large" effect
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1: BETWEEN-SUBJECTS DATA (independent-sample t-test)
%
% NB: There is a small difference between 1a and 1b/1c. I don't understand
%     why this happens. 1b agrees with G*Power 3 and is a slightly
%     different method than 1a
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
%     does this). This agrees with G*Power 3.
%
%   method = 'between_g'
%   data1  = data for group 1
%   data2  = data for group 2
%
% 1c. Using t-test results (we do not run TTEST2.M for you)
%
%  mentioned on the FieldTrip email list (citing Rosnow & Rosenthal, Effect
%  Sizes for Experimenting Psychologists, Canadian Journal of Experimental
%  Psychology, 2003, 57:3, 221-237.)
%
%   method = 'between_t'
%   t      = t-value
%   df     = degrees of freedom
%
% 1d. Using t-test results (we run TTEST2.M for you). Assumes equal
%     variances. If you want unequal variance, use method 1b.
%
%   method = 'between_t_run'
%   data1  = data for group 1
%   data2  = data for group 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2: WITHIN-SUBJECTS DATA (paired-sample t-test)
%
% Divide the difference of the means by the pooled standard deviation,
% corrected for correlation between the two measures (using CORR.M)
%
%   method = 'within'
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
%   method = 'single'
%   data1  = data for group 1
%   const  = single constand number to which to compare data1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also: TTEST2, CORR, STD
%

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
  
  % Divide the difference of the means by the pooled standard deviation,
  % corrected for correlation
  
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
  
elseif strcmp(method,'single')
  % 3: SINGLE-SAMPLE ONLY, comparison to a constant
  
  data1 = varargin{1};
  const = varargin{2};
  
  d = abs(mean(data1) - const) / std(data1);
  
end

end

