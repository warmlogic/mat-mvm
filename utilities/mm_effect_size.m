function [d] = mm_effect_size(method,varargin)
%MM_EFFECT_SIZE - Calculate effect size (Cohen's d)
%
%  [d] = mm_effect_size(method,varargin)
%
% A number of calculation methods are possible, and you need to carefully
% choose the method depending on the data. I HOPE these are all correct,
% but I do see a small discrepancy betwen 1a and 1b/1c.
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
% NB: There is a small difference 1a and 1b/1c. 1a agrees with G*Power 3.
%
% 1a. Divide the difference of the means by the pooled standard
%     deviation (independent-samples, equal sample sizes).
%
%   method = 'between'
%   data1  = data for group 1
%   data2  = data for group 2
%
% 1b. Using t-test results (we do not run the t-test for you)
%
%  mentioned on the FieldTrip email list (citing Rosnow & Rosenthal, Effect
%  Sizes for Experimenting Psychologists, Canadian Journal of Experimental
%  Psychology, 2003, 57:3, 221-237.)
%
%   method = 'between_t'
%   t      = t-value
%   df     = degrees of freedom
%
% 1c. Using t-test results (we run TTEST2.M for you). Assumes equal
%     variances.
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


% if strcmp(method,'between_equal')
%   % 1a. Using the root mean squared standard deviation
%   % (independent-samples, equal sample sizes)
%   
%   data1 = varargin{1};
%   data2 = varargin{2};
%   
%   d = abs(mean(data1) - mean(data2)) / ...
%     sqrt((std(data1)^2 + std(data2)^2) / 2);


if strcmp(method,'between')
  % 1a. Using the pooled standard deviation (independent-samples, equal or
  % unequal sample sizes)
  
  data1 = varargin{1};
  data2 = varargin{2};
  
  d = abs(mean(data1) - mean(data2)) / ...
    sqrt(...
    (((length(data1) - 1) * std(data1)^2) + ((length(data2) - 1) * std(data2)^2)) / ...
    (length(data1) + length(data2) - 2)...
    );
  
elseif strcmp(method,'between_t')
  % 1b. Using t-test results and pooled standard deviation
  % (independent-samples, unequal sample sizes)
  
  t = varargin{1};
  df = varargin{2};
  
  d = 2 * abs(t) / sqrt(df);
  
elseif strcmp(method,'between_t_run')
  % 1c. run the ttest2
  
  data1 = varargin{1};
  data2 = varargin{2};
  
  [h,p,ci,stats] = ttest2(data1,data2,0.05,'both');
  d = 2 * abs(stats.tstat) / sqrt(stats.df);
  
elseif strcmp(method,'within')
  % METHOD 2: WITHIN-SUBJECTS DATA (paired-sample t-test)
  
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
  
  d = abs(mean(data1) - mean(data2)) / ...
    sqrt((std(data1)^2 + std(data2)^2) - (2 * rho * std(data1) * std(data2)));
  
elseif strcmp(method,'single')
  % METHOD 3: SINGLE-SAMPLE ONLY, comparison to a constant
  
  data1 = varargin{1};
  const = varargin{2};
  
  d = abs(mean(data1) - const) / std(data1);
  
end

end

