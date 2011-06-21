function [stat] = statistics_glm(cfg, data, design);

% This is from RO: http://bugzilla.fcdonders.nl/show_bug.cgi?id=504
%
% This is a helper function that computes a statistical test based on the
% general linear mdel (GLM), given the biological data and the design
% matrix. 
% 
% Required configuration option:
% cfg.method  = 't_test' or 'glm'
%
% The configuration should also contain: 
% cfg.alpha   = 0.05
% cfg.repmeas = 'no' (default) or 'yes'
%
% This function is called by STATISTICS_WRAPPER, which again is called by
% either TIMELOCKSTATISTICS, FREQSTATISTICS or SOURCESTATISTICS.

% Copyright (C) 2005, Markus Bauer
%
% $Log: statistics_glm.m,v $
% Revision 1.6  2008/11/25 13:27:57  estmee
% Documentation update
%
% Revision 1.5  2008/09/26 15:27:14  sashae
% checkconfig: checks if the input cfg is valid for this function
%
% Revision 1.4  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.3  2006/06/07 12:58:36  roboos
% added default value for cfg.alpha
%
% Revision 1.2  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.1  2005/06/03 10:44:50  roboos
% new implementation by marbau, also contains t_test
%

ft_defaults

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'required',    {'method'});

% set the defaults
if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end
  
if ~isfield(cfg,'repmeas')
  warning('assume independent samples - no repeated measurement design');
  cfg.repmeas = 'no';
end;

sizreplic = size(design);
[numobs,numrepl] = size(data);

if (strcmp(cfg.method,'t_test'))

  %find out whether this is appropriate design for t-test, i.e. only two
  %conditions
  contrast = unique(design);
  if(length(contrast)>2)  error('contrasts specified in a wrong way !'); end;

  %find homogeneous conditions and average those
  relconditions{1} = find(design==(contrast(1)));numrepl1=length(relconditions{1});
  relconditions{2} = find(design==(contrast(2)));numrepl2=length(relconditions{2});
  condition1 = data(:,relconditions{1});  %/;
  condition2 = data(:,relconditions{2});  %/length(relconditions{2});

  contrastdiff = [];
  if (contrast(2)>contrast(1)), contrast(2)=1;contrast(1)=-1;  else contrast(1)=-1;contrast(2)=1;    end;

  if strcmp(cfg.repmeas,'yes')
    if (numrepl1~=numrepl2)
      error('unequal samples - no paired t-test possible');
    end;
    contrastdiff = condition1*contrast(1) + condition2*contrast(2);
    contrastdiffmean = (nanmean(contrastdiff'))'; %individual subjects mean value over all conditions, for each tpmt,chan
    stderrmean = (nanstd(contrastdiff'))';
    df = numrepl1-1;
  else
    contrastdiffmean = ((nanmean(condition1'))')*contrast(1) + ((nanmean(condition2'))')*contrast(2);
    stderrmean = (((nanstd(condition1'))')*numrepl1+((nanstd(condition2'))')*numrepl2)/(numrepl1+numrepl2-2)*sqrt((1/numrepl1)+(1/numrepl2));
    df = numrepl1+numrepl2-2;
  end;

  %calculate the t-test here
  stat.t_value = contrastdiffmean./stderrmean;
  stat.p_value(:,:) = tcdf(stat.t_value,repmat(df,[numobs,1]));
  stat.signif = zeros(size(stat.p_value));

  %include a conditional statement for two-tailed testing here

  signif = find(stat.p_value>(1-cfg.alpha));
  stat.signif(signif) = 1;
  
elseif (strcmp(cfg.method,'glm'))

  %in order to prevent the variable for testing the intercept to corrupt
  %result of the R^2 for criterion with predictor variable of interest,
  %the total mean of each observation is subtracted - such that the
  %obligatory intercept variable won't contribute !
  %     data = data - repmat(mean(data,2),[1,size(data,2)]);
  data = data' - repmat(nanmean(data'),[size(data,2),1]);
  data = data';

  if strcmp(cfg.repmeas,'yes')
    if size(design,1) > 2
      error('repeated measurement design currently not supported for more than one predictor variable');
    end;
    factorsteps = unique(design(:));
    nconditions = factorsteps;
    for i=1:nconditions
      dumindx(:,i) = find(design==factorsteps(i))';
    end;
    for i=1:size(dumindx,1)
      meanrt(:,i) = (nanmean(data(:,dumindx(i,:))'))';
    end;
    dummy = ones(size(data,2),1);
    for i=1:nconditions
      dummy(dumindx(:,i)) = meanrt(:,i);
    end;
    meanrt = dummy';
    regressor = [design;meanrt;ones([1,size(design,2)])];
  else
    regressor = [design;ones([1,size(design,2)])];
  end;
  %warning off MATLAB:divideByZero;
  for obs=1:numobs
    %throw away missing values - df correction automatically taken care of
    datadum = data(obs,:);datadum = datadum(:);
    q = find(isnan(datadum));
    regressdum = regressor';
    datadum(q,:)=[];regressdum(q,:)=[];
    clear q;
    for zz = 1:size(regressdum,2)
      q = find(isnan(regressdum(:,zz)));
      datadum(q,:)=[];regressdum(q,:)=[];
    end;
    clear q;

    if ~isempty(datadum)
      [b,bint,r,rint,stats] = regress(datadum,regressdum,cfg.alpha);
      beta = b'.*std(regressdum,0,1)/(std(datadum,0,1));
      stat.b_value(obs,:) = b;
      stat.beta_value(obs,:) = beta;
      stat.R_value(obs,:) = stats(1);
      stat.F_value(obs,:) = stats(2);
      stat.p_value(obs,:) = stats(3);
    else
      stat.b_value(obs,:) = [NaN,NaN];
      stat.beta_value(obs,:) = [NaN,NaN];
      stat.R_value(obs,:) = NaN;
      stat.F_value(obs,:) = NaN;
      stat.p_value(obs,:) = NaN;
    end;
  end;
  if strcmp(cfg.repmeas,'yes')
    stat.b_value = stat.b_value(:,1:(end-2));
    stat.beta_value = stat.beta_value(:,1:(end-2));
  else
    stat.b_value = stat.b_value(:,1:(end-1));
    stat.beta_value = stat.beta_value(:,1:(end-1));
  end;
end;

