function [stat] = statfun_glm(cfg, dat, design)

% this statfun does a regression of the design matrix on the data the
% design matrix/data are assumed to be well behaved, i.e. either the
% design contains a constant term, or the data is de-meaned
%
% This is from JM: http://bugzilla.fcdonders.nl/show_bug.cgi?id=504
%
% Note that the code is highly experimental.
%
% What you need to specify prior to calling ft_freqstatistics /
% ft_timelockstatistics is the following:
%
% cfg.statistic = 'glm';
%
% cfg.glm.statistic = 'T' 'beta' or 'R' (T-statistic, beta-weight,
% or Residuals)
%
% cfg.glm.contrast = vector specifying the linear contrast between regressors,
% should have length N, where N is the number of regressors in your
% design matrix
%
% cfg.design = design matrix NxM = number of regressors times number
% of observations
%
% The function is designed to treat each trial as an observation.
%
% I'll let you play with it, and change the status of this bug to fixed (because
% essentially it's a bug no more)
%

if ~isfield(cfg.glm, 'contrast'),  cfg.glm.contrast  = [];  end
if ~isfield(cfg.glm, 'statistic'), cfg.glm.statistic = 'F'; end
if ~isfield(cfg.glm, 'ivar0'),     cfg.glm.ivar0     = [];  end
if isfield(cfg, 'uvar') && ~isempty(cfg.uvar), error('cfg.uvar is not allowed');  end
%FIXME why?

[nvox, nobs] = size(dat);

%FIXME do preconditioning on this to speed up randomizations
dat     = standardise(dat,2);
dat     = dat./sqrt(nobs);

for k = 1:size(design,1)
  if ~all(design(k,:)==design(k,1))
    design(k,:) = standardise(design(k,:),2);
  end
end
%design  = standardise(design,2);
design  = design./sqrt(nobs);

if strcmp(cfg.glm.statistic, 'T') && ~isempty(cfg.glm.contrast),
  doT    = 1;
  doF    = 0;
  dobeta = 0;
  doR    = 0;
  lambda = cfg.glm.contrast;
elseif strcmp(cfg.glm.statistic, 'F'), 
  doT    = 0;
  doF    = 1;
  dobeta = 0;
  doR    = 0;
  ivar0  = cfg.glm.ivar0;
elseif strcmp(cfg.glm.statistic,'beta'),
  doT    = 0;
  doF    = 0;
  dobeta = 1;
  doR    = 0;
  lambda = cfg.glm.contrast;
elseif strcmp(cfg.glm.statistic, 'R'),
  doT    = 0;
  doF    = 0;
  dobeta = 0;
  doR    = 1;
end


%core business
cdx    = design*design';
invcdx = inv(cdx);
tmpdot = design*dat'; %dot product of full model
beta   = invcdx*tmpdot; %betas for full model
R      = 1-sum((beta'*design).^2,2)';

if doT,
  denom = lambda*invcdx*lambda';
  T     = (lambda*beta)./sqrt(denom.*R./(nobs-2));
  stat.stat = T(:);
end

if doF,
  if ~isempty(cfg.glm.ivar0),
    dx = design(cfg.glm.ivar0,:);

    nx = size(dx,1);
    invcdxR = inv(dx*dx'); %inverse covariance of design of reduced model
    tmpdotR = dx*dat';     %dot product of reduced model
    betaR   = invcdxR*tmpdotR; %betas of reduced model
    Rred    = 1-sum([betaR'*dx].^2,2)';
    F       = ((Rred-R)./nx)./(R./(nobs-nx-2));
  else
    error('don''t know what to do ... yet');
  end
  stat.stat = F(:);
end

if dobeta,
  if ~isempty(lambda)
    b = lambda*beta;
    stat.stat = b(:);

  else
    b = beta;
    stat.stat = b';

  end
end

if doR,
  stat.stat = R(:);
end