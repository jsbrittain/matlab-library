function [p1,p2,stats] = analysisDriftPermJointProb( metric, nulldistr, varargin )
%[p1,p2] = analysisDriftPermJointProb( metric, nulldistr )
%
% Empirical test based on joint distributions
%
%  If multivariate data provided then a PCA is applied to the null
%  distributions and marginal probabilities are multiplied.
%
%[p1,p2] = analysisDriftPermJointProb( metric, nulldistr, decorrelation_method )

% Check for empty inputs
if (isempty(metric))
    p1 = NaN;
    p2 = NaN;
    return;
end;

% Input parser
p = inputParser;
addRequired(p,'metric',@(x) (isnumeric(x) && isvector(x)));
addRequired(p,'nulldistr',@(x) (isnumeric(x) && (size(metric,1)==size(x,1))));
addParameter(p,'side','right',@ischar);
parse(p,metric,nulldistr,varargin{:});
% Unpack optional parameters
p1side = lower( p.Results.side );

% Sort distribution for each variable
snd = nan(size(nulldistr));
for k = (1:size(nulldistr,1))
    snd(k,:) = sort(nulldistr(k,:));
end;

% One-sided test -- numerical
switch ( p1side )
    case { 'left' },
        %p1 = 1-mvtcdf(metric-mean(nulldistr,2),cov(nulldistr'),size(nulldistr,2)-1); return;
        jcdf = repmat(metric',size(nulldistr,2),1) >= (snd');
    case { 'right' },
        %p1 = mvtcdf(metric-mean(nulldistr,2),cov(nulldistr'),size(nulldistr,2)-1); return;
        jcdf = repmat(metric',size(nulldistr,2),1) <= (snd');
    otherwise
        error('Unknown side specified for test -- should be left/right');
end;
jcdf = sum(jcdf,2)==size(jcdf,2);                           % <-- AND for multivariate test
p1 = sum(jcdf)/numel(jcdf);
if ( p1 == 0 )
    p1 = 1/numel(jcdf);
end;

% % Two-sided test -- numerical (assumes independence; does not match Hotelling's T^2)
% jR = repmat(abs(metric)',size(nulldistr,2),1) <= (snd');
% pR =   sum( jR(:) )/numel(jR);
% jL = repmat(-abs(metric)',size(nulldistr,2),1) <= (snd');
% pL = 1-sum( jL(:) )/numel(jL);
% p2 = pL + pR;
% if (p2==0)
%     p2 = 1/numel(jcdf);
% end;

% Hotelling's T2 (parametric two-sided test)
n = size(nulldistr,2); p = size(nulldistr,1);
popcov = n*cov(nulldistr.');
t2 = n*(mean(nulldistr,2)-metric)'*inv(popcov)*(mean(nulldistr,2)-metric);
fstat = t2*(n-p)/(p*(n-1));
p2 = 1-fcdf( fstat, p, n-p );        % <-- t-stat ( unknown covar )
%p2 = 1-chi2cdf(t2,p);                              % <-- chi-2 ( known covar )

% Summarise T2 test statistics
stats.p2.type = 'Hotelling''s T^2';
stats.p2.t2 = t2;
stats.p2.n = n;
stats.p2.p = p;
stats.p2.pval = p2;
stats.p2.str = sprintf('%s, T2 = %.3f, F(%g,%g) = %.3f, p = %.3f',stats.p2.type,t2,p,n-p,fstat,p2);
