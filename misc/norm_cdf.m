function cdf=normcdf(x,mu,Sigma);
%function cdf=normcdf(x,mu,Sigma);
%
% Multivariate Normal Cumulative Distribution Function
%
% UNIVARIATE ONLY AT THIS STAGE
%
%function cdf=cdf(x,mu,Sigma);

if (length(mu)~=1)
    error(' Multivariate CDF not supported at this time.');
end;

cdf=(1+erf((x-mu)/(Sigma*sqrt(2))))/2;
