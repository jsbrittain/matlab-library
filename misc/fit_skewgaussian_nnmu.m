function [c,gof]=fit_skewgaussian_nnmu(tt,y)
%
%
% Non-negative mu variant
%
% Parameters (alphabetical)
%       A       Amplitude scaling (->pdf)
%       mu      Mean
%       sd      Std.dev.
%       skew    Alpha-skew level
%

% Sample estimates for start point
correction=sum(y);          %
p=y/correction;             % Correct for sampling interval
sp_mu=sum(tt.*p);
sp_sd=sum((tt-sp_mu).^2.*p);
sp_alpha=sum((tt-sp_mu).^3.*p);
startpoint=[1 sp_mu, sp_sd, sp_alpha];

s = fitoptions( 'Method','NonlinearLeastSquares', ...
                'Lower',[0, -Inf, 0 -Inf],          ...
                'Upper',[Inf, Inf, Inf Inf],          ...
                'Startpoint',startpoint           ...
              );

strnormpdf='(1/(sd*sqrt(2*pi)))*exp(-(x-mu)^2/(2*sd^2))';  % PDF-norm
strnormcdf='(1/2)*(1+erf((skew*(x-mu)/sd)/sqrt(2)))';  % CDF-norm with skew-alpha
strskewgauss=['2*' strnormpdf '*' strnormcdf '*A'];

f=fittype(strskewgauss,'options',s);
[c,gof]=fit(tt,y,f);
