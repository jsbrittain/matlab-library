function csn=csn(x,mu,Sigma,alpha);%D,v,delta);
%
% Closed Skew Normal (CSN) Probability Distribution Function
%
% Univariate case
%
%

% Closed skew normal distribution
%csn=(1/cdf(0,v,delta+D*Sigma*transpose(D)))*pdf(x,mu,Sigma).*cdf(D*(x-mu),v,delta);   % Naveau (2005)
csn=2*pdf(x,0,1).*cdf(alpha*x,0,1);