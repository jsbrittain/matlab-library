function [S,c,nu,Sk]=mtmadap(x,E,V,conf,nfft)
% Syntax: [S,c,nu,Sk]=mtmadap(x,E,V,conf,nfft);
% Mtmadap produces a Thomson adaptive multiple-taper 
% spectral estimate of the time series x, using the
% dpss (tapers) in E and their associated eigenvalues in V.
% The averaged spectral estimate is returned in S, and the
% individual estimates are returned in Sk. 
% Confidence intervals (returned in c) are computed using 
% a chi-squared approach for a confidence level specified
% in 'conf', if not specified, this defaults to .95.
% The FFTs will be zero-padded to 'nfft' points. The default
% is to use the next power of 2.
% Only real time series are supported.
% For details, see Thomson 1982, Park et al. 1987.,
% Percival and Walden 1993.
%
% Written by Eric Breitenberger, version date 10/1/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%

x=x(:);
N=length(x);
W=length(V)/2;
k=2*W-1; % By convention, the first 2W eigenvalues/vectors are stored 
V=V(1:k);

if nargin==3, conf=.95; nfft=2.^(fix(log(N-1)/log(2))+1); end
if nargin==4
  if rem(conf,1), nfft=conf; conf=.95;
  else,  nfft=2.^(fix(log(N-1)/log(2))+1); end
end

% Compute the windowed dfts and the
%   corresponding spectral estimates:
Y=zeros(nfft,k); 
Sk=zeros(nfft,k); 
for i=1:k
  wk=E(:,i).*x;
  Y(:,i)=fft(wk,nfft);
  Sk(:,i)=abs(Y(:,i)).^2;
end

% Select the proper points from fft:
if rem(nfft,2)==0, M=nfft/2+1; else M=(nfft+1)/2; end
%Y=Y(1:M,:);
Sk=Sk(1:M,:);

% Set up the iteration to determine the adaptive weights: 

sig2=x'*x/N;                % Power
S=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate
Stemp=zeros(M,1);
S1=zeros(M,1);

% Set tolerance for acceptance of spectral estimate:
% The algorithm converges so fast that results are
% usually 'indistinguishable' after about three iterations.

% This version uses the equations from P&W pp 368-370

tol=.0005*sig2/M;
i=0;
a=sig2*(1-V);

% Do the iteration:
while sum(abs(S-S1)/M)>tol
  i=i+1;
  % calculate weights
  b=(S*ones(1,k))./(S*V'+ones(M,1)*a'); 
  % calculate new spectral estimate
  wk=(b.^2).*(ones(M,1)*V');
  S1=sum(wk'.*Sk')./ sum(wk');
  S1=S1';
  Stemp=S1; S1=S; S=Stemp;  % swap S and S1
end

nu=2*sum(wk').^2./sum(wk'.^2);
lim=(1-conf)/2;
lim=[lim 1-lim];
c=wilhil(nu, lim);
c(:,1)=nu'.*S./c(:,1);
c(:,2)=nu'.*S./c(:,2);

