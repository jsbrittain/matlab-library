function [S,c,Sk]=mtmlin(x,E,V,conf,nfft)
% Syntax: [S,c,Sk]=mtmlin(x,E,V); [S,c,Sk]=mtmlin(x,E,V,conf,nfft)
% Mtmlin produces a linearly weighted Thomson multiple-taper 
% spectral estimate of the time series x, using the
% dpss (tapers) in E and their associated eigenvalues in V.
% Confidence limits (using a chi-squared approach) are
% computed for the level 'conf', which defaults to .95.
% The averaged spectral estimate is returned in S, and the
% confidence limits are in the 2 by N matrix c. The individual
% spectral estimates are returned in Sk. 
% For details, see e. g. Thomson 1982, Park et al. 1987.,
% Percival and Walden 1993.
%
% Written by Eric Breitenberger, version date 10/1/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%

x=x(:);
N=length(x);
W=length(V)/2;
k=2*W-1; % By convention, the first 2W eigenvalues/vectors are stored 

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


% Compute the averaged estimate: simple arithmetic
% averaging is used. The Sk can also be weighted 
% by the eigenvalues, as in Park et al. Eqn. 9.;
% note that that eqn. apparently has a typo. as
% the weights should be V and not 1/V.
S=zeros(M,1);
for i=1:k
  S=S+Sk(:,i);
%  S=S+Sk(:,i)*V(i); % Park estimate
end
S=S/k;

% Calculate confidence limits
nu=2*k;
lim=(1-conf)/2;
lim=[lim 1-lim];
conf=wilhil(nu, lim); 
conf=nu./conf;
c(:,1)=conf(1)*S;
c(:,2)=conf(2)*S;

