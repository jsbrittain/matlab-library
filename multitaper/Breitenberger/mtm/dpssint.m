  function [En,V] = dpssint(N, W, M, int)
% Syntax: [En,V]=dpssint(N,W); [En,V]=dpssint(N,W,M,'spline');
%  Dpssint calculates discrete prolate spheroidal
%  sequences for the parameters N and W. Note that
%  W is normally 2, 5/2, 3, 7/2, or 4 - not i/N. The 
%  dpss are interpolated from previously calculated 
%  dpss of order M (128, 256, 512, or 1024). 256 is the 
%  default for M. The interpolation can be 'linear' 
%  or 'spline'. 'Linear' is faster, 'spline' the default.
%  Linear interpolation can only be used for M>N. 
%  Returns:
%              E: matrix of dpss (N by 2W)
%              V: eigenvalue vector (2W)
% 
% Errors in the interpolated dpss are very small but should be 
% checked if possible. The differences between interpolated
% values and values from dpsscalc are generally of order
% 10ee-5 or better. Spline interpolation is generally more
% accurate. Fractional errors can be very large near
% the zero-crossings but this does not seriously affect
% windowing calculations. The error can be reduced by using
% values for M which are close to N.
%
% Written by Eric Breitenberger, version date 10/3/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%

if     nargin==2,
  M=256; int='spline';
elseif nargin==3,
  if isstr(M), int=M; M=256; 
  else, int='spline';, end
end

if int=='linear' & N>M
  error('Linear interpolation cannot be used for N>M. Use splining instead.')
end

if rem(W,1)==0
  filename=['s' int2str(M) '_' int2str(W) '.mat']; 
else
  filename=['s' int2str(M) '_' int2str(10*W) '.mat']; 
end

if exist(filename)==2
  eval(['load ' filename])
else
  filename
  error('Requested DPSS file not found.')
end


k=2*W;  % Use only first k values
E=E(:,1:k);
V=V(1:k);
En=zeros(N,k);
x=1:M;

% The scaling for the interpolation:
% This is not necessarily optimal, and 
% changing s can improve accuracy.
 
s=M/N;
midm=(M+1)/2;
midn=(N+1)/2;
delta=midm-s*midn;
xi=linspace(1-delta, M+delta, N);

% Interpolate from M values to N
% Spline interpolation is a bit better,
% but takes about twice as long.
% Errors from linear interpolation are 
% usually smaller than errors from scaling.

if int=='linear'
  for i=1:k, En(:,i)=interp1(x,E(:,i),xi); end
elseif int=='spline'
  for i=1:k, En(:,i)=spline(x,E(:,i),xi)'; end
end

% Re-normalize the eigenvectors
En=En./(ones(N,1)*sqrt(sum(En.*En)));

