  function [E,V] = dpsscal2(N,W)
% DPCSSCAL2 - Slepian Sequence calculation for short sequences
%  Syntax:  [E,V] = dpsscal2(N,W);
%
%  Dpsscal2 calculates discrete prolate spheroidal sequences for the
%  parameters N and W. The tridiagonal method of Slepian is used.
%  Note that W is of form 2, 5/2, 3, 7/2, 4, ... and not 2/N, 5/2N, 3/N, etc.
%
%  DPSSCAL2 is very inefficient computationally for large N. If N is greater
%  than about 200 or so (this depends on the platform), DPSSCALC should
%  be used instead. 
%
%  Returns:
%              E: matrix of dpss (N by 2W)
%              V: eigenvalue vector (2W)
%
% See Percival and Walden, Chapter 8.
%
% Written by Eric Breitenberger, version date 2/6/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%

W=W/N;
% Generate the diagonal
d=((N-1-2*(0:N-1)).^2)*.25*cos(2*pi*W);

% Generate the off-diagonals
off=(1:N-1).*(N-1:-1:1)/2;

% Generate the tri-diagonal matrix
B = diag(d) + diag(off,1) + diag(off, -1);

% Get the eigenvectors/values of B, then sort
[E, V]=eig(B);
[v, index]=sort(diag(V));
clear B V v
E=E(:,index(N:-1:1));
% Keep only the needed eigenvectors
E=E(:,1:2*N*W);

% Now calculate the desired eigenvalues by plugging
% the eigenvectors back into the defining equation:
d(1)=2*W;
d(2:N)=sin(2*pi*W*(1:N-1))./(pi*(1:N-1));

A=toeplitz(d);
V=E'*A*E;
V=diag(V);

% Normalize the eigenvectors
E=E./(ones(N,1)*sqrt(sum(E.*E)));

% Polarize symmetric dpss
d=mean(E);
for i=1:2:2*N*W
  if d(i)<0, E(:,i)=-E(:,i); end
end

% Polarize anti-symmetric dpss
for i=2:2:2*N*W
  if E(2,i)<0, E(:,i)=-E(:,i); end
end

