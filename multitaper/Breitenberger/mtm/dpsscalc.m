  function [E,V] = dpsscalc(N, W)
% DPSSCALC - calculate the Discrete Prolate Spheroidal (Slepian) Sequences 
%  Syntax:  [E,V]=dpsscalc(N, W);
%
%  Dpsscalc calculates discrete prolate spheroidal sequences for the 
%  parameters N and W, where N is the sequence length.
%  Note that W is of form 2, 5/2, 3, 7/2, 4, ... and not 2/N, 5/2N, 3/N, etc.
%
%  Returns:
%              E: matrix of dpss (N by 2W)
%              V: eigenvalue vector (2W)
%
% DPSSCALC uses the tridiagonal method of Slepian to get the eigenvectors,
% and then uses these eigenvectors with the sinc matrix to get the 
% eigenvalues. Spectrum slicing is applied to the tridiagonal matrix to get 
% the eigenvalues, then partial recursion and reflection are used to get the 
% eigenvectors. Finally, the eigenvalues of the sinc matrix are computed.
%  
% Note that DPSSCAL2 will be faster than DPSSCALC for N<200 or so.
%
% See Percival and Walden, Chapter 8.
%
% Written by Eric Breitenberger, version date 2/6/96.
% Please send comments and suggestions to eric@gi.alaska.edu

n=2*W;   % The last eigenvalue/vector pair is marginally useful - if you
         % really want to tweak for speed, set n=2*W-1;
W=W/N;

% Generate the diagonal of the tridiagonal matrix:
a=((N-1-2*(0:N-1)).^2)*.25*cos(2*pi*W);

% Generate the off-diagonals:
b=(1:N-1).*(N-1:-1:1)/2; 

% An empirical fit to the eigenvalue distribution is used to narrow down
% the initial search range:
offset=(n-1)/10*(22.5*(W)+15);
hi=N^2/4;  % Upper Gerschgorin limit
lo=hi-offset ;

% Set up to localize eigenvalues
K=2*n;                % number of intervals for initial search                   
int=(hi-lo)/K;        % split Gerschgorin interval into K intervals of width 'int'.
lam=lo+int*(0:K-1);   % initial search values, in increasing order

% First iteration: find an interval containing all n eigenvalues:
keepgoing=1;
while keepgoing
  p=trislice(a,b,lam);
  i=find(p>n);
  if p(1)==n % range containing all eigenvalues has been found - exit loop.
    keepgoing=0;
  elseif p(1)<n  % range is too small - increase it by int/2
    lo=lam(1)-int/2;
    int=(hi-lo)/K;
    lam=lo+int*(0:K-1);
  else           % range is too big - reduce it.
    lo=lam(i(length(i)))+int/2;
    int=(hi-lo)/K;
    lam=lo+int*(0:K-1);
  end 
end

% Second iteration - reduce to intervals containing one eigenvalue each:
% This loop is just for insurance - eigenvalues *should* be localized already.
d=diff([p 0]);
loc=find(d==-1); % locate the eigenvalue

while length(loc)~=n % 
  K=2*K;
  int=(hi-lo)/K;
  lam=lo+int*(0:K-1);
  p=trislice(a,b,lam);
  d=diff([p 0]);
  loc=find(d==-1); 
end

% Third iteration: nail down the eigenvalues by parallel M-section:
% Choice of M is subjective, but M=10 seems fairly good.
M=10;
p=p(loc);
m=0:M-1;
tol=M*2*eps;     % Tolerance

while int>tol
  int=int/M;
  temp=lam(loc);
  for j=1:n
    lam((j-1)*M+1:j*M)=temp(j)+(int*m);
  end
  p=trislice(a,b,lam);
  d=diff([p 0]);
  loc=find(d==-1); 
end
lam=lam(loc);

% Now calculate the eigenvectors by recursion: To reduce accumulated
% error and speed computation, only the first half is computed directly.
% The second half is computed by reflection.

if rem(N,2)==0
  midpt=N/2;
  righthalf=midpt+1:N;
  lefthalf=midpt:-1:1;
else
  midpt=(N-1)/2+1;
  righthalf=midpt+1:N;
  lefthalf=midpt-1:-1:1;
end

% calculation of dpss:
E=zeros(N,n);
E(1,:)=ones(1,n);
E(2,:)=(lam-a(1))/b(1);
for i=3:midpt
  E(i,:)=((lam-a(i-1)).*E(i-1,:) - b(i-2)*E(i-2,:))/b(i-1);
end

% sort lam and E:
lam=lam(n:-1:1);
E=E(:,n:-1:1);

% fill in the symmetric dpss:
E(righthalf,1:2:n)=E(lefthalf,1:2:n);

% fill in the anti-symmetric dpss:
E(righthalf,2:2:n)=-E(lefthalf,2:2:n);


% Normalize the eigenvectors
E=E./(ones(N,1)*sqrt(sum(E.*E)));

% Polarize symmetric dpss
d=mean(E);
for i=1:2:2*N*W-1
  if d(i)<0, E(:,i)=-E(:,i); end
end

% Polarize anti-symmetric dpss
for i=2:2:2*N*W-1
  if E(2,i)<0, E(:,i)=-E(:,i); end
end

% Now calculate the desired eigenvalues by plugging the
% eigenvectors back into the defining sinc matrix:
d(1)=2*W;
d(2:N)=sin(2*pi*W*(1:N-1))./(pi*(1:N-1));

% Two methods can be used here to compute the eigenvalues:
% on my system, (P-100 with 16MB), N>700 or so starts heavy swapping
% if the direct multiplication is used, so the looping is faster for
% these large N. The parameter can be changed depending on the
% platform and available memory.

if N>700 % use less memory intensive looping
  S=zeros(N,n);
  for i=1:N
    row=d([i:-1:2 1:N-i+1]);
    S(i,:)=sum((row'*ones(1,n)).*E);
  end
  V=diag(E'*S);
else   % more memory intensive but faster
  A=toeplitz(d);
  V=diag(E'*A*E);
end
