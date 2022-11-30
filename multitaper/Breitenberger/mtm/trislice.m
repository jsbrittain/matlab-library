function p=trislice(a,b,lam)
% TRISLICE - Spectrum slicing for a tridiagonal matrix
% Syntax: p=trislice(a,b,lam);
%
% Input:  a - vector of length N, containing the diagonal elements of A.
%         b - vector of length N-1, containing the off-diagonal elements of A.
%       lam - eigenvalue estimate(s). lam may be a vector.
%
% Output: p - number of eigenvalues of A greater than lam. If lam
%              is a vector, pi has n elements, which correspond to the 
%              elements of lam. 
%

a=a(:);
lam=lam(:);
N=length(a);
n=length(lam);
p1=zeros(n,N);
won=ones(n,1);
b=b(:).^2;

p1(:,1)=won*a(1)-lam;

z=find(p1(:,1)==0);
%p1(z,1)=ones(length(z),1)*eps; % replace any zeros with eps.

for i=2:N
  p1(:,i)=won*a(i)-lam-won*b(i-1)./p1(:,i-1);
%  z=find(p1(:,i)==0);
%  p1(z,i)=ones(length(z),i)*eps; % replace any zeros with eps.
end

p=zeros(1,n);
for i=1:n
  p(i)=length(find(p1(i,:)>0));
end

