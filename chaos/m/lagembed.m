function X=lagembed(x,l,m)
%function X=lagembed(x,l,m)
%
% Multivariate lag-embedding
%
% Input Parameters
%    x      Input column vectors (N x k)
%    l      Lag-embed time
%    m      Embedding dimension
%
%function X=lagembed(x,l,m)

dim=size(x,2);                      % Input dimension
N=(length(x)-(m-1)*l);              % Length of embedding vector
X=zeros(N,m,dim);                   % Allocate memory
for k=(1:dim)
    for i=(1:N)
        X(i,:,k)=x(i+(0:l:(m-1)*l),k);
    end;
end;
