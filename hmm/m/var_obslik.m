function [b,yb] = var_obslik(y,A,sigma)
%[b,yb] = var_obslik(y,ar,s2)
%
% Observation probability given AR process coefficients
%
%       A          Coefficient set
%                   (  )
%                  where 1..Q are the AR coefficient sets
%
%[b,yb] = var_obslik(y,A,sigma)

N = size(y,2);
m = size(A,1);
p = size(A,3);
Q = size(A,4);

b = zeros(Q,N);
yb = zeros(Q,N);

for j=(1:Q)
    q{j} = reshape( A(:,:,:,j), m, m*p ).';
    q{j} = q{j}(:);
end;

for t=((p+2):N)
    % (p-dimension) time-series vector
    Ft = kron( eye(m), reshape( -(y(:,(t-1):(-1):(t-p))), 1, m*p ) );
    for j=(1:Q)
        yhat=Ft*q{j};
        b(j,t) = exp( -0.5*(y(:,t)-yhat)'*inv(sigma)*(y(:,t)-yhat) );
        %b(j,t) = b(j,t) / ( (2*pi)^(m/2) * det(sigma)^(1/2) );
    end;
end;
