function [b,yb] = ar_obslik(y,ar,s2)
%[b,yb] = ar_obslik(y,ar,s2)
%
% Observation probability given AR process coefficients
%
%       ar        Coefficient set
%                   ( p x Q )
%                  where 1..Q are the AR coefficient sets
%
%[b,yb] = ar_obslik(y,ar,s2)

N = length(y);
p = size(ar,1);
Q = size(ar,2);

b = zeros(Q,N);
yb = zeros(Q,N);

for t=((p+2):N)
    % (p-dimension) time-series vector
    Ft = -[y((t-1):(-1):(t-p))];
    for j=(1:Q)
        yb(j,t) = Ft*ar(:,j);
        b(j,t) = exp(-((y(t)-yb(j,t)).^2)/(2*s2(j)));% / sqrt(2*pi*s2(j));
    end;
end;
