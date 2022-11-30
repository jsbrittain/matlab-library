function [x,x1,x2] = modelS2( z, target, clamptrials )
% Multi-rate model

% Specify parameters
Af = z(1);
Bf = z(2);
As = z(3);
Bs = z(4);

% Reserve memory
x = zeros(1,length(target));
x1 = x; x2 = x; e = x;

% Target
F = target;

% Recurse trials
for n = (2:length(target))
    % Update parameters
    x1(n) = Af*x1(n-1) + Bf*e(n-1);
    x2(n) = As*x2(n-1) + Bs*e(n-1);
    x(n)  = x1(n) + x2(n);
    % Error-clamp trial, measure adaptation
    if ( clamptrials(n) )
        F(n) = x(n);
    end
    e(n) = F(n) - x(n);
end
