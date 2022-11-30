function [a,s2,wf] = arw_burg( y, p, lambda )
%function [a,s2,wf] = arw_burg( y, p, lambda )
%
%
% Reference
%   Frequency Warped Burg's Method for AR-Modeling
%   Roth, Kauppinen, Esquef, Valimaki
%   IEEE Workshop on Appls of Sig Proc to Aud and Acoust. Oct 19-22, 2003.
%
%function [a,s2,wf] = arw_burg( y, p, lambda )

if (~exist('lambda'))
    lambda = [];
end;
if (isempty(lambda))
    lambda = 0.723;
    disp([ ' Lambda = ' num2str( lambda ) ]);
end;

% Initialise variables
f(:,1) = y;
b(:,1) = y;
wb = zeros(length(y),p+1);
wk = zeros(1,p+1);
a = nan(p+1,p+1); a(1,:) = 1;

% JSB Addition - Variance
a11 = -mean(y(1:end-1).*y(2:end))/mean(y.^2);
s2 = (1-abs(a11)^2)*mean(y.^2);

% Increment model order
for l = (1:p)
    % Warped backwards prediction error
    for n = (l+1:length(y))
        wb(n,l+1) = b(n-1,l) - lambda*(b(n,l)-wb(n-1,l+1));
    end;
    
    % Warped reflection coefficients
    wk(l+1) = -2*sum(f((l+1):end,l).*wb((l+1):end,l+1))/(sum((f((l+1):end,l)).^2+(wb((l+1):end,l+1)).^2));
    % Autoregressive coefficients
    for m = (1:(l-1))
        a(m+1,l+1) = a(m+1,l) + wk(l+1)*a(l-m+1,l);
    end;
    a(l+1,l+1) = wk(l+1);
    % Forward and backward prediction errors
    f(:,l+1) = f(:,l) + wk(l+1)*wb(:,l+1);
    b(:,l+1) = wb(:,l+1) + wk(l+1)*f(:,l);
    
    % JSB Addition - Variance
    s2(l+1) = (1-abs(a(l+1,l+1))^2)*s2(l);
end;
% Frequency map
wf = inline(['atan2((1-' num2str(lambda) '^2)*sin(x),((1+' num2str(lambda) '^2)*cos(x)-2*' num2str(lambda) '));']);
