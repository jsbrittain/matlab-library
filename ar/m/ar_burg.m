function [a,s2] = ar_burg( y, p )
%function [a,s2] = ar_burg( y, p )
%
%
% Reference
%   Frequency Warped Burg's Method for AR-Modeling
%   Roth, Kauppinen, Esquef, Valimaki
%   IEEE Workshop on Appls of Sig Proc to Aud and Acoust. Oct 19-22, 2003.
%
%function [a,s2] = ar_burg( y, p )


% JSB Addition - Variance
a11 = -mean(y(1:end-1).*y(2:end))/mean(y.^2);
s2 = (1-abs(a11)^2)*mean(y.^2);

f(:,1) = y;
b(:,1) = y;
k = zeros(1,p+1);
a = nan(p+1,p+1); a(1,:) = 1;
for l = (1:p)
    % Reflection coefficients
    k(l+1) = -2*sum(f((l+1):end,l).*b(l:(end-1),l))/(sum((f((l+1):end,l)).^2+(b(l:(end-1),l)).^2));
    % Autoregressive coefficients
    for m = (1:(l-1))
        a(m+1,l+1) = a(m+1,l) + k(l+1)*a(l-m+1,l);
    end;
    a(l+1,l+1) = k(l+1);
    % Forward and backward prediction errors
    f(:,l+1) = f(:,l) + k(l+1)*b([end 1:end-1],l);
    b(:,l+1) = b([end 1:end-1],l) + k(l+1)*f(:,l);
    
    % JSB Addition - Variance
    s2(l+1) = (1-abs(a(l+1,l+1))^2)*s2(l);
end;
