function [y,hy] = blflc( s, frange, G, mu, dt )
%function [y] = blflc( ... )
%
% Band-Limited Fourier Linear Combiner (batch implementation)
%
% Input parameters
%   s       Vector input signal
%   omega0  Base frequency (fixed)
%   M       Harmonics
%   mu      Adaptation coefficient (for instance, 0.01)
%
% Output
%   y       Estimated signal
%   hy       ### Not implemented ###  Complex estimate of tracked signal at BASE frequency ONLY
%
% Reference
%   Veluvolu et al. (2007) Bandlimited Multiple Fourier Linear Combiner for
%   Real-time Tremor Compensation. Proc 29th Ann Int Conf IEEE EMBS, Aug
%   23-26, 1007.
%
%function [y] = blflc( ... )

f0 = frange(1);
f1 = frange(2);
L = (f1-f0)*G;

% FLC
K = length(s);
x = zeros(2*L,K); w = x;
y = zeros(1,K); e = y;
for k = (1:K)
    for r = 1:2*L
        if (r<=L)
            x(r,k) = sin(2*pi*(f0+(r-1)/G)*k*dt);
        else
            x(r,k) = cos(2*pi*(f0+((r-L)-1)/G)*k*dt);
        end;
    end;
    y(k) = w(:,k)'*x(:,k);
    e(k) = s(k) - y(k);
    w(:,k+1) = w(:,k) + 2*mu*x(:,k)*e(k);
    
    [~,ix(k)] = max( abs(x((L+1):end,k) + 1i*x(1:L,k)) );
end;

for k = 1:K
    hy(k) = (w(ix(k)+L,k) + 1i*w(ix(k),k)) ./ (x(ix(k)+L,k) + 1i*x(ix(k),k));
end;
