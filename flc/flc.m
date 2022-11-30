function [y,hy] = flc( s, omega0, M, mu, dt )
%function [y,hy] = flc( s, omega0, M, mu, dt )
%
% Fourier Linear Combiner (batch implementation)
%
% Input parameters
%   s       Vector input signal
%   omega0  Base frequency (fixed)
%   M       Harmonics
%   mu      Adaptation coefficient (for instance, 0.01)
%
% Output
%   y       Estimated signal
%   hy      Complex estimate of tracked signal at BASE frequency ONLY
%
% Reference
%   Veluvolu et al. (2007) Bandlimited Multiple Fourier Linear Combiner for
%   Real-time Tremor Compensation. Proc 29th Ann Int Conf IEEE EMBS, Aug
%   23-26, 1007.
%
%function [y,hy] = flc( s, omega0, M, mu )

% FLC
K = length(s);
x = zeros(2*M,K); w = x;
y = nan(size(s));
for k = (1:K)
    for r = 1:2*M
        if (r<=M)
            x(r,k) = sin(r*omega0*k*dt);
        else
            x(r,k) = cos((r-M)*omega0*k*dt);
        end;
    end;
    y(k) = w(:,k)'*x(:,k);
    e(k) = s(k) - y(k);
    w(:,k+1) = w(:,k) + 2*mu*x(:,k)*e(k);
end;
w=w(:,1:end-1);

hy = (w(1+M,:) + 1i*w(1,:)) ./ (x(1+M,:) + 1i*x(1,:));
hy = reshape(hy,size(y));
