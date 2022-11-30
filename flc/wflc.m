function [y,hy,w0] = wflc( s, omega0, M, mu0, mu1, muhat, dt )
%function [y,hy] = wflc( s, omega0, M, mu0, mu1, muhat, dt )
%
% Weighted Fourier Linear Combiner (batch implementation)
%
% Input parameters
%   s       Vector input signal
%   omega0  Base frequency (fixed)
%   M       Harmonics
%   mu      Adaptation coefficient (scaled by sampling rate)
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
%function [y,hy] = wflc( s, omega0, M, mu0, mu1, muhat, dt )

% WFLC
K = length(s);
x = zeros(2*M,K); w = x; what = w;
w0 = zeros(1,K); e = w0; ehat = e; y = e;

w0(1) = omega0*dt;
sumw0 = w0(1);

for k = (1:K)
    
    sumw0 = sumw0 + w0(k);
    for r = 1:2*M
        if (r<=M)
            x(r,k) = sin(r*sumw0);
        else
            x(r,k) = cos((r-M)*sumw0);
        end;
    end;
    
    e(k) = s(k) - w(:,k)'*x(:,k);
    
    w0sum = 0;
    for r = 1:M
        w0sum = w0sum + r*( w(r,k)*x(M+r,k) - w(M+r)*x(r,k) );
    end;
    w0(k+1) = w0(k) + 2*mu0*e(k)*w0sum;
    
    w(:,k+1) = w(:,k) + 2*mu1*x(:,k)*e(k);
    
    % Amplitude tracking correction (otherwise delayed)
    y(k) = what(:,k)'*x(:,k);
    ehat(k) = s(k) - y(k);
    what(:,k+1) = what(:,k) + 2*muhat*x(:,k)*ehat(k);
    
end;
w=w(:,1:end-1);
what=what(:,1:end-1);

hy = (what(1+M,:) + 1i*what(1,:)) ./ (x(1+M,:) + 1i*x(1,:));
