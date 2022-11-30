function [T,F] = ht_coi( wlparam )
%
% Returns the e-folding time (e^-2) for temporal and spectral aspects.
%
% Compatible with multitaper analyses
%

% Normalise weights
wk = wlparam.E/sum(wlparam.E);

%% Temporal

dt = 0.001;
nu = (-15:dt:15);             % Time vector
N = length(nu);

scale = wlparam.wlopt{1};
K = length(wk);
h = ht_hermites(N,K,1/dt,scale);

% Energy profile of composite (time-domain) analytic function
mag = zeros(length(nu),1);      % Magnitude
for k = (1:length(wk))
    mag = mag + wk(k)*abs( h(:,k) ).^2;
end;
mag = mag / max(mag);

% Determine e-folding time
left  = nu(find(diff(mag>exp(-2))>0,1,'first'));
right = nu(find(diff(mag>exp(-2))<0,1,'first'));
T = (right-left)/2;

%% Spectral

% High-resolution frequency response
H = []; tt = (1:N);
freqs = (-100:0.1:100);
for k = (1:K)
    fh = [];
    for n = (1:length(freqs))
        fh(n) = sum( (h(tt,k).').*exp(-1i*2*pi*freqs(n)*tt/wlparam.rate) );
    end;
    H(:,k) = abs(fh).^2;
end;
mag = sum( abs(H).^2, 2 );
mag = mag/max(mag(:));

% Determine e-folding time
left  = freqs(find(diff(mag>exp(-2))>0,1,'first'));
right = freqs(find(diff(mag>exp(-2))<0,1,'first'));
F = (right-left)/2;
