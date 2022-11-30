function [Te,Fe,Fmax] = mw_coi_old( wlparam )
%
% Returns the e-folding time (e^-2) for temporal and spectral aspects.
%
% Compatible with multitaper analyses
%

% Average time-domain wavelet response using same weighting
nu = -15:0.01:15;
mag = zeros(length(nu),1);
wk = wlparam.dk/sum(wlparam.dk);
for k = (1:length(wk))
    mag = mag + wk(k)*abs(wl_morse(nu,wlparam.wlopt{1},wlparam.wlopt{2},k-1)).^2;
end;
mag = mag / max(mag);

% Determine e-folding time
left  = nu(find(diff(mag>exp(-2))>0,1,'first'));
right = nu(find(diff(mag>exp(-2))<0,1,'first'));
Te = (right-left)/2;

% Average frequency-domain wavelet response using same weighting
omega = -15:0.01:15;
mag = zeros(length(omega),1);
for k = (1:length(wk))
    mag = mag + wk(k)*abs(wlf_morse(omega.',1,wlparam.dt,wlparam.wlopt{1},wlparam.wlopt{2},k-1)).^2;
end;
mag = mag / max(mag);

% Determine e-folding time
left  = omega(find(diff(mag>exp(-2))>0,1,'first'));
right = omega(find(diff(mag>exp(-2))<0,1,'first'));
Fe = (right-left)/2/2/pi;

% Determine maximum analysis frequency
f0 = 1/wlparam.fourier_factor;
Fmax = f0/wlparam.dt/2/(Fe+f0);
