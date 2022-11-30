function [Te,Fe,Fmax] = wl_coi_old( wlparam )
%
% Returns the e-folding time (e^-2) for temporal and spectral aspects.
%
% NUMERICAL METHOD REQUIRES TIME/FREQ AXIS TO BE COVERED - SLOW
%
% Compatible with multitaper analyses
%

% Instigate wavelet functions
wl_wavelet  = str2func(['wl_'  wlparam.mother]);
wlf_wavelet = str2func(['wlf_' wlparam.mother]);

% Average time-domain wavelet response using same weighting
nu = (-15:0.01:15);
mag = abs(wl_wavelet(nu,wlparam.wlopt{:})).^2;
mag = mag / max(mag);

% Determine e-folding time
left  = nu(find(diff(mag>exp(-2))>0,1,'first'));
right = nu(find(diff(mag>exp(-2))<0,1,'first'));
Te = (right-left)/2;

% Average frequency-domain wavelet response using same weighting
omega = -15:0.01:15;
mag = abs(wlf_wavelet(omega.',1,wlparam.dt,wlparam.wlopt{:})).^2;
mag = mag / max(mag);

% Determine e-folding time
left  = omega(find(diff(mag>exp(-2))>0,1,'first'));
right = omega(find(diff(mag>exp(-2))<0,1,'first'));
Fe = (right-left)/2/2/pi;       % half diff then f = omega/2/pi

% Determine maximum analysis frequency
f0 = (1/wlparam.fourier_factor);
Fmax = f0/wlparam.dt/2/(Fe+f0);
