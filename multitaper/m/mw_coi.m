function [Te,Fe,Fmax] = mw_coi( wlparam )
%
% Returns the e-folding time (e^-2) for temporal and spectral aspects.
%
% Compatible with multitaper analyses
%


[Te,Fe,Fmax] = mw_coi_old( wlparam );

return;






% Instigate wavelet functions
wl_wavelet  = str2func(['wl_'  wlparam.mother '_composite_energy']);
wlf_wavelet = str2func(['wlf_' wlparam.mother '_composite_energy']);

% Average time-domain wavelet response using same weighting
wk = wlparam.dk/sum(wlparam.dk);
K = length(wk);
fnc = @(x) abs( (abs(wl_wavelet(x,wlparam.wlopt{1},wlparam.wlopt{2},wk)).^2)/(abs(wl_wavelet(0,wlparam.wlopt{1},wlparam.wlopt{2},wk)).^2) - exp(-2) );
Te = wlparam.dt*fminsearch( fnc, 0 );       % Scale since single search exposure for wl_morse returns for dt=1

% Average time-domain wavelet response using same weighting
nu = -15:0.01:15;
mag = zeros(length(nu),1);

for k = (1:length(wk))
    mag = mag + wk(k)*abs(wl_morse(nu,wlparam.wlopt{1},wlparam.wlopt{2},k-1)).^2;
end;
mag = mag / max(mag);

% Determine e-folding time
left  = -nu(find(diff(mag>exp(-2))>0,1,'first'));
right =  nu(find(diff(mag>exp(-2))<0,1,'first'));
Te = (right+left)/2;

% Average frequency-domain wavelet response using same weighting
omega = -15:0.01:15;
mag = zeros(length(omega),1);
for k = (1:length(wk))
    mag = mag + wk(k)*abs(wlf_morse(omega.',1,wlparam.dt,wlparam.wlopt{1},wlparam.wlopt{2},k-1)).^2;
end;
mag = mag / max(mag);

% Determine e-folding time
left  = -omega(find(diff(mag>exp(-2))>0,1,'first'));
right =  omega(find(diff(mag>exp(-2))<0,1,'first'));
Fe = (right+left)/2;

% Determine maximum analysis frequency
f0 = (1/wlparam.fourier_factor);
Fmax = f0/wlparam.dt/2/(Fe+f0);
