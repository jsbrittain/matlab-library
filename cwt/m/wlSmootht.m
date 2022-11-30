function Wbar=wlSmootht(W,trange,rate);
%function Wbar=wlSmootht(W,trange,rate);
%
% Wavelet smoothing in the time domain
%

% Smoothing parameters
tmin=trange(1);
tmax=trange(2);
tt=[(tmin*rate/1000):(tmax*rate/1000)];

% Perform smoothing between min and max ranges
Wbar=mean(W(:,tt),2);
