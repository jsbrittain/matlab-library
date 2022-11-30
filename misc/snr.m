function snrdb=snr(signal,noise);
%function snrdb=snr(signal,noise);
%
% Function to compute the signal-to-noise ratio in dB using a total power
% ratio.
%
% Input parameters
%       signal      Time-series of clean signal
%       noise       Time-series of additive noise
%
% Output parameter
%       snrdb       Signal-to-noise ratio in dB
%
%function snrdb=snr(signal,noise);

% Calculate signal power
Ps=mean(abs(signal).^2);
Pn=mean(abs(noise).^2);

% Determine SNR in dB
snrdb=10*log10(Ps./Pn);
