function [psi]=wlf_morlet(omega,s,dt,w0);
%function [psi]=wlf_morlet(omega,s,dt,w0);
%
% Wavelet: Morlet
% Fourier transform
%
% Common parameters
%   s       Scale
%   omega   Frequencies
%   dt      Time resolution
%
% Additional parameters
%   w0      central frequency of the mother wavelet

% Calculate wavelet response in the frequency domain
sw=s*omega;
psi0=(pi.^(-1/4)).*(omega>0).*exp(-((sw-w0).^2)/2);
norm=sqrt(2*pi*s/dt);

% Normalised wavelet
psi=norm*psi0;
