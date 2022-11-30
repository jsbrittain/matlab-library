function [psi]=wlf_dog(omega,s,dt,m);
%function [psi]=wlf_dog(omega,s,dt,m);
%
% Wavelet: DOG
% Fourier transform
%
% Common parameters
%   s       Scale
%   omega   Frequencies
%   dt      Time resolution
%
% Additional parameters
%   m       derivative (of Gaussian)

% Calculate wavelet response in the frequency domain
sw=s*omega;
psi0=((-(i^m))/sqrt(gamma(m+.5)))*(sw.^m).*exp(-(sw.^2)/2);
norm=sqrt(2*pi*s/dt);

% Normalised wavelet
psi=norm*psi0;
