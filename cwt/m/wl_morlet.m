function psi=wl_morlet(nu,w0)
%function psi=wl_morlet(nu,w0)
%
%  Wavelet: Morlet
%  Time-domain
%
% Common parameters
%   nu      Time
%
% Additional parameters
%   w0      Central frequency of the mother wavelet
%
% NOTE: Normalisation not included like in the frequency domain definition.
%       To normalise multiply by "sqrt(dt/scale)" for each scale.

psi=(pi^(-1/4)).*exp(1i*w0*nu).*exp(-(nu.^2)/2);
