function psi=wl_morse_composite_energy(nu,b,g,wk)
%function psi=wl_morse_composite_energy(nu,beta,gamma,wk);
%
%  Wavelet: Composite Morse energy spectrum (Time-domain)
%
% Common parameters
%   nu      Time (secs)
%
% Additional parameters
%   beta    Beta parameter
%   gamma   Gamma parameter
%   K       Wavelet order
%
% NOTE: OrthoNORMAL wavelets
%

K = length(wk);
psi = wk(1)*abs( wl_morse(nu,b,g,0) ).^2;
for k = (1:(K-1))
    psi = psi + wk(k+1)*abs( wl_morse(nu,b,g,k) ).^2;
end;
