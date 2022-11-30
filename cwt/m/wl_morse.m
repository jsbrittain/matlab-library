function psi=wl_morse(nu,b,g,k)
%function psi=wl_morse(nu,beta,gamma,k);
%
%  Wavelet: Morse
%  Time-domain
%
% Common parameters
%   nu      Time (secs)
%
% Additional parameters
%   beta    Beta parameter
%   gamma   Gamma parameter
%   k       Wavelet order
%
% NOTE: OrthoNORMAL wavelets
%

% Determine data parameters
N=length(nu);
halflength=floor(N/2);
if (length(nu)>1)
    dt=nu(2)-nu(1);
else
    dt = 1;
end;

% Construct wavelet in the frequency domain
omega=[(2*pi*(0:N-halflength-1)/(N*dt)) (-(2*pi*((0:halflength-1)))/(N*dt))]';
WF=wlf_morse(omega,1,dt,b,g,k);

% Inverse FFT to obtain time-domain version
psi=fftshift(ifft(WF));
