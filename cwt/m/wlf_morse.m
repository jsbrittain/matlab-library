function [psi]=wlf_morse(omega,s,dt,b,g,k,dk)
%function [psi]=wlf_morse(omega,s,dt,beta,gamma,k,[dk])
%
% Wavelet: Morse
% Fourier transform
%
% Common parameters
%   omega   Frequencies
%   s       Scale
%   dt      Time resolution
%
% Additional parameters
%   beta    beta >= 1
%   gamma   gamma > (beta-1)/2
%   k       Morse order (i.e. Zeroth is Klauder)
%
% NB: (Torrence & Compo 1998) normalisation used so that sum(abs(psi)^2)=N.

% Check input parameters
if ((nargin<6) || (nargin>7))
    error(' Incorrect number of wavelet parameters');
end;

% Calculate wavelet response in the frequency domain
sw=s*omega;
r=(2*b+1)/g;
c=r-1;
A=sqrt(pi*g*(2^r)*gamma(k+1)/gamma(k+r));
% Generate Laguerre function
L=0;
x=2*(sw.^g);
for m=0:k
    L=L + ((-1)^m)*(gamma(k+c+1)/(gamma(c+m+1)*gamma(k-m+1)))*((x.^m)/factorial(m));
end;
% Specify wavelet in the frequency domain
psi0=sqrt(2)*A*(sw.^b).*exp(-(sw.^g)).*L;   % Analytic wavelet
psi0(sw<=0)=0;                              %  |

% Normalised wavelet (unit energy per scale, subsequently sum(abs(psi)^2)=N)
norm=sqrt(s/dt);                            % As (Torrence & Compo) but using Hz
psi=norm*psi0;
