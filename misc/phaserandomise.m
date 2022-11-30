function y = phaserandomise( x )
%function y = phaserandomise( x )
%
% Phase randomise time-series
%
% Works by taking the Fourier transform of the whole data segment,
% randomising all phase arguments, forcing conjugate symmetry in the
% Cartesian representation then inverse Fourier transforming.
%
% Input parameter
%       x       Time-series (vector input)
%
% Output parameter
%       y       Phase-randomised time-series
%
%function y = phaserandomise( x )

% Check input parameters
if (~isvector( x ))
    error(' Vector inputs only!');
end;

% Data parameters
N = length(x);

% FFT
fx=fft(x);

% Separate amplitude and phase ( Cartesian -> Polar form )
amp = abs(fx);
ph = angle(fx);                     % NB: Not used

% Phase randomise
ph = 2*pi*rand(size(ph))-pi;        % [ -pi, pi ]

% Polar -> Cartesian
fx = amp.*exp(1i*ph);

% Enforce conjugate symmetry ( for real-valued time-series )
fx(2:ceil(N/2)) = conj( rot90(fx(floor(N/2+2):N),2) );       % Accounts for odd/even sample sizes

% Inverse FFT ( real operator discards computational discrepancies )
y = real(ifft(fx));
