function y = fft_notch( x, rate, freq )
%function fft_notch( x, rate, freq )
%
% FFT-notch, zero +ve (and -ve) frequencies specified in the `freq' vector.
%
% Parameter `freq' should be specified as row tuples, i.e.
%       freq = [ 49 51 ];
%
%
%function fft_notch( x, rate, freq )

N = length(x);
freqs = (ceil(-N/2):ceil(N/2-1))*rate/N;
fx = fftshift(fft( x ));

for k = (1:size(freq,1))
    % +ve freqs
    fn = dsearchn( freqs', freq(k,:)' );
    fx(fn(1):fn(2)) = 0;
    % -ve freqs
    fn = dsearchn( freqs', -freq(k,:)' );
    fx(fn(2):fn(1)) = 0;
end;

% Reconstruct time-series
y = real( ifft(ifftshift( fx )) );
