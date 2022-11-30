function [ dat0, F1 ] = mt_zeropad_freq( dat, N )

% Zero pad Fourier transform (zero mean and energy normalise segment first)
F1 = fft(dat)/sqrt(length(dat));        % Fourier transform
zcount = N-length(F1);                  % Zero pad count

% Zero pad in freq domain
F1=[ F1(1:ceil(length(F1)/2)) zeros(1,zcount) F1(ceil(length(F1)/2+1):end) ];

% Apply power correction factor for zero padding
F1=F1*sqrt(N/length(dat));

% Recover time-domain signal
dat0 = real(ifft(F1)*sqrt(length(F1)));
