
clear all

% Parameters
N = 1e4-1;          % Trial length (samples)
M = 45;             % Trial count
rate = 100;         % Sampling rate (used to construct frequency vector)
taper = 2;          % Taper time-series (0=none, 1=Gaussian, 2=hamming)
smoothing = 1;      % Smoothing over frequencies (reduces variance)

% Generate data
switch ( 2 )
    case 1,     % Common sinusoid with noise
        x = 1 + randn(N,M) + repmat(sin(2*pi*20*(0:N-1)/100)',1,M);
        y = 1 + randn(N,M) + repmat(sin(2*pi*20*(0:N-1)/100)',1,M);
    case 2,     % Common frequency band with noise
        [bh,ah] = butter(3,10/(rate/2),'high');
        [bl,al] = butter(3,15/(rate/2),'low');
        for k = (1:M)
            x(:,k) = filtfilt(bh,ah,filtfilt(bl,al,randn(N,1)));
            y(:,k) = x(:,k) + randn(N,1);
            x(:,k) = x(:,k) + randn(N,1);
        end;
end;

% Select taper
switch ( taper )
    case 0,     % None
        w = ones(N,1);
    case 1,     % Gaussian
        w = pdf('normal',1:N,N/2,N/4)';
    case 2,     % Hamming
        w = hamming( N );
end;
w = w/sum(w);   % Normalise taper weight

% Generate single-trial periodograms
for k = (1:M)
    % Fourier transform
    fx = fft(w.*x(:,k));                        % Apply taper to time-series THEN zero-pad (not implemented) THEN transform
    fy = fft(w.*y(:,k));                        %  |
    % Periodograms
    Pxx(:,k) = fx.*conj(fx);%abs(fx).^2;        % Real-valued    (complex parts of Fourier transform cancel)
    Pyy(:,k) = fy.*conj(fy);%abs(fy).^2;        % Real-valued    (complex parts of Fourier transform cancel)
    Pxy(:,k) = fx.*conj(fy);                    % Complex-valued (constructed from 2 different time-series)
end;

% Generate Spectra (trial-average periodograms)
Sxx = mean( fftshift(Pxx), 2 );                 % fftshift rotates form +ve freqs followed by -ve freqs
Syy = mean( fftshift(Pyy), 2 );                 %  to -ve freqs followed by +ve freqs (i.e. the way we
Sxy = mean( fftshift(Pxy), 2 );                 %  are used to viewing it).

% Frequency smoothing
Sxx = smooth( Sxx, smoothing );
Syy = smooth( Syy, smoothing );
Sxy = smooth( Sxy, smoothing );

% Frequency vector
freqs = (-floor(N/2):ceil(N/2-1))/N*rate;

% Coherence 95% confidence limit
if (smoothing == 1)
    R95 = (1-0.05^(1/(M-1)));
end;

% Plot results
figure;
subplot(321);
    plot( freqs, 10*log10(Sxx), 'k');
    ylabel('dB');
    title('AUTOSPECTRUM x');
subplot(322);
    plot( freqs, 10*log10(Syy), 'k');
    title('AUTOSPECTRUM y');
    ylabel('dB');
subplot(323);
    plot( freqs, 10*log10(abs(Sxy)), 'k');
    title('MAGNITUDE CROSS-SPECTRUM');
    ylabel('dB');
subplot(324);
    plot( freqs, angle(Sxy), 'k');
    title('PHASE CROSS-SPECTRUM');
    ylabel('\phi');
subplot(325);
    plot( freqs, abs(Sxy).^2./(Sxx.*Syy), 'k');
    if (exist('R95'))
        hold('on'); plot(xlim,R95*[1 1],'k--');
    end;
    title('COHERENCE');
    ylabel('|R_{xy}|^2');
    
% Coherence bias
1/M
mean( abs(mean(Sxy,2)).^2./(mean(Sxx,2).*mean(Syy,2)) );
