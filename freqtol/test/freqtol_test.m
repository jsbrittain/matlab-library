N = 1e6;
f = 5.0;
rate = 1e3;

switch ( 1 )
    case 1      % Pure sine-wave (fixed phase)
        sigma = 3;
        dat = sin( 2*pi*f*(1:N)/rate ) + sigma*randn(1,N);
    case 2      % Filtered white noise
        sigma   = 1;
        [bh,ah] = butter(3,(f-1)/(rate/2),'high');
        [bl,al] = butter(3,(f+1)/(rate/2),'low');
        dat = filtfilt( bh, ah, filtfilt( bl, al, randn(1,N) ) ) + randn(1,N);
    case 3      % Phase-incrementing sine-wave
        sigmaph = 100;
        sigma   = 0;
        dat = cos( cumsum( 2*pi*f*(ones(1,N)+sigmaph*randn(1,N)) )/1e3 ) + sigma*randn(1,N);
end

ft_freq = f;
ft_bandwidth = 4;
precision = 0.10;
rejectmode = 0;
analysis = freqtolerance( dat, ft_freq, ft_bandwidth, rate, precision, rejectmode );

iqr(analysis.deltaIF)
