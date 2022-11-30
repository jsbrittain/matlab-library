%addpath(genpath('c:/users/brittain/dropbox/matlab'))
addpath('../m');

%%

N = 1e6;
f = 5.0;
rate = 1e3;

switch ( 3 )
    case 1,     % Pure sine-wave (fixed phase)
        sigma = 3;
        dat = sin( 2*pi*f*(1:N)/rate ) + sigma*randn(1,N);
    case 2,     % Filtered white noise
        sigma   = 1;
        [bh,ah] = butter(3,(f-1)/(rate/2),'high');
        [bl,al] = butter(3,(f+1)/(rate/2),'low');
        dat = filtfilt( bh, ah, filtfilt( bl, al, randn(1,N) ) ) + randn(1,N);
    case 3,     % Phase-incrementing sine-wave
        sigmaph = 100;
        sigma   = 0;
        dat = cos( cumsum( 2*pi*f*(ones(1,N)+sigmaph*randn(1,N)) )/1e3 ) + sigma*randn(1,N);
end;

ft_freq = f;
ft_bandwidth = 4;
precision = 0.10;
rejectmode = 0;
analysis = freqtolerance( dat, ft_freq, ft_bandwidth, rate, precision, rejectmode );

figure;
subplot(4,2,[1 3]);
    plot( analysis.freqs, analysis.ifchangeMed ); hold('on');
    plot( analysis.freqs, analysis.ifchangeMed-analysis.ifchangeStd )
    plot( analysis.freqs, analysis.ifchangeMed+analysis.ifchangeStd )
    plot( analysis.freqs, analysis.ifchangeMed-analysis.ifchangeStd./sqrt(analysis.ifchangeCount), 'r' )
    plot( analysis.freqs, analysis.ifchangeMed+analysis.ifchangeStd./sqrt(analysis.ifchangeCount), 'r' )
    plot( xlim, [0 0], 'k' );
    xlim([0 10]);
    title(analysis.fit.freqtol);

    analysis.ifchangeMed( analysis.ifchangeCount < 20 ) = NaN;
    plot( analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0)), polyval(polyfit(analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0)),analysis.ifchangeMed(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0)),1),analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0))), 'g--', 'linewidth', 3 );
    plot( analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0)), polyval(polyfit(analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0)),analysis.ifchangeMed(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0)),1),analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0))), '--', 'color', [1 1 1]*0.8, 'linewidth', 3 );
    plot( analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0)), polyval(polyfit(analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0)),analysis.ifchangeMed(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0)),1),analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0))), 'g--', 'linewidth', 3 );
    plot( analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0)), polyval(polyfit(analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0)),analysis.ifchangeMed(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed<0)),1),analysis.freqs(~isnan(analysis.ifchangeMed) & (analysis.ifchangeMed>0))), '--', 'color', [1 1 1]*0.8, 'linewidth', 3 );

subplot(4,2,5);
    freqs = (1:precision:10);
    phconsist = zeros(1,length(freqs));
    for k = (1:length(freqs))
        phconsist(k) = abs(mean( dat.*exp(1i*2*pi*freqs(k)*(1:length(dat))/rate) ));
    end;
    plot( freqs, phconsist );
    xlim([0 10]);
    
if ( exist('mt_sp0','file') )
    subplot(4,2,7);
        [sp11,sp22,sp12,params]=mt_sp0(dat',dat',2000,rate,'W0.5 f10');
        plot( params.freqs, sp11 );
        xlim([0 10]);
end

%% SNR comparison

N = 1e5;
f = 5.0;
rate = 1e2;

% Filtered white noise
[bh,ah] = butter(3,(f-1)/(rate/2),'high');
[bl,al] = butter(3,(f+1)/(rate/2),'low');
%signal = zscore( filtfilt( bh, ah, filtfilt( bl, al, randn(1,N) ) ) );
signal = zscore( sin(2*pi*f*(1:N)/rate) );
var_signal = var( signal );

snrlist = (-100:2:100);     % dB
iters   = 1;
trifit  = zeros(length(snrlist),iters);

ft_freq      = f;
ft_bandwidth = 4;
precision    = 0.10;
rejectmode   = 0;

for k = (1:length(snrlist))
    [ k length(snrlist) ]
    for n = (1:iters)
        
        dat = zscore( signal + sqrt( var_signal*(10^(-snrlist(k)/10)) )*randn(1,N) );
        analysis = freqtolerance( dat, ft_freq, ft_bandwidth, rate, precision, rejectmode );
        
        trifit(k,n) = analysis.fit.piecefit(4)-analysis.fit.piecefit(3);
    end;
end;

%%

figure;
plot( snrlist, median(trifit,2) ); hold('on');
plot( snrlist, quantile(trifit',0.25) );
plot( snrlist, quantile(trifit',0.75) );
