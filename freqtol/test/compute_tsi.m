function analysis = compute_tsi( dat, fc, fbandwidth, rate, rejectmode )
%function analysis = compute_tsi( dat, fc, fbandwidth, rate, rejectmode )
%
% TSI is given as iqr(analysis.deltaIF)
%
%function analysis = compute_tsi( dat, fc, fbandwidth, rate, rejectmode )

if (~isvector(dat))
    error(' Input must be univariate vector time-series.');
end
dat=dat(:);

% Detrend
[bh,ah] = butter(3,(2.5)/(rate/2),'high');
dat = filtfilt(bh,ah,filtfilt(bh,ah,dat));

% Filter
if (~isempty(fc))
    [bh,ah] = butter(3,(fc-fbandwidth/2)/(rate/2),'high');
    [bl,al] = butter(3,(fc+fbandwidth/2)/(rate/2),'low');
    datf = filtfilt(bh,ah,filtfilt(bl,al,dat));
else
    datf = dat;
end
mag = smooth( abs(hilbert(datf)), 1.0*rate );
mag = mag(:);
noise = abs( smooth(abs(hilbert(dat)),1.0*rate) - mag );
noise = noise(:);
snr = 20*log10( mag./noise );
snr = snr(:);
time = (1:length(dat))/rate;
time = time(:);

% Threshold crossing
crossingtimes = find(diff(datf>0)>0);
crossingtimes = crossingtimes(:);
instperiod = diff(time(crossingtimes));
instfreq = 1./instperiod;
instmag = mag(crossingtimes);
instsnr = snr(crossingtimes);
deltaIF = diff(instfreq);

% Reject data
switch ( rejectmode )
    case 0, keep = true(size(instmag));
    case 1, keep = instmag > median(instmag);
    case 2, keep = instmag < median(instmag);
    case 3, keep = instmag > quantile(instmag,0.05);
    otherwise
        warning('Unknown reject mode specified!');
end

%
analysis.instmag = instmag(keep);
analysis.instsnr = instsnr(keep);
analysis.instfreq = instfreq(keep(1:end-1));
analysis.deltaIF = deltaIF(keep(1:end-2));
