function analysis = freqtolerance_timeseries( dat, fc, fbandwidth, rate, precision, rejectmode, varargin )
%function analysis = freqtolerance_timeseries( dat, fc, fbandwidth, rate, precision, rejectmode, <param-value> pairs )
%
%
%
%function analysis = freqtolerance_timeseries( dat, fc, fbandwidth, rate, precision, rejectmode, <param-value> pairs )

if (~exist('precision','var'))
    precision = [];
end
if (isempty(precision))
    precision = 10;
end

if (~isvector(dat))
    error(' Input must be univariate vector time-series.');
end
if (size(dat,2)>size(dat,1))
    dat = dat.';
end

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

% Subset select data and pass on to analysis routine
instmag = instmag(keep);
instsnr = instsnr(keep);
instfreq = instfreq(keep(1:end-1));
deltaIF = deltaIF(keep(1:end-2));

% Add SNR to parameter-value pairs
varargin{end+1} = 'instsnr';
varargin{end+1} = instsnr;

% Analyse crossing times
analysis = freqtolerance_instfreq( instfreq, instmag, precision, deltaIF, varargin{:} );

% Additional parameters as needed
analysis.crossingtimes      = crossingtimes(keep);
analysis.instperiod         = instperiod(keep(1:end-1));
