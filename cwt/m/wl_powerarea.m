function P = wl_powerarea(time,freqs,sp,tspan,fspan)
%
% Extract power of time-frequency area designated by tspan, fspan
%
% Spectrogram expected in ( freq x time ) format
%

% Check spectrogram structure against time/freq vectors
sp = squeeze(sp);
if ((length(time)~=size(sp,2)) || (length(freqs)~=size(sp,1)))
    error(' Spectrum does not match time-frequency vectors.');
end;

% Check input parameters (default to full range)
if (isempty(tspan))
    tspan = time([1 end]);
end;
if (isempty(fspan))
    fspan = freqs([1 end]);
end;

% Extract time/frequency region indicies
trange = dsearchn( time', tspan' );
frange = sort( dsearchn( freqs', fspan' ) );

% Normalisation factors
dt = time(2) - time(1);
df = abs( diff(freqs(frange(1):frange(2))) );

% Take mid-points of time/freq grid
sp = sp( frange(1):frange(2), trange(1):trange(2) );
sp = sp(1:(end-1),:) + diff(sp,[],1)/2;
sp = sp(:,1:(end-1)) + diff(sp,[],2)/2;

% Determine power
P = mean( mean( sp/dt, 2)./df', 1 );
