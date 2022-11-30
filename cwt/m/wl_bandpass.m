function P = wl_bandpass(freqs,sp,fspan)
%
% Extract power of time-frequency area designated by tspan, fspan
%
% Spectrogram expected in ( freq x time ) format
%

% Check spectrogram structure against time/freq vectors
sp = squeeze(sp);
if (length(freqs)~=size(sp,1))
    error(' Spectrum does not match time-frequency vectors.');
end;

% Check input parameters (default to full range)
if (isempty(fspan))
    fspan = freqs([1 end]);
end;

% Extract time/frequency region indicies
frange = sort( dsearchn( freqs', fspan' ) );

% Normalisation factors
df = abs( diff(freqs(frange(1):frange(2))) );

% Take mid-points of time/freq grid
sp = sp( frange(1):frange(2), : );
% sp = sp(1:(end-1),:) + diff(sp,[],1)/2;
% sp = sp(:,1:(end-1)) + diff(sp,[],2)/2;

% Determine power
%P = mean( sp./repmat(df',1,size(sp,2)), 1 );
P = mean( sp, 1 );
