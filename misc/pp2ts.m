function ts = pp2ts( pp, fs, N )
%function pp2ts( pp, fs, N )
%
% Function to convert a point-process (pp) vector into a time-series (sr)
%
% Inputs
%   pp      Point-process data (vector of times [secs])
%   fs      Sample rate of new time-series vector
%   N       Length (in samples) of new time-series vector
%            (default [] = maximal time in point-process)
%
%function pp2ts( pp, fs, N )

pp = round( pp*fs );
pp(pp<1) = [];
if ( isempty(N) )
    N = max(pp);
end
if ( max(pp) > N )
    warning('Point-process indexing is beyond specified time-series length... extending vector.');
    N = max(pp);
end

ts = zeros( 1, N );
ts(pp) = 1;
