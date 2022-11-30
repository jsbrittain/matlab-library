function [f,amp,ph] = mt_freqfit( x, fs, flist )
%function [f,amp,ph] = mt_freqfit( x, fs, flist )
%
% Parameters
%       x       Time-series
%       fs      Sample rate
%       flist   Interrogation frequecncy list
%
%function [f,amp,ph] = mt_freqfit( x, fs, flist )

% Check inputs
if (nargin<3)
    error(' All parametes must be specified.');
end;
if (~isvector(x))
    error(' Input x must be a vector.');
end;
if (size(x,1)>1)
    x = x.';
end;

% Determine frequency
x = x - mean(x);
for fn = (1:length(flist))
    yc = sum( x.*cos(2*pi*flist(fn).*(0:(length(x)-1))/fs) );
    ys = sum( x.*sin(2*pi*flist(fn).*(0:(length(x)-1))/fs) );
    amp(fn) = abs( yc + 1i*ys );
    ph(fn) = angle( yc + 1i*ys );
end;
[~,ix] = max( amp );
f = flist(ix);
