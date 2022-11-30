function [f,amp,ph,K] = mt_freqfit_disjoint( x, fs, flist, duration )
%function [f,amp,ph,K] = mt_freqfit_disjoint( x, fs, flist, duration )
%
% Parameters
%       x           Time-series
%       fs          Sample rate
%       flist       Interrogation frequecncy list
%       duration    Segment duration (msecs)
%
%function [f,amp,ph,K] = mt_freqfit_disjoint( x, fs, flist, duration )

% Check inputs
if (nargin<4)
    error(' All parametes must be specified.');
end;

% Epoch
dur = duration*fs/1000;
x = x(1:dur*floor(length(x)/dur));
x = reshape( x, dur, floor(length(x)/dur) );
K = size(x,2);

% Periodograms
f = zeros(1,K);
amp = zeros(length(flist),K);
ph = amp;
for k = (1:K)
    [f(k),amp(:,k),ph(:,k)] = mt_freqfit( x(:,k), fs, flist );
end;

% Determine peak of average amplitude spectrum
amp = mean( amp, 2 );
ph = angle( mean( exp(1i*ph), 2 ) );
[~,ix] = max( amp );
f = flist( ix );
