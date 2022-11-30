function kftf_psp(sp,params,freq)
%function kftf_psp(sp,params,freq)
%
% Hybrid Kalman-Fourier analysis Time-Frequency
% Coherence plotting routine
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       freq        Plot frequency
%
%function kftf_psp(sp,params,freq)

% Check input parameters
if (nargin~=3)
    error(' Incorrect number of input parameters');
end;

% Determine frequency location
fpos=dsearchn(params.freqs',freq);

% Determine 95% confidence limit
Lprime=params.Lprime(fpos,:); Lprime=transpose(Lprime(ones(size(sp,4),1),:));
warning off MATLAB:divideByZero
R95=(1-0.05.^(1./(Lprime-1)));
warning on MATLAB:divideByZero

% Surface plot
trials=size(sp,3); segs=size(sp,4);
%surface([0:(segs-1)]/segs,[1:trials],squeeze(sp(fpos,4,:,:))); hold on;
surface([0:(segs-1)]/segs,[1:trials],double(squeeze(sp(fpos,4,:,:)))); hold on;
surface([0:(segs-1)]/segs,[1:trials],R95,'edgecolor','none','facecolor','k','facealpha',0.75);
xlabel('Offset'); ylabel('Trials'); axis('tight'); view([140 50]);
title(['Absolute-relative time coherence (' int2str(params.freqs(fpos)) 'Hz) ' params.what]);
