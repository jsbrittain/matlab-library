function sig=wlCohSig(nco,ncy);
%function sig=wlCohSig(nco,ncy);
%
% Significance threshold value for single trial wavelet coherence
% (P=0.005) for various nco, ncy values.
%
% Input parameters
%   nco     No. of cycles of the wavelet
%   ncy     No. of cycles contained in the integration window
%           Window: [t-delta/2 : t+delta/2]
%
% Supported range
%   nco     4-8
%   ncy     4-10
%
% To be used in conjunction with the single-trial coherence estimates
% as described in Lachaux et al. (2002).  Values taken from the same
% articles, Table 1, p.170.
%
% Ref:  Lachaux J-P, et al. (2002) Estimating the time-course of coherence
%       between single-trial brain signals: an introduction to wavelet
%       coherence. Neurophysiol. Clin. 32:157-174
%


sig=[ 0.88 0.80 0.75 0.73 0.64 0.60 0.60;
      0.90 0.88 0.81 0.76 0.75 0.72 0.67;
      0.94 0.91 0.88 0.82 0.77 0.76 0.71;
      0.95 0.93 0.90 0.87 0.84 0.79 0.78;
      0.97 0.96 0.93 0.90 0.88 0.86 0.84 ];

if ((nco>8) | (nco<4)) | ((ncy>10) | (ncy<4))
    disp('Significance level not supported: Returning zero');
    sig=0;
    return;
end;
sig=sig(nco-3,ncy-3);
