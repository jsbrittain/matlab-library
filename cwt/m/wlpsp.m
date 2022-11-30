function wlpsp(time,W,wlparam,opt);
%function wlpsp(time,W,wlparam,opt);
%
% Five figure summary plot depicting two autospectra, cross-spectra,
% coherence and phase.
%
% Input parameters
%       time        Time vector for abcissa
%       W           Wavelet spectral matrix (3-D)
%       wlparam     Wavelet parameters
%       opt         (opt) Options string
%                       c   -   Plot with colour bars
%
%function wlpsp(time,W,wlparam,opt);

% Check input parameters
if (nargin<3)
    error(' Not enough input parameters');
end;
if (nargin<4)
    opt='';
end;
if (length(time)~=size(W,2))
    error(' Time vector of inconsistent length with regard wavelet spectra');
end;

% Determine confidence limits
R95=[];
if (isfield(wlparam,'trigcount'))
    R95=(1-0.05^(1/(wlparam.trigcount)));
end;

% Draw summary plot
figure;
subplot(3,2,1);
wlpsp_a(time,W(:,:,1),wlparam,opt);
title('Autospectra channel 1');
subplot(3,2,2);
wlpsp_a(time,W(:,:,2),wlparam,opt);
title('Autospectra channel 2');
subplot(3,2,3);
wlpsp_a(time,abs(W(:,:,3)),wlparam,opt);
title('Crossspectra');
subplot(3,2,4);
wlpsp_coh(time,W(:,:,4),wlparam,R95,opt);
title('Coherence');
subplot(3,2,5);
wlpsp_ph(time,W(:,:,5),wlparam,opt);
title('Phase');
