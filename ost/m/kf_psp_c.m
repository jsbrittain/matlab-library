function kf_psp_c(sp,params,label);
%function kf_psp_c(sp,params,[label]);
%
% Kalman-Fourier plotting routine for cross-spectra ch.1,2
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       label       (opt) Additional label
%
%function kf_psp_c(sp,params,[label]);

% Construct label if absent
if (~exist('label'))
    label=kf_label(params);
end;

% Extract spectra
flen=size(sp,1);
fab=10*log10(abs(reshape(double(sp(:,3,:)),flen,params.M)));

% Determine x-axis and label
if (isfield(params,'xorder'))
    xorder=params.xorder;
    xwhat=params.xlabel;
else
    xorder=1:params.M;
    xwhat='Trials';
end;

% Plot spectra
surface(xorder,params.freqs,fab,'edgecolor','none');
xlim([min(xorder) max(xorder)]);
ylim([params.freqs(1) params.freqs(end)]);
zlim([min(min(fab)) max(max(fab))]);
view([140 50]);
xlabel(xwhat);
ylabel('Freq (Hz)');
zlabel('dB');
title(['fc ' label]);
set(gca,'ydir','normal');
