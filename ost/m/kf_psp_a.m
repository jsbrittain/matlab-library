function kf_psp_a(sp,params,label);
%function kf_psp_a(sp,params,[label]);
%
% Kalman-Fourier plotting routine for auto-spectra ch.1
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       label       (opt) Additional label
%
%function kf_psp_a(sp,params,[label]);

% Construct label if absent
if (~exist('label'))
    label=kf_label(params);
end;

% Extract spectra
flen=size(sp,1);
faa=10*log10(reshape(double(sp(:,1,:)),flen,params.M));

% Determine x-axis and label
if (isfield(params,'xorder'))
    xorder=params.xorder;
    xwhat=params.xlabel;
else
    xorder=1:params.M;
    xwhat='TRIAL';
end;

% Plot spectra
if (false)
    surface(xorder,params.freqs,faa,'edgecolor','none');
    xlim([min(xorder) max(xorder)]);
    ylim([params.freqs(1) params.freqs(end)]);
    zlim([min(min(faa)) max(max(faa))]);
    view([140 50]);
else
    imagesc(xorder,params.freqs,faa);
end;
xlabel(xwhat);
ylabel('FREQ (Hz)');
zlabel('dB');
title(['fa ' label]);
set(gca,'ydir','normal');
