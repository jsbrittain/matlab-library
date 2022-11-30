function kf_psp_ph(sp,params,label);
%function kf_psp_ph(sp,params,[label]);
%
% Kalman-Fourier plotting routine for phase
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       label       (opt) Additional label
%
%function kf_psp_ph(sp,params,[label]);

% Construct label if absent
if (~exist('label'))
    label=kf_label(params);
end;

% Extract phase
flen=size(sp,1);
phase=angle(reshape(double(sp(:,3,:)),flen,params.M));

% Determine x-axis and label
if (isfield(params,'xorder'))
    xorder=params.xorder;
    xwhat=params.xlabel;
else
    xorder=1:params.M;
    xwhat='TRIAL';
end;

% Plot phase
%surface(xorder,params.freqs,phase,'edgecolor','none'); view([140 50]);
imagesc(xorder,params.freqs,phase);
xlim([min(xorder) max(xorder)]);
ylim([params.freqs(1) params.freqs(end)]);
zlim([-pi pi]);
xlabel(xwhat);
ylabel('FREQ (Hz)');
title(['ph ' label]);
set(gca,'ydir','normal');
