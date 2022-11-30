function kf_psp_coh(sp,params,label);
%function kf_psp_coh(sp,params,[label]);
%
% Kalman-Fourier plotting routine for coherence
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       label       (opt) Additional label
%
%function kf_psp_coh(sp,params,[label]);

% Construct label if absent
if (~exist('label'))
    label=kf_label(params);
end;

% Extract coherence
coh=squeeze(double(sp(:,4,:)));

% Determine 95% confidence limit
if (~params.jackknife)
    % Analytic limits
    warning off MATLAB:divideByZero
    R95=(1-0.05.^(1./(double(params.Lprime)-1)));
    warning on MATLAB:divideByZero
else
    % Jackknife limits
    R95=tanh(squeeze(2*1.96*sqrt(params.jk(:,4,:)))).^2;
end;

% Determine x-axis and label
if (isfield(params,'xorder'))
    xorder=params.xorder;
    xwhat=params.xlabel;
else
    xorder=1:params.M;
    xwhat='TRIAL';
end;

% Plot coherence
if (false)
    surface(xorder,params.freqs,coh,'edgecolor','none');
    hold on;
    surface(xorder,params.freqs,R95,'edgecolor','none','facecolor','k','facealpha',0.75);
    zlim([0 max(max(coh))]);
    view([140 50]);
else
    imagesc(xorder([1 end]),params.freqs([1 end]),coh);
end;
xlim([min(xorder) max(xorder)]);
ylim([params.freqs(1) params.freqs(end)]);
xlabel(xwhat);
ylabel('FREQ (Hz)');
title(['coh ' label]);
set(gca,'ydir','normal');
