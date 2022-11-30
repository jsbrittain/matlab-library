function kfslice_psp_ph(sp,params,freq,label);
%function kfslice_psp_ph(sp,params,freq,[label]);
%
% Kalman-Fourier plotting routine for phase
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       label       (opt) Additional label
%
%function kfslice_psp_ph(sp,params,freq,[label]);

% Construct label if absent
if (~exist('label'))
    label=kf_label(params);
end;

% Extract phase
fpos=dsearchn(params.freqs',freq);
phase=angle(squeeze(sp(fpos,3,:)));

% Determine x-axis and label
if (isfield(params,'xorder'))
    xorder=params.xorder;
    xwhat=params.xlabel;
else
    xorder=1:params.M;
    xwhat='Trials';
end;

% Plot phase
ax(1)=newplot;
if (isfield(params,'jk'))
    h=plot(xorder,phase-1.96*squeeze(sqrt(params.jk(fpos,5,:))), ...
           xorder,phase+1.96*squeeze(sqrt(params.jk(fpos,5,:)))  );
    set(h,'color',0.6*[1 1 1]);         % Solid gray line
    hold on;
end;
plot(xorder,phase,'k');
hold off;
xlim([min(xorder) max(xorder)]);
ylim([-pi pi]);
xlabel(xwhat);
title(['ph ' int2str(round(params.freqs(fpos))) 'Hz ' label]);

% Add second axis
if (isfield(params,'y2'))
    y2colour=[10 118 74]/255;
    set(ax(1),'box','off'); ylabel('ph');
    ax(2)=axes('position',get(ax(1),'position'));
    plot(xorder,params.y2,'color',y2colour);
    set(ax(2),'YAxisLocation','right','color','none','box','off','ycolor',y2colour);
    xlim([min(xorder) max(xorder)]);
    set(ax(1),'activepositionproperty','position');
    set(ax(2),'activepositionproperty','position');
    if (isfield(params,'y2label'))
        ylabel(params.y2label);
    end;
    subplot(ax(1));
end;
