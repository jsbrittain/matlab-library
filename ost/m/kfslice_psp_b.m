function kfslice_psp_b(sp,params,freq,label)
%function kfslice_psp_b(sp,params,freq,[label])
%
% Kalman-Fourier plotting routine for a frequency slice through
% auto-spectra channel 1
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       freq        Slice frequency
%       label       (opt) Additional label
%
%function kfslice_psp_b(sp,params,freq,[label])

% Construct label if absent
if (~exist('label'))
    label=kf_label(params);
end;

% Extract spectra
fpos=dsearchn(params.freqs',freq);
ch=2;
faa=10*log10(squeeze(sp(fpos,ch,:)));

% Determine 95% conf limits
cl95_lower=[]; jkcl95_lower=[];
switch (params.method)
    case 0          % null-smoothing - use periodogram noise (provided)
        cl95=(10*log10(exp(1)))*(1.96*squeeze(sqrt(double(params.Psp(fpos,ch,:)))));
        cl95_lower=faa-cl95; cl95_upper=faa+cl95;
    case 2          % kflog - use returned variance
        cl95=(10*log10(exp(1)))*(1.96*squeeze(sqrt(double(params.Psp(fpos,ch,:)))));
        cl95_lower=faa-cl95; cl95_upper=faa+cl95;
        if (isfield(params,'jk'))
            jkcl95=(10*log10(exp(1)))*(1.96*squeeze(sqrt(params.jk(fpos,ch,:))));
            jkcl95_lower=faa-jkcl95; jkcl95_upper=faa+jkcl95;
        end;
    case 5          % kfroot - use returned variance
        cl95_lower=real(10*log10(squeeze(((sp(fpos,ch,:).^(1/params.root))-1.96*sqrt(params.Psp(fpos,ch,:)))).^params.root));
        cl95_upper=real(10*log10(squeeze(((sp(fpos,ch,:).^(1/params.root))+1.96*sqrt(params.Psp(fpos,ch,:)))).^params.root));
        if (isfield(params,'jk'))
            jkcl95_lower=real(10*log10(squeeze(((sp(fpos,ch,:).^(1/params.root))-1.96*sqrt(params.jk(fpos,ch,:))).^params.root)));
            jkcl95_upper=real(10*log10(squeeze(((sp(fpos,ch,:).^(1/params.root))+1.96*sqrt(params.jk(fpos,ch,:))).^params.root)));
        end;
end;

% Determine x-axis and label
if (isfield(params,'xorder'))
    xorder=params.xorder;
    xwhat=params.xlabel;
else
    xorder=1:params.M;
    xwhat='Trials';
end;

% Plot spectra
ax(1)=newplot; hold('on');
if (~isempty(cl95_lower))
    if (true)
        % Red lines
        h=plot(xorder,cl95_lower,xorder,cl95_upper);
        set(h,'color',[0.8 0.2 0.2]);      % Deep red line
    else
        % Gray fill
        h=fill([xorder fliplr(xorder) xorder(1)],[cl95_lower; flipud(cl95_upper); cl95_lower(1)],'k');
        set(h,'edgecolor','none');
        set(h,'facecolor',0.75*[1 1 1]);
    end;
end;
if (~isempty(jkcl95_lower))
    h=plot(xorder,jkcl95_lower,xorder,jkcl95_upper);
    set(h,'color',0.6*[1 1 1]);         % Solid gray line
end;
h=plot(xorder,faa,'k');               % Plot spectra last to overlay conf limits
set(h,'linewidth',1.5);
hold off;
xlim([min(xorder) max(xorder)]);
xlabel(xwhat);
ylabel('dB');
title(['fa ' int2str(round(params.freqs(fpos))) 'Hz ' label]);

% Add second axis
if (isfield(params,'y2'))
    y2colour=[10 118 74]/255;
    set(ax(1),'box','off'); ylabel('dB');
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
