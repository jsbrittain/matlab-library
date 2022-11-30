function df_psp_coh(coh11,coh22,coh12,params,fmax,smoothing_kernel)
%function df_psp_coh(coh11,coh22,coh12,params,maxf,smoothing_kernel)
%
% Summary plot of dual-frequency coherence between two processes
%
% Input parameters
%       coh11       Dual-frequency auto-coherence ch.1
%       coh22       Dual-frequency auto-coherence ch.2
%       coh12       Dual-frequency cross-coherence ch.1,2
%       params      Parameters structure
%
%function df_psp_coh(coh11,coh22,coh12,params,maxf,smoothing_kernel)

if (~exist('fmax','var'))
    fmax=[];
end;
if (~exist('smoothing_kernel','var'))
    smoothing_kernel = [];
end;

if (isempty(fmax))
    if isfield(params,'fmax')
        fmax=params.fmax;
    else
        fmax = params.freqs(end);
    end;
end;
if (isempty(smoothing_kernel))
    smoothing_kernel = 1;
end;
smoothing_kernel = smoothing_kernel/sum(smoothing_kernel(:));
if (~isfield(params,'amplitudecorrelation'))
    amplitudecorrelation = false;
else
    amplitudecorrelation = params.amplitudecorrelation;
end;

% Determine 95% significance level
freqs=params.freqs;
Kprime=params.L;
R95=1-0.05.^(1./(Kprime-1));

% Plot DF coherence
figure;
cutoff=find(abs(freqs-fmax)==min(abs(freqs-fmax)));
freqs = freqs(1:cutoff);
coh11 = coh11(1:cutoff,1:cutoff);
coh22 = coh22(1:cutoff,1:cutoff);
coh12 = coh12(1:cutoff,1:cutoff);
plotconflimit=true;

subplot(2,2,1);
imagesc(freqs,freqs,conv2(coh11,smoothing_kernel,'same'));
colorbar;
if ( plotconflimit && ( isscalar(smoothing_kernel) ))
    hold on;
    [~,h]=contour(freqs,freqs,(coh11>=R95),1,'k');
    set(h,'linewidth',4);
end;
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
set(gca,'ydir','normal');
set(gca,'PlotBoxAspectRatio',[1 1 1]);
if (amplitudecorrelation)
    title('DF correlation ch.1');
else
    title('DF coherence ch.1');
end;

subplot(2,2,2);
imagesc(freqs(1:cutoff),freqs(1:cutoff),conv2(coh22(1:cutoff,1:cutoff),smoothing_kernel,'same'));
colorbar;
if ( plotconflimit && ( isscalar( smoothing_kernel ) ))
    hold on;
    [~,h]=contour(freqs(1:cutoff),freqs(1:cutoff),(coh22(1:cutoff,1:cutoff)>=R95),1,'k');
    set(h,'linewidth',4);
end;
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
set(gca,'ydir','normal');
set(gca,'PlotBoxAspectRatio',[1 1 1]);
if (amplitudecorrelation)
    title('DF correlation ch.2');
else
    title('DF coherence ch.2');
end;

subplot(2,2,3);
imagesc(freqs(1:cutoff),freqs(1:cutoff),conv2(coh12(1:cutoff,1:cutoff),smoothing_kernel,'same'));
colorbar;
if ( plotconflimit && ( isscalar(smoothing_kernel) ))
    hold on;
    [~,h]=contour(freqs(1:cutoff),freqs(1:cutoff),(coh12(1:cutoff,1:cutoff)>=R95),1,'k');
    set(h,'linewidth',4);
end;
hold('on'); plot( xlim, ylim, 'w', 'linewidth', 2 );
ylabel('Freq ch.2 (Hz)'); xlabel('Freq ch.1 (Hz)');
set(gca,'ydir','normal');
set(gca,'PlotBoxAspectRatio',[1 1 1]);
if (amplitudecorrelation)
    title('DF correlation ch.1,2');
else
    title('DF coherence ch.1,2');
end;

subplot(2,2,4);
if (amplitudecorrelation)
    plot(freqs,diag(coh12),'k',freqs,zeros(size(freqs)),'k-');
    title('Stationary correlation ch.1,2');
else
    plot(freqs,abs(diag(coh12)),'k',freqs,R95*ones(size(freqs)),'k--');
    title('Stationary coherence ch.1,2');
end;
xlim(freqs([1 end])); xlabel('Freq (Hz)');

% Set clim
cmax = max([max(coh11-diag(diag(coh11))) max(coh22-diag(diag(coh22)))]);
ymax = max(coh12(:));
if (amplitudecorrelation)
    subplot(2,2,1); set(gca,'clim',max(get(gca,'clim'))*[-1 1]);
    subplot(2,2,2); set(gca,'clim',max(get(gca,'clim'))*[-1 1]);
    subplot(2,2,4); ylim([-ymax ymax]);
else
    subplot(2,2,1); set(gca,'clim',[0 cmax]);
    subplot(2,2,2); set(gca,'clim',[0 cmax]);
    subplot(2,2,4); ylim([0 ymax]);
end;
%subplot(2,2,3); set(gca,'clim',[0 ymax]);
