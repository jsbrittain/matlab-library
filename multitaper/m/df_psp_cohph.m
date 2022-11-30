function df_psp_cohph(coh11,coh22,coh12,ph11,ph22,ph12,params)
%function df_psp_cohph(coh11,coh22,coh12,ph11,ph22,ph12,params)
%
% Summary plot of dual-frequency coherence between two processes
%
% Input parameters
%       coh11       Dual-frequency auto-coherence ch.1
%       coh22       Dual-frequency auto-coherence ch.2
%       coh12       Dual-frequency cross-coherence ch.1,2
%       sp11
%       sp22
%       sp12
%       params      Parameters structure
%
%function df_psp_cohph(coh11,coh22,coh12,ph11,ph22,ph12,params)

% Determine 95% significance level
freqs=params.freqs;
maxf=params.fmax;
Kprime=params.L;
R95=1-0.05^(1/(Kprime-1));

% Plot DF coherence
figure;
cutoff=find(abs(freqs-maxf)==min(abs(freqs-maxf)));
plotconflimit=logical(1);
subplot(4,2,1);
imagesc(freqs(1:cutoff),freqs(1:cutoff),flipud(coh11(1:cutoff,1:cutoff)));
colorbar;
if (plotconflimit)
    hold on;
    contour(freqs(1:cutoff),freqs(1:cutoff),flipud(coh11(1:cutoff,1:cutoff)>=R95),1,'k');
end;
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
image_reverseyaxis;
set(gca,'PlotBoxAspectRatio',[1 1 1]);
title('DF coherence ch.1');
subplot(4,2,3);
imagesc(freqs(1:cutoff),freqs(1:cutoff),flipud(coh22(1:cutoff,1:cutoff)));
colorbar;
if (plotconflimit)
    hold on;
    contour(freqs(1:cutoff),freqs(1:cutoff),flipud(coh22(1:cutoff,1:cutoff)>=R95),1,'k');
end;
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
image_reverseyaxis;
set(gca,'PlotBoxAspectRatio',[1 1 1]);
title('DF coherence ch.2');
subplot(4,2,5);
imagesc(freqs(1:cutoff),freqs(1:cutoff),flipud(coh12(1:cutoff,1:cutoff)));
colorbar;
if (plotconflimit)
    hold on;
    contour(freqs(1:cutoff),freqs(1:cutoff),flipud(coh12(1:cutoff,1:cutoff)>=R95),1,'k');
end;
ylabel('Freq ch.1 (Hz)'); xlabel('Freq ch.2 (Hz)');
image_reverseyaxis;
set(gca,'PlotBoxAspectRatio',[1 1 1]);
title('DF coherence ch.1,2');
subplot(4,2,7);
plot(freqs,abs(diag(coh12)),'k',freqs,R95*ones(size(freqs)),'k--');
xlim(freqs([1 end]));
title('Stationary coherence ch.1,2');
xlabel('Freq (Hz)');

% Plot DF phase
cutoff=find(abs(freqs-maxf)==min(abs(freqs-maxf)));
subplot(4,2,2);
imagesc(freqs(1:cutoff),freqs(1:cutoff),flipud(ph11(1:cutoff,1:cutoff)));
colorbar;
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
image_reverseyaxis;
set(gca,'PlotBoxAspectRatio',[1 1 1]);
title('DF phase ch.1');
subplot(4,2,4);
imagesc(freqs(1:cutoff),freqs(1:cutoff),flipud(ph22(1:cutoff,1:cutoff)));
colorbar;
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
image_reverseyaxis;
set(gca,'PlotBoxAspectRatio',[1 1 1]);
title('DF phase ch.2');
subplot(4,2,6);
imagesc(freqs(1:cutoff),freqs(1:cutoff),flipud(ph12(1:cutoff,1:cutoff)));
colorbar;
ylabel('Freq ch.1 (Hz)'); xlabel('Freq ch.2 (Hz)');
image_reverseyaxis;
set(gca,'PlotBoxAspectRatio',[1 1 1]);
title('DF phase ch.1,2');
subplot(4,2,8);
plot(freqs,diag(ph12),'k',freqs,R95*ones(size(freqs)),'k--');
xlim(freqs([1 end]));
title('Stationary phase ch.1,2');
xlabel('Freq (Hz)');
