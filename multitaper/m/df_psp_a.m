function df_psp_a(DF11,DF22,DF12,params,maxf,smoothing_kernel)
%function df_psp_a(DF11,DF22,DF12,params,maxf,smoothing_kernel)
%
% Summary plot of dual-frequency spectra between two processes
%
% Input parameters
%       DF11        Dual-frequency auto-spectra ch.1
%       DF22        Dual-frequency auto-spectra ch.2
%       DF12        Dual-frequency cross-spectra ch.1,2
%       params      Parameters structure
%
%function df_psp_a(df11,df22,df12,params,maxf,smoothing_kernel)

if (~exist('maxf','var'))
    maxf=[];
end;
if (~exist('smoothing_kernel','var'))
    smoothing_kernel = 1;
end;
smoothing_kernel = smoothing_kernel/sum(smoothing_kernel);

% Extract stationary spectra
S11=diag(DF11);
S22=diag(DF22);
freqs=params.freqs;

% Plot dual-frequency spectra
figure;
subplot(3,2,1);
    imagesc(freqs,freqs,convolve2(log10(abs(DF11)),smoothing_kernel,'same'));
    set(gca,'ydir','normal');
    xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    if (~isempty(maxf)), xlim([0 maxf]); ylim([0 maxf]); end;
    title('Auto-DF-ch.1');
subplot(3,2,2);
    plot(freqs,10*log10(diag(DF11)));
    xlim(freqs([1 end]));
    title('Stationary Auto-DF ch.1');
    if (~isempty(maxf)), xlim([0 maxf]); end;
    ylabel('dB');
subplot(3,2,3);
    imagesc(freqs,freqs,convolve2(log10(abs(DF22)),smoothing_kernel,'same'));
    set(gca,'ydir','normal');
    xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    if (~isempty(maxf)), xlim([0 maxf]); ylim([0 maxf]); end;
    title('Auto-DF-ch.2');
subplot(3,2,4);
    plot(freqs,10*log10(diag(DF22)));
    xlim(freqs([1 end]));
    if (~isempty(maxf)), xlim([0 maxf]); end;
    title('Stationary Auto-DF ch.2');
    ylabel('dB');
subplot(3,2,5);
    imagesc(freqs,freqs,convolve2(log10(abs(DF12)),smoothing_kernel,'same'));
    ylabel('Freq ch.2 (Hz)'); xlabel('Freq ch.1 (Hz)');
    set(gca,'ydir','normal');
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    if (~isempty(maxf)), xlim([0 maxf]); ylim([0 maxf]); end;
    title('|Cross-DF|');
subplot(3,2,6);
    plot(freqs,10*log10(abs(diag(DF12))));
    xlim(freqs([1 end]));
    if (~isempty(maxf)), xlim([0 maxf]); end;
    title('Stationary |Cross-DF|');
    xlabel('Freq (Hz)');
    ylabel('dB');
