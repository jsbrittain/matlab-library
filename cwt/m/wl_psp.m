function wl_psp(time,sp11,sp22,sp12,wlparam)
%function wl_psp(time,W,wlparam)
%function wl_psp(time,sp11,sp22,sp12,wlparam)

if ( nargin == 5 )
    W = cat(3,sp11,sp22,sp12,abs(sp12).^2./(sp11.*sp22),angle(sp12));
end;

figure;
h = subplot(2,2,1);
    %scalogram(time,log10(W(:,:,1)),wlparam);
    imagesc( time, wlparam.freqs, log10(W(:,:,1)) );
    colorbar;
    title('Autospectra channel 1');
h(2) = subplot(2,2,2);
    %scalogram(time,log10(W(:,:,2)),wlparam);
    imagesc( time, wlparam.freqs, log10(W(:,:,2)) );
    colorbar;
    title('Autospectra channel 2');
h(3) = subplot(2,2,3);
    %scalogram(time,log10(W(:,:,3)),wlparam);
    imagesc( time, wlparam.freqs, log10(abs(W(:,:,3))) );
    colorbar;
    title('| Crossspectra |');
h(4) =subplot(2,2,4);
    %scalogram(time,W(:,:,4),wlparam);
    imagesc( time, wlparam.freqs, W(:,:,4) );
    colorbar;
    title('Coherence');
linkaxes(h,'xy');
