function xf_psp( xf11, xf22, xf12, params, maxf )

if (~exist('maxf','var'))
    maxf = [];
end;

figure;

subplot(221);
    imagesc( params.freqs, params.freqs, xf11 );
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    colorbar;
    if (~isempty(maxf)), xlim([0 maxf]); ylim([0 maxf]); end;
    title('Channel 1');
    set(gca,'clim',max(max(abs( (xf11-diag(diag(xf11))) )))*[-1 1]);
    
subplot(222);
    imagesc( params.freqs, params.freqs, xf22 );
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    colorbar;
    if (~isempty(maxf)), xlim([0 maxf]); ylim([0 maxf]); end;
    title('Channel 2');
    set(gca,'clim',max(max(abs( (xf22-diag(diag(xf22))) )))*[-1 1]);
    
subplot(223);
    imagesc( params.freqs, params.freqs, xf12 );
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    colorbar;
    if (~isempty(maxf)), xlim([0 maxf]); ylim([0 maxf]); end;
    ylabel('Channel 1'); xlabel('Channel 2');
    hold('on'); plot( xlim, ylim, 'k', 'linewidth', 2 );
    set(gca,'clim',max(abs(get(gca,'clim')))*[-1 1]);
 