function bi_psp( p11,bs11,p22,bs22,p12,bs12,bs21,params,bicoh11,bicoh22,bicoh12,bicoh21, fmax )
%function bi_psp( p11,bs11,p22,bs22,p12,bs12,bs21,params,bicoh11,bicoh22,bicoh12,bicoh21, fmax )
%
%
%function bi_psp( p11,bs11,p22,bs22,p12,bs12,bs21,params,bicoh11,bicoh22,bicoh12,bicoh21, fmax )

if (~exist('fmax'))
    fmax = [];
end;

ylims = [params.bifreqs([2 end])];
if (~isempty(fmax))
    ylims(2) = dsearchn(params.freqs',fmax);
end;


figure;

subplot(4,5,1);
    plot( 10*log10(p11), params.freqs, 'k' );
    ylim(ylims);
    ylabel('Channel 1');

subplot(4,5,[2 3]);
    imagesc(params.bifreqs,params.bifreqs, 10*log10( abs(bs11) ) );
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    colorbar;

subplot(4,5,[4 5]);
    imagesc(params.bifreqs,params.bifreqs,abs(bicoh11));
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    set(gca,'clim',[0 max(get(gca,'clim'))]);
    colorbar;

    
    
subplot(4,5,5+1);
    plot( 10*log10(p22), params.freqs, 'k' );
    ylim(ylims);
    ylabel('Channel 2');

subplot(4,5,5+[2 3]);
    imagesc(params.bifreqs,params.bifreqs,10*log10( abs(bs22) ));
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    colorbar;

subplot(4,5,5+[4 5]);
    imagesc(params.bifreqs,params.bifreqs,abs(bicoh22));
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    set(gca,'clim',[0 max(get(gca,'clim'))]);
    colorbar;

    
    
    
subplot(4,5,10+1);
    plot( 10*log10(abs(p12)), params.freqs, 'k' );
    ylim(ylims);
    ylabel('Channel 1,2');

subplot(4,5,10+[2 3]);
    imagesc(params.bifreqs,params.bifreqs,10*log10( abs(bs12) ) );
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    colorbar;

subplot(4,5,10+[4 5]);
    imagesc(params.bifreqs,params.bifreqs,abs(bicoh12));
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    set(gca,'clim',[0 max(get(gca,'clim'))]);
    colorbar;


    
    
subplot(4,5,15+1);
    plot( 10*log10(abs(p12)), params.freqs, 'k' );
    ylim(ylims);
    ylabel('Channel 2,1');

subplot(4,5,15+[2 3]);
    imagesc(params.bifreqs,params.bifreqs,10*log10( abs(bs21) ) );
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    colorbar;

subplot(4,5,15+[4 5]);
    imagesc(params.bifreqs,params.bifreqs,abs(bicoh21));
    set(gca,'ydir','normal');
    axis equal; axis square; axis tight;
    xlim(ylims); ylim(ylims);
    set(gca,'clim',[0 max(get(gca,'clim'))]);
    colorbar;
