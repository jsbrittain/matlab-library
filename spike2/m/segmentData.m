function hdat = segmentData( hdat, head, seglen )

segk = seglen*head(1).rate;
segcount = floor(size(hdat,1)/segk);
hdat = hdat( 1:segcount*segk, : );

hdat = reshape( hdat, [ segk segcount size(hdat,2) ] );
