function err = trilinearstep_err( coeffs, x, y, ystd )

if (~exist('ystd'))
    ystd = [];
end;

if (isempty(ystd))
    err = nansum( abs(y - trilinearstep( coeffs, x )).^2 );
else
    y(ystd==0) = NaN;
    %ystd(ystd==0) = 1;
    err = nansum( ((y - trilinearstep( coeffs, x )).^2)./ystd );
end;
