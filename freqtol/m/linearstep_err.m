function err = linearstep_err( coeffs, x, y, ystd )

if (~exist('ystd','var'))
    ystd = [];
end;

if (isempty(ystd))
    err = nansum( abs(y - linearstep( coeffs, x )).^2 );
else
    y(ystd==0) = NaN;
    %ystd(ystd==0) = 1;
    err = nansum( ((y - linearstep( coeffs, x )).^2)./(ystd.^2) );
end;
