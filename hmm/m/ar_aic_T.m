function [p,porder] = ar_aic_T(y,prange,T)

N=length(y);
pmin = prange(1); pmax=prange(2);

% Determine measurement noise (static moving window estimates)
disp(['Estimating model order over period T = ' num2str(T) ' samples ...']);
npts = (0:T:(N-T));
porder = zeros(length(npts),(pmax-pmin+1));
for n = (1:length(npts))
    pn = 1;
    for p = (pmin:pmax)
        m = arx( y( npts(n) + (1:T) ), p );
        porder(n,pn) = aic(m);
        pn = pn + 1;
    end;
end;
[~,pordermin] = min(porder,[],2);
p = ceil(median(pordermin)) + pmin - 1;

disp(['Resolved model order p = ' num2str(p)]);
