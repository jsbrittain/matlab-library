function be=mwMaxBeta(ga,A);
%function be=mwMaxBeta(ga,A);
%
% Multi-wavelet routine to determine the maximum possible beta
%
%function be=mwMaxBeta(ga,A);

maxb=1000;

for be=1:maxb
    r=(2*be+1)/ga;
    C=ga*A*gamma(r)^2/(2*pi*gamma(r+1-1/ga)*gamma(r+1/ga))+1;
    if (isnan(C))
        break;
    end;
end;
be=be-1;
if (be==0)
    error(' Unable to determine beta with these parameters');
end;
