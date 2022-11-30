function x=invbetainc(y,a,b);
%function x=invbetainc(y,a,b);
%
% Inverse Incomplete Beta Function
%
% Uses function minimisation to numerically determine the inverse of the
% incompelte beta function (betainc - Matlab routine).
%
%function x=invbetainc(y,a,b);

% Pass numeric routine 'y' as the cutoff parameter
[x,fval,exitflag]=fminbnd('betainc_cutoff',0,1,[],a,b,y);

% Check results
if (exitflag<=0)
    warning(' Inverse incomplete beta function did not converge');
end;
