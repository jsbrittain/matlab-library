function [m,stats]=changepoint(X,plevel,pint)

if (~exist('pint'))
    pint=[];
end;

[m,stats]=cp_evaluatechange(X,plevel,pint,0,1,[]);

if (nargout<2)
    clear('stats');
end;
