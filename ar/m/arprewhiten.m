function y = arprewhiten( x, p )

if (~exist('p'))
    p = [];
end;
if (isempty(p))
    p = 1;
end;

if (~exist('arfit'))
    error(' Routine requires third party ARFIT');
end;

if (p~=1)
    error(' Only supports AR-1 models at moment (due to simple reconstruction)');
end;

[w,A,C,SBC,FPE,th] = arfit(x,p,p);

y = [ x(1); x(2:end)-A*x(1:end-1) ];
