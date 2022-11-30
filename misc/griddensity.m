function [I,xx,yy] = griddensity( x, y, N, M, conv_kernel, uniqueonly )
%function [I,xx,yy] = griddensity( x, y, N, M, conv_kernel, uniqueonly )

if (~exist('N','var'))
    N = [];
end;
if (~exist('M','var'))
    M = [];
end;
if (~exist('conv_kernel','var'))
    conv_kernel = [ 0 1 0; 1 1 1; 0 1 0 ];
end;
if (~exist('uniqueonly','var'))
    uniqueonly = false;
end;

if (isempty(N))
    N = 100;
end;
if (isempty(M))
    M = N;
end;

xx = min(x):(max(x)-min(x))/(N-1):max(x);
yy = min(y):(max(y)-min(y))/(M-1):max(y);

xn = x - min(x);
xn = xn/max(xn);
xn = round(xn*(N-1))+1;

yn = y - min(y);
yn = yn/max(yn);
yn = round(yn*(N-1))+1;

I = zeros(N,M);
for n = (1:length(x))
    if (~isnan(xn(n)+yn(n)))
        I(xn(n),yn(n)) = I(xn(n),yn(n)) + 1;
    end;
end;
if (uniqueonly)
    I = double(I>0);
end;

% Plot
if (nargout==0)
    conv_kernel = conv_kernel/sum(conv_kernel(:));
    contourf( xx, yy, convolve2( I, conv_kernel, 'same' )' );
    hold on;
    plot( x, y, 'w.' );
end;
