function [Zhist,xgrid,ygrid,Zmean,Zssd] = image_density( x, y, bincount, kernel, quantity )
%function [Z,xgrid,ygrid] = image_density( x, y, bincount, kernel, quantity )
%
% quantity is an optional vector of length(x) that takes the average value
% per pixel instead of the histogram measure as normal.
%
% NB: `quantity' changes the metric from a histogram measure to an average
% measure, therefore providing quantity = ones(size(x)) is NOT the same as
% setting quantity = [];
%
% Option returns:
%   Zhist - histogram
%   Zmean - mean based on quantity vector
%   Zssd  - ssd based on quantity vector
%
%function [Z,xgrid,ygrid] = image_density( x, y, bincount, kernel, quantity )

binnormalise = true;
if (~exist('bincount','var'))
    bincount = [];
end
if (~exist('kernel','var'))
    kernel = [];
end
if (~exist('quantity','var'))
    quantity = [];
end

if (isempty(kernel))
    kernel = 1;
end
if (isempty(bincount))
    bincount = 20;
end
if (isempty(quantity))
    binnormalise = false;
    quantity = ones(size(x));
end

% If bincount is a 2-column vector then (x,y) grids set
if ( ~isscalar( bincount ) )
    minx = min(bincount(:,1));
    maxx = max(bincount(:,1));
    x(x<minx) = NaN; x(x>maxx) = NaN;
    miny = min(bincount(:,2));
    maxy = max(bincount(:,2));
    y(y<miny) = NaN; y(y>maxy) = NaN;
    bincount = size(bincount,1);
else
    minx = min(x); maxx = max(x);
    miny = min(y); maxy = max(y);
end
rangex = maxx - minx;
rangey = maxy - miny;

kk = (~isnan(x)) & (~isnan(y));
x = x(kk); y = y(kk); quantity = quantity(kk);

xk = round( (bincount)*(x-minx)/rangex )+1; xk(xk==(bincount+1)) = 1;
xk(xk<1) = 1; xk(xk>length(xk)) = length(xk)-1;
yk = round( (bincount)*(y-miny)/rangey )+1; yk(yk==(bincount+1)) = 1;
yk(yk<1) = 1; yk(yk>length(yk)) = length(yk)-1;
histZ = zeros(size(xk));

xgrid = linspace(minx,maxx,bincount);
ygrid = linspace(miny,maxy,bincount);

Zmean = zeros( bincount );
Zssd = Zmean;       % squared sum of differences; better to iterate than var
Zhist = Zmean;
for k = (1:length(x))
    if isnan(yk(k)) || isnan(xk(k))
        continue;
    end
    
    % Iterative mean
    Zmeanp = Zmean(yk(k),xk(k));
    Zmean(yk(k),xk(k)) = (Zhist(yk(k),xk(k))*Zmean(yk(k),xk(k)) + quantity(k))/(Zhist(yk(k),xk(k))+1);
    Zssd(yk(k),xk(k)) = Zssd(yk(k),xk(k)) + (quantity(k)-Zmeanp)*(quantity(k)-Zmean(yk(k),xk(k)));
    Zhist(yk(k),xk(k)) = Zhist(yk(k),xk(k)) + 1;
end
Zssd = sqrt(Zssd./Zhist);         % Sample std dev

if (exist('convolve2','file')~=0)
    Z = convolve2( Zhist, kernel/sum(kernel(:)), 'wrap' );
else
    Z = conv2( Zhist, kernel/sum(kernel(:)), 'same' );
end

if (nargout==0)
    imagesc( [minx maxx], [miny maxy], Z );
    set(gca,'ydir','normal');
end
