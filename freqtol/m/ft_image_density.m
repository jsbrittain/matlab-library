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
    bincount = 20;
end;
if (~exist('kernel','var'))
    kernel = 1;
end;
if (~exist('quantity','var'))
    quantity = [];
end;

if (isempty(quantity))
    binnormalise = false;
    quantity = ones(size(x));
end;

xk = round( (bincount)*(x-min(x))/range(x) )+1; xk(xk==(bincount+1)) = 1;
yk = round( (bincount)*(y-min(y))/range(y) )+1; yk(yk==(bincount+1)) = 1;
histZ = zeros(size(xk));

xgrid = linspace(min(x),max(x),bincount);
ygrid = linspace(min(y),max(y),bincount);

Zmean = zeros( bincount );
Zssd = Zmean;       % squared sum of differences; better to iterate than var
Zhist = Zmean;
for k = (1:length(x))
    if isnan(yk(k)) || isnan(xk(k))
        continue;
    end;
    
    % Iterative mean
    Zmeanp = Zmean(yk(k),xk(k));
    Zmean(yk(k),xk(k)) = (Zhist(yk(k),xk(k))*Zmean(yk(k),xk(k)) + quantity(k))/(Zhist(yk(k),xk(k))+1);
    Zssd(yk(k),xk(k)) = Zssd(yk(k),xk(k)) + (quantity(k)-Zmeanp)*(quantity(k)-Zmean(yk(k),xk(k)));
    Zhist(yk(k),xk(k)) = Zhist(yk(k),xk(k)) + 1;
end;
Zssd = sqrt(Zssd./Zhist);         % Sample std dev

if (exist('convolve2','file')~=0)
    Z = convolve2( Zhist, kernel/sum(kernel(:)), 'wrap' );
else
    Z = conv2( Zhist, kernel/sum(kernel(:)), 'same' );
end;

if (nargout==0)
    imagesc( [min(x) max(x)], [min(y) max(y)], Z );
    set(gca,'ydir','normal');
end;
