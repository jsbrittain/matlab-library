function T=transferentropy(x,y,gridpts)
%function T=transferentropy(x,y,[gridpts])
%
% Compute transfer entropy between time-series channels `y -> x` by
% discretising both series to a finite grid (default 24 levels) over 3*sd
% of each channel (outliers rejected).
%
% Presently assumes `gridpts' is a linear vector and therefore only
% actually considers it's min, max and length properties.
%
%function T=transferentropy(x,y,[gridpts])

% Check input parameters
if (nargin<3)
    if (length(x)~=length(y))
        gridpts=y;
        y=x;
    end;
end;
if (~exist('gridpts'))
    gridpts=[];
end;
if (isempty(gridpts))
    gridpts=24;
end;

if (length(gridpts)==1)

    % Power (and stdev) normalise data
    x=x-mean(x);
    x=x/sqrt(mean(x.^2));
    y=y-mean(y);
    y=y/sqrt(mean(y.^2));
    
    % Discretise time-series to grid
    x=round(gridpts*(x+3)/6);       % Discretise onto +/- 3 sd grid
    x(x<1)=NaN; x(x>gridpts)=NaN;   % Dismiss outliers
    y=round(gridpts*(y+3)/6);
    y(y<1)=NaN; y(y>gridpts)=NaN;
    
else
    
    % Discretise onto linear vector grid

    % Signal x
    x=x-min(gridpts); x=x/max(gridpts);
    x=round(x*length(gridpts));
    x(x<1)=NaN; x(x>length(gridpts))=NaN;   % Dismiss outliers
    
    % Signal y
    y=y-min(gridpts); y=y/max(gridpts);
    y=round(y*length(gridpts));
    y(y<1)=NaN; y(y>length(gridpts))=NaN;   % Dismiss outliers
    
    gridpts=length(gridpts);
    
end;

% Reserve working memory space
px1xy = zeros(gridpts,gridpts,gridpts);
px1x = zeros(gridpts,gridpts);
pxy = zeros(gridpts,gridpts);
px = zeros(gridpts,1);
py = zeros(gridpts,1);

% Traverse time-series
for n=(1:length(x)-1)
    % p(x_n)
    if (~isnan(x(n)))
        px(x(n))=px(x(n)) + 1;
    end;
    % p(y_n)
    if (~isnan(y(n)))
        py(y(n))=py(y(n)) + 1;
    end;
    % p(x_n+1,x_n)
    if ((~isnan(x(n+1))) && (~isnan(x(n))))
        px1x(x(n+1),x(n))=px1x(x(n+1),x(n)) + 1;
    end;
    % p(x_n,y_n)
    if ((~isnan(x(n))) && (~isnan(y(n))))
        pxy(x(n),y(n))=pxy(x(n),y(n)) + 1;
    end;
    % p(x_n+1,x_n,y_n)
    if ((~isnan(x(n+1))) && (~isnan(x(n))) && (~isnan(y(n))))
        px1xy(x(n+1),x(n),y(n))=px1xy(x(n+1),x(n),y(n)) + 1;
    end;
end;

% Ensure probability distributions sum=1
px=px/sum(px);
py=py/sum(py);
pxy=pxy/sum(pxy(:));
px1x=px1x/sum(px1x(:));
px1xy=px1xy/sum(px1xy(:));

% Calculate transfer entropy statistic

warning off

T=px1xy.*log2(px1xy.*reshape(repmat(px',gridpts,gridpts),gridpts,gridpts,gridpts)./permute(reshape(repmat(pxy,1,gridpts),gridpts,gridpts,gridpts),[3 1 2])./reshape(repmat(px1x,1,gridpts),gridpts,gridpts,gridpts));
T=T(~isnan(T) & ~isinf(T));
T=sum(T);

warning on
