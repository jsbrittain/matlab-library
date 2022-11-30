function [H,Hmax,p]=shannonentropy(x,y,gridpts)
%function [H,Hmax,p]=shannonentropy(x,[y,[gridpts]])
%
% Compute Shannon entropy between 2 processes (x,y) on a grid of `gridpts'
% levels.
%
% `gridpts' can be scalar or vector. ensure vector length not equal to `x'
% length to avoid ambiguation with `y'.
%
%function [H,Hmax,p]=shannonentropy(x,y,[gridpts])

% Check input parameters
if (nargin<3)
    if (nargin<2)
        y=x; gridpts=[];
    else
        if (length(x)~=length(y))
            gridpts=y;
            y=x;
        end;
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

    % Empirical distribution
    p=hist(x,gridpts);
else
    % Empirical distribution (and remove endpoints)
    dp=gridpts(2)-gridpts(1);
    p=hist(x,[gridpts(1)-dp gridpts gridpts(end)+dp]);
    p=p(2:end-1);
end;

% Compute entropy
p(p==0)=1;
p=p/sum(p);
H=-sum(p.*log(p));

% Calculate maximum entropy
bincount=length(p);
Hmax=-sum(repmat(1/bincount,bincount,1).*log(repmat(1/bincount,bincount,1)));

% Check output arguments
if (nargout<3)
    clear('p');
end;
if (nargout<2)
    clear('hmax');
end;
