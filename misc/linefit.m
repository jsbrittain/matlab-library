function [stats,h,x,y] = linefit( x, y, varargin )
%
% Line-fit with confidence intervals
%
% This function takes parameter-value par arguments
%
% Reference
%   http://people.stfx.ca/bliengme/ExcelTips/RegressionAnalysisConfidence2.htm
%   ( a copy of this site with support material is available in the "library/support/linefit" folder)
%

% Created JSB 20/10/2011
% Last modified JSB 24/02/2016

% Input parser
p = inputParser;
addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric);
addParameter(p,'Display',true,@islogical);
addParameter(p,'NumberItems',false,@islogical);
addParameter(p,'MarkerStyle','.',@ischar);
addParameter(p,'MarkerSize',6,@isnumeric);
addParameter(p,'LineWidth',0.5,@isnumeric);
addParameter(p,'Legend','on',@(x) ischar(x) && ismember(lower(x),{'on','off'}))
addParameter(p,'LegendLocation','NorthEast',@ischar);
addParameter(p,'Rank','none',@ischar);       % ascend/descend/none <- order makes no difference to stats, only plot
addParameter(p,'PredictionIntervals',0,@isnumeric); % Currently plots confidence intervals, for prediction intervals set this value to '1'
addParameter(p,'ForceIdentityFit',false,@islogical); % Force fir to identity line (graphical only)
parse(p,x,y,varargin{:});
% Move parsed parameters into function workspace
Display = p.Results.Display;
NumberItems = p.Results.NumberItems;
DisplayLegend = strcmpi(p.Results.Legend,'on');
LegendLocation = p.Results.LegendLocation;
MarkerSize = p.Results.MarkerSize;
LineWidth = p.Results.LineWidth;
MarkerStyle = p.Results.MarkerStyle;
Rank = lower( p.Results.Rank );
PredictionIntervals = p.Results.PredictionIntervals;
ForceIdentityFit = p.Results.ForceIdentityFit;

% Ensure column vectors
x = x(:);
y = y(:);
assert(length(x)==length(y),'Input vectors x and y must be the same length.');

% Rank (fractional ranking method)
if ( ~find(strcmp({'none','ascend','descend'},Rank)) )
    error('Unknown ranking specified');
end;
% Ascend/descend makes no difference to stats, only scatter plot
if ( ~strcmp(Rank,'none') )
    % Fractional ranking ( x )
    [~,ix] = sort(x,Rank); [~,ix] = sort(ix,'ascend'); % Sort of indices must be ascend
    uq = unique(x); for k = (1:length(uq)), ix(x==uq(k)) = mean(ix(x==uq(k))); end; x = ix;
    % Fractional ranking ( y )
    [~,ix] = sort(y,Rank); [~,ix] = sort(ix,'ascend'); % Sort of indices must be ascend
    uq = unique(y); for k = (1:length(uq)), ix(y==uq(k)) = mean(ix(y==uq(k))); end; y = ix;
end;

% Default parameters
alpha = 0.05;
dxfactor = 10;

% Remove NaNs
strIndex = cellfun(@num2str,num2cell(1:length(x)),'UniformOutput',false);
ix = (~isnan(x)) & (~isnan(y)) & (~isinf(x)) & (~isinf(y));
x = x(ix); y = y(ix); strIndex = strIndex(ix);
xmin = min(x);
xmax = max(x);

% Data parameters
n = length(x);
dx = (xmax-xmin)/n/dxfactor;

% Linear regression ( y = mx + c )
C = cov(x,y);                   % Covariance matrix
m = C(1,2)/C(1);                % Gradient
c = mean(y) - m*mean(x);        % Intercept

if ( ForceIdentityFit )
    m = 1;
    c = 0;
    DisplayLegend = false;
    warning('Forcing identity fit --- graphical only --- also, check method!');
end

% Additional values
sxy = std( (m*x+c) - y );       % Std err of estimate
xbar = mean(x);                 % Mean x
tcrit = tinv((1-alpha/2),n-2);    % Critical t-value
ssx = sum((x-xbar).^2);         % Sum-of-squares diff from mean of x

% Determine confidence interval
x0 = (xmin:dx:xmax);
cl = tcrit*sxy.*sqrt(PredictionIntervals+1/n+(x0-xbar).^2/ssx);

% Fill output structure
stats.n = n;
stats.cov = C;
stats.m = m;
stats.c = c;
stats.sxy = sxy;
stats.xbar = xbar;
stats.ssx = ssx;
stats.x0 = x0;
stats.cl = cl;
stats.alpha = alpha;
stats.tcrit = tcrit;
stats.r = stats.cov(1,2)/sqrt(prod(diag(stats.cov)));
stats.R2 = stats.r.^2;

% Plot output
if ( Display )
    h(1) = plot( x, y, 'Color', 'k', 'LineStyle', 'None', 'Marker', MarkerStyle, 'MarkerSize', MarkerSize );    % Scatter plot
    hold('on');
    if ( NumberItems )
        for k = (1:length(x))
            text( x(k), y(k), [' ' strIndex{k}] );
        end;
    end;
    h(2) = plot( x0([1 end]), m*x0([1 end])+c, 'k-', 'LineWidth', LineWidth ); % Line of best fit
    %if ( ~strcmp(Rank,'none') )
        h(3) = plot( x0, (m*x0+c)-cl, 'r-', 'LineWidth', LineWidth );
        h(4) = plot( x0, (m*x0+c)+cl, 'r-', 'LineWidth', LineWidth );
    %end;
end;

% Do regression analysis if available
if (exist('regress','file')==2)
    [b,bint,r,rint,statsR] = regress(x,[ones(size(y)) y]);
    stats.regress.b=b;
    stats.regress.bint=bint;
    stats.regress.r=r;
    stats.regress.rint=rint;
    stats.regress.R2=statsR(1);
    stats.regress.F=statsR(2);
    stats.regress.p=statsR(3);
    stats.regress.errvar=statsR(4);
    if ( Display )
        h2 = plot(nan,nan,nan,nan,nan,nan);
        if ( DisplayLegend )
            if ( stats.regress.p < 0.001 )
                pstr = sprintf('p < 0.001');
            else
                pstr = sprintf('p = %1.3f',stats.regress.p);
            end
            legend(h2,sprintf('R^2 = %1.3f',stats.regress.R2),sprintf('F(%g,%g) = %1.3f',1,stats.n-2,stats.regress.F),pstr,'Location',LegendLocation);
        end
        if ( stats.regress.p > 0.05 )
            delete(h([2:end]));
            h([2:end]) = [];
            %set(h(2),'linestyle','--');
        end;
    end;
end;
