function T=te(x,y)
%function T=te(x,y)
%
% Compute transfer entropy between time-series channels `y -> x`
%
%function T=te(x,y)

% Reserve working memory space
gridpts = max(x);
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
    if ((~isnan(x(n+1))) & (~isnan(x(n))))
        px1x(x(n+1),x(n))=px1x(x(n+1),x(n)) + 1;
    end;
    % p(x_n,y_n)
    if ((~isnan(x(n))) & (~isnan(y(n))))
        pxy(x(n),y(n))=pxy(x(n),y(n)) + 1;
    end;
    % p(x_n+1,x_n,y_n)
    if ((~isnan(x(n+1))) & (~isnan(x(n))) & (~isnan(y(n))))
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
