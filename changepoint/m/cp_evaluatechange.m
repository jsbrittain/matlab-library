function [m,stats]=cp_evaluatechange(X,p,pint,m0,level,stats)

% Sanity check
if (length(X)<2)
    m=[];
    return;
end;

% Evaluate whether changepoint exists in data
[pchange,S]=cp_bootstrap(X);
if (pchange>p)
    % Find changepoint
    mse=zeros(length(X)-1,1);
    for m=(1:length(X)-1)
        mse(m)=sum((X(1:m)-mean(X(1:m))).^2)+sum((X((m+1):end)-mean(X((m+1):end))).^2);
    end;
    m=find(mse==min(mse),1,'first');
     %cl=cp_conflim(X,pint);
    % Populate stats structure
    stats(1).S=S;
    stats(1).m=m;
    stats(1).pth=p;
    stats(1).pchange=pchange;
    stats(1).level=level;
     %stats(1).cl=cl; 
    % Split data
    X1=X(1:m); X2=X((m+1):end);
    m=m0+m(1);
    % Recursively evaluate subsections for changepoints
    [m1,stats1]=cp_evaluatechange(X1,p,pint,0,(level+1),stats);
    [m2,stats2]=cp_evaluatechange(X2,p,pint,m,(level+1),stats);
    stats=[stats1 stats stats2];
    % Return changepoints
    m=[m1 m m2];
else
    m=[];
    stats=[];
end;

% Remove redundant arguments
if (nargout<2)
    clear('stats');
end;
