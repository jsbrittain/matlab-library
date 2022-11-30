X=[10.7 13.0 11.4 11.5 12.5 14.1 14.8 14.1 12.6 16.0 11.7 10.6 10.0 11.4 7.9 9.5 8.0 11.8 10.5 11.2 9.2 10.1 10.4 10.5];

% Ref:
%  http://www.variation.com/cpa/tech/changepoint.html

p_level=(1-0.1);
p_int=0.05;

% CP on data
[m,stats]=changepoint(X,p_level);
%[m,stats]=changepoint(X,p_level,p_int);

return;

% CP on differences
D=diff(X); D=D(1:2:end); % (subsample to avoid correlated elements)
[m,S]=changepoint(D,p_level);

% Outlier (CP on data)
X(6)=25;
[m,S]=changepoint(X,p_level);

% Outlier (CP on ranked data)
[x,IX]=sort(X);
[m,S]=changepoint(IX,p_level);
