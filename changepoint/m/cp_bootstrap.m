function [p,S]=cp_bootstrap(X,Sdiff)

% Determine (mean subtracted) cumulative sum and range
S=[0 cumsum(X-mean(X))];
Sdiff=max(S)-min(S);

% Bootstrap
B=1000; N=0;
for ind=(1:B)
    % Randomly permute vector (sample without replacement)
    X0=X(randperm(length(X)));
    S0=[0 cumsum(X0-mean(X0))];
    S0diff=max(S0)-min(S0);
    % Counting statistic ( S0diff < Sdiff )
    if (S0diff < Sdiff)
        N=N+1;
    end;
end;

% Determine confidence limit based on bootstraps
p=N/B;
