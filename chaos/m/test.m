% test.m

% Construct surrogate data
switch (1)
    case 1,
        x=[cos(0:0.1:32)' randn(321,1)];
        l=1; m=2; w1=10; w2=100; Pref=0.05;
    case 2,
        x=[cos((0:0.1:32)') cos((0:0.1:32)')];
        l=1; m=2; w1=10; w2=100; Pref=0.05;
    case 3,
        % Henon system
        %N=5096; B=0.3; C=1.0;
        N=5096; B=0.1; C=0.6;
        x=rand(1,2); y=x;%rand(1,2);
        for i=(2:N-1)
            x(i+1)=1.4-x(i)^2+0.3*x(i-1);
            y(i+1)=1.4-(C*x(i)+(1-C)*y(i))*y(i)+B*y(i-1);
        end;
        x=[x(1001:end)' y(1001:end)'];
        % Embedding parameters
        l=1; m=10; w1=100; w2=410; Pref=0.05;
    otherwise
        error(' No data specified');
end;

% Lag embedding
X=lagembed(x,l,m);

% Determine epsilon from probability distance
tic
[E,H,Skij,Ski]=pdist_esearch(X,w1,w2,Pref,false);
toc

% Check epsilson computation
P=pdist(X,E,w1,w2);

% Mutual information (time-dependent) for bivariate processes
if (size(x,2)==2)
    Mi=log2(Ski/Pref);
end;
