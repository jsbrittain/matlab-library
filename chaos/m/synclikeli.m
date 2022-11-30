function [Skij,Ski]=synclikeli(X,E,H,w1,w2)
%
%
%

N=size(X,1);                % Time points
m=size(X,2);                % Embedding dimension
dim=size(X,3);              % Channels

% Synchronisation likelihood (k,i,j)
Skij=zeros(m,N,N);          % Allocate memory
for k=(1:m)
    for i=(1:N)
        for j=(1:N)
            if (sqrt(sum((X(i,:,k)-X(j,:,k)).^2,2)) < E(k,i))
                Skij(k,i,j)=(H(i,j)-1)/(m-1);
            end;
        end;
    end;
end;

% Synchronisation likelihood (k,i)
Ski=zeros(m,N);             % Allocate memory
for k=(1:dim)
    for i=(1:N)
        % Probability distance (Stam 2005, Eq.2)
        jj=(1:N);
        jj=jj(w1<abs(i-jj)); jj=jj(abs(i-jj)<w2);
        Ski(k,i)=sum(Skij(k,i,jj))/(2*(w2-w1));
    end;
end;
