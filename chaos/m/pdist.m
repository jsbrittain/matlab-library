function P=pdist(X,epsilon,w1,w2)

N=size(X,1);                % Time points
m=size(X,2);                % Embedding dimension
dim=size(X,3);              % Channels

if (length(epsilon)==1)
    E=epsilon*ones(dim,N);
else
    E=epsilon;
end;
clear('epsilon');

P=zeros(dim,N);             % Allocate memory
for k=(1:dim)
    for i=(1:N)
        % Probability distance (Stam 2005, Eq.2)
        jj=(1:N); jj=jj(w1<abs(i-jj)); jj=jj(abs(i-jj)<w2);
        Xi=repmat(X(i,:,k),length(jj),1);
        P(k,i)=sum(((E(k,i)-sqrt(sum((Xi-X(jj,:,k)).^2,2)))>0),1)/(2*(w2-w1));
    end;
end;
