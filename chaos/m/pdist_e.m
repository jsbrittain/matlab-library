function E=pdist_e(X,w1,w2,Pref)

N=size(X,1);                % Time points
m=size(X,2);                % Embedding dimension
dim=size(X,3);              % Channels

eplist=(0.01:0.001:5);

E=zeros(dim,N);             % Allocate memory
e=zeros(length(eplist),1);
for k=(1:dim)
    for i=(1:N)
        % Probability distance (Stam 2005, Eq.2)
        jj=(1:N);
        jj=jj(w1<abs(i-jj)); jj=jj(abs(i-jj)<w2);
        Xi=repmat(X(i,:,k),length(jj),1);
        
        % Epsilon search
        sq=sqrt(sum((Xi-X(jj,:,k)).^2,2));
        for ind=(1:length(eplist))
            e(ind)=sum(heaviside2(eplist(ind)-sq),1);
        end;
        E(k,i)=eplist(dsearchn(e,Pref));
        
    end;
end;
