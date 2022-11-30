function H=hdist(X,E,w1,w2)
%
% Redundant - much faster if computed via pdist_esearch
%

N=size(X,1);                % Time points
m=size(X,2);                % Embedding dimension
dim=size(X,3);              % Channels

H=zeros(N);                 % Allocate memory
for i=(1:N)
    for j=(1:N)
        if ((w1<abs(i-j)) && (abs(i-j)<w2))
            for k=(1:dim)
                H(i,j)=H(i,j)+heaviside2(E(k,i)-sqrt(sum((X(i,:,k)-X(j,:,k)).^2,2)));
            end;
        end;
    end;
end;
