function X = hmm_rownorm(X)

X = X ./ repmat(sum(X,2),1,size(X,2));
