function I = mi( pxy )
%function I = mi( pxy )
%
% Mutual Information
%
% takes histogram
%   p(x,y)
% computes
%   I(X;Y)
%
%function I = mi( pxyz )

% Marginal probabilities
pxy = pxy/sum(pxy(:));   % Ensure normalisation
px = sum(pxy,2);
py = sum(pxy,1);

% Mutual information
I = sum(sum( pxy.*log(pxy./repmat(px,[1 size(pxy,2)])./repmat(py,[size(pxy,1) 1])) ,1),2);
