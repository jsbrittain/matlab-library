function b=mt_dk(Sk,sig2,lambda)
%function b=mt_dk(Sk,sig2,lambda)
%
% Generate adaptive weights for expansion coefficients
% (Algorithm: Percival et al. 1998, p.368-370)
%
% Input parameters
%   Sk      Eigenspectra (freqs,k)
%   sig2    Variance of data
%   lambda  Slepian eigenvalues
%
% Output parameter
%   b       Adaptive weights for expansion coefficients
%
% Reference
%   Percival et al. (1998) Spectral Analysis for Physical Applications:
%      Multitaper and Conventional Univariate Techniques.  Cambridge
%      University Press, reprint with corrections.  Lib QA320.P434
%      ISBN 0-521-43541-2 (paperback).
%
%function b=mt_dk(Sk,sig2,lambda)

% Check inputs
if (~isreal(sig2))
    error(' Variance must be real-valued');
end;

% Initiate variables
fcount=size(Sk,1);
K=length(lambda);
S=(Sk(:,1)+Sk(:,2))/2;
tol=.0005*sig2/fcount;      % Tolerance
i=0;                        % Iteration count
a=sig2*(1-lambda);          % Broad-band bias

% Perform iteration
S1=zeros(fcount,1);
while sum(abs(S-S1)/fcount)>tol
	i=i+1;
	% Calculate weights
	b=(S*ones(1,K))./(S*lambda'+ones(fcount,1)*a'); 
	% Calculate new spectral estimate
	dk=(b.^2).*(ones(fcount,1)*lambda');
	S1=sum(dk'.*Sk')./ sum(dk');                        % (370a)
	S1=S1';
	Stemp=S1; S1=S; S=Stemp;  % swap S and S1
end;
