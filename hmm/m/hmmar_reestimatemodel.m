function [prior,A,hmmar,hmms2]=hmmar_reestimatemodel(y,prior,A,step,hmmar,hmms2)
%
% Re-estimate HMM-AR model given initial estimate of (prior, A),
% measurement sequence and a set of AR coefficients.
%
% Input parameters
%       y           Measurement sequence
%                     ( 1 x N )
%       prior       
%       A           
%       step        Scale estimate every `step' samples (0=none; 1=default)
%       hmmar       
%       hmms2       
%
% Output parameters
%       prior, A, hmmar, hmms2      Updated parameters
%

% Check input parameters
validateaccuracy = false;

% Determine data parameters
N = length(y);
Q = length(prior);
p = size(hmmar,1);

% Ensure each state has an associated variance estimate
if ((p>1) && (length(hmms2)==1))
    hmms2 = hmms2 * ones(1,p);
end;

% Reserve variable space
X = zeros(N,p);

% Compute probability of measurement being in each state
[b,yb] = ar_obslik(y,hmmar,hmms2);
logp0 = chmm_probobs(b(:,(p+2:N)),prior,A,step);

% Compute matrix X
for t=((p+2):N)
    X(t,:) = -(y((t-1):(-1):(t-p)));
end;

% Iteration threshold (percentage change log p)
th=1e-6;
maxiter=100;

% Iterate re-estimation procedure
n=1; logp=0;
while (((1-logp/logp0)>th) && (n<maxiter))
    
    % Iterate model parameters
    [prior,A,gammati] = chmm_estimatemodel(b(:,(p+2:N)),prior,A,step);
    
    % Update AR coefficients (See Penny & Roberts 1999, p.487)
    for j=(1:Q)
        % Phrase AR coefficient update in least-squares form
        C = sparse(diag(sqrt(gammati(:,j))));
        Y = y(p+2:N-1).';
        Xtilde = C*X(p+2:N-1,:);
        Ytilde = C*Y;
        
        % Solve least-squares using (more efficient) SVD
        %hmmar(:,j) = inv((Xtilde.')*Xtilde)*(Xtilde.')*Ytilde;
        [U,S,V] = svd(Xtilde,'econ');
        hmmar(:,j) = V*pinv(S)*(U.')*Ytilde;
        
        hmms2(j) = sum((gammati(:,j).').*((y(p+2:N-1)-yb(j,p+2:N-1)).^2))./sum(gammati(:,j));
    end;
    
    % Track fitting
    [b,yb] = ar_obslik(y,hmmar,hmms2);
    logp0 = logp;
    logp  = chmm_probobs(b(:,(p+2:N)),prior,A,step);
    
    % Display iteration progress
    disp(['Step ' num2str(n) ' (log p=' num2str(logp) ', ' num2str((1-logp/logp0)/(100*th)) ')']);
    
    % Increment loop counter
    n=n+1;
    
end;
