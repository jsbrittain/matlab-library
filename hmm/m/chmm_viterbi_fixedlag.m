function [q,sigmati] = chmm_viterbi_fixedlag(b,prior,A,lag)
%function [q,[sigmati]] = chmm_viterbi_fixedlag(b,prior,A,lag)
%
% Continuous observation (pdf) Hidden-Markov Model (HMM)
%
% Viterbi algorithm
%
% Returns the most likely set of states for a given model and observation
% set (dynamic programming approach; non-unique).
%
% Input parameters
%       b             Observation likelihood
%                       ( M x N )
%       prior         Initial state distribution
%                       ( Q x 1 )
%       A             Transition probability matrix
%                       ( Q x Q )
%       lag           Lag (samples)
%                       ( scalar )
%
% Output argument
%       q             State-sequence
%                       ( 1 x N )
%
% This routine has been modified to propagate log P throughout the
% calculation. Also, scaling is applied to the sequence probability at each
% sample to facilitate continous update (log P propagation reset). Scaling
% alone is sufficient, but can be muted with an internal switch, thus log
% propagation also remains.
%
% Reference
%   Rabiner, K. R. (1989) A Tutorial on Hidden Markov Models and Selected
%   Applications in Speech Recognition. Proceedings of the IEEE, 77(2):
%   257-286.
%
%function [q,[sigmati]] = chmm_viterbi_fixedlag(b,prior,A,lag)

% Determine data parameters
N=size(b,2);
Q=size(A,1);

% Reserve memory
q = zeros(1,N);
sigmati = zeros(N,Q);
phiti = zeros(N,Q);
logA = log(A);

% Initialise
for i=(1:Q)
    sigmati(1,i) = log(prior(i)) + log(b(i,1));
    phiti(1,i) = 0;
    [dummy,q(1)] = max(sigmati(1,:));
end;

% Induction
for t=(2:N)
    for j=(1:Q)
        [sigmati(t,j),phiti(t,j)] = max( sigmati(t-1,:) + logA(:,j).' );
        sigmati(t,j) = sigmati(t,j) + log(b(j,t));
    end;
    if (true)
        % Normalise sigmati: P()
        sigmati(t,:) = log(exp(sigmati(t,:))/sum(exp(sigmati(t,:))));
    end;
    
    % Fixed-lag backtracking
    [dummy,qk] = max(sigmati(t,:));
    % Path (state-sequence) backtracking
    if (lag>0)
        for k=((t-1):(-1):max(1,t-lag))
            qk = phiti(k+1,qk);
        end;
    end;
    q(max(1,t-lag)) = qk;
end;
