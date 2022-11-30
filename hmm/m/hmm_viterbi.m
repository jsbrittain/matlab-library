function [q,logP] = hmm_viterbi(O,prior,A,B)
%function [q,[logP]] = hmm_viterbi(O,prior,A,B)
%
% Viterbi algorithm
%
% Returns the most likely set of states for a given model and observation
% set (dynamic programming approach; non-unique).
%
% Input parameters
%       O             Observation sequence
%                       ( 1 x N )
%       prior         Initial state distribution
%                       ( Q x 1 )
%       A             Transition probability matrix
%                       ( Q x Q )
%       B             Emission probability matrix
%                       ( Q x M )
%
% Output argument
%       q             State-sequence
%                       ( 1 x N )
%
% This routine has been modified to propagate log P throughout the
% calculation.
%
% Reference
%   Rabiner, K. R. (1989) A Tutorial on Hidden Markov Models and Selected
%   Applications in Speech Recognition. Proceedings of the IEEE, 77(2):
%   257-286.
%
%function [q,[logP]] = hmm_viterbi(O,prior,A,B)

% Determine data parameters
N=length(O);
Q=size(A,1);

% Reserve memory
sigmati = zeros(N,Q);
phiti = zeros(N,Q);

% Initialise
for i=(1:Q)
    %sigmati(1,i) = prior(i)*B(i,O(1));
    sigmati(1,i) = log(prior(i)) + log(B(i,O(1)));
    phiti(1,i) = 0;
end;

% Induction
for t=(2:N)
    for j=(1:Q)
        %sigmati(t,j) = max( sigmati(t-1,:).*(A(:,j).')*B(j,O(t)) );
        [sigmati(t,j),phiti(t,j)] = max( sigmati(t-1,:) + log(A(:,j).') );
        sigmati(t,j) = sigmati(t,j) + log(B(j,O(t)));
    end;
end;

% Termination
[logP,q(N)] = max(sigmati(end,:));

% Path (state-sequence) backtracking
for t=((N-1):(-1):1)
    q(t) = phiti(t+1,q(t+1));
end;

% Check output arguments
if (nargout<2)
    clear('P');
end;
