function [prior,A,B,gammati,epsilontij] = hmm_estimatemodel(O,prior,A,B,method,precision)
%function [prior,A,B,[gammati,[epsilontij]]] =  = hmm_estimatemodel(O,prior,A,B,[method])
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
% Output arguments
%       prior         Initial state distribution
%       A             Transition probability matrix
%       B             Emission probability matrix
%
% Reference
%   Rabiner, K. R. (1989) A Tutorial on Hidden Markov Models and Selected
%   Applications in Speech Recognition. Proceedings of the IEEE, 77(2):
%   257-286.
%
%function [prior,A,B,[gammati,[epsilontij]]] =  = hmm_estimatemodel(O,prior,A,B,[method])

% Determine input parameters
if  (~exist('method','var'))
    method = [];
end;
if (~exist('precision','var'))
    precision = [];
end;

% Default parameters
if (isempty(method))
    method = 1;
end;
if (isempty(precision))
    precision = 0.001;
end;

% Determine data parameters
N=length(O);
Q=size(A,1);
M=size(B,2);

% Reserve memory space
epsilontij = zeros(N-1,Q,Q);
gammati = zeros(N-1,Q);

% Fetch initial fwd-bwd coefficients
[~,alphati,betati] = hmm_probobs(O,prior,A,B,method);

% Iterate estimation procedure
logp0 = 1; logp = 10; rep = 0; th = precision/100;
while ( (abs(logp-logp0)/abs(logp0))>th )     % Iterate until change is less than the given %age
    logp0 = logp;
    rep = rep + 1;
    if (rep>200)
        warning('Over 200 iterations... returning aborted HMM fit.');
        break;
    end;
    
    % Calculate epsilon ( P(q_t=S_i, q_{t+1}=S_j | O, (prior, A, B)) )
    for t=(1:(N-1))
        for i=1:Q
            for j=1:Q
                epsilontij(t,i,j) = alphati(t,i)*A(i,j)*B(j,O(t+1))*betati(t+1,j);
            end;
        end;
        
        epsilontij(t,:,:) = epsilontij(t,:,:) / sum(sum(epsilontij(t,:,:)));
        gammati(t,:) = sum(epsilontij(t,:,:),3);
    end;
    
    % Update parameters
    prior = gammati(1,:);
    for i=(1:Q)
        sumgammati = sum(gammati(:,i));
        for j=(1:Q)
            A(i,j) = sum(epsilontij(:,i,j),1) / sumgammati;
        end;
        for k=(1:M)
            B(i,k) = sum(gammati(O(1:(end-1))==k,i)) / sumgammati;
        end;
    end;
    
    % Update display
    [logp,alphati,betati] = hmm_probobs(O,prior,A,B,method);
    disp([' Iteration ' num2str(rep) ' (log p=' num2str(logp) ', change = ' num2str(100*abs(logp-logp0)/abs(logp0)) ' % [threshold ' num2str(precision) ' %]).']);

end;

% Check output arguments
if (nargout<5)
    clear('epsilontij');
    if (nargout<4)
        clear('gammati');
    end;
end;
