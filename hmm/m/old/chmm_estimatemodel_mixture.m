function [prior,A,gammati,epsilontij] = chmm_estimatemodel(O,prior,A,M,method)
%function [prior,A,[gammati,[epsilontij]]] =  = chmm_estimatemodel(O,prior,A,M,[method])
%
% Continuous observation (pdf) Hidden-Markov Model (HMM)
%
% Viterbi algorithm
%
% Returns the most likely set of states for a given model and observation
% set (dynamic programming approach; non-unique).
%
% Input parameters
%       O             Observation likelihood
%                       ( 1 x N )
%       prior         Initial state distribution
%                       ( Q x 1 )
%       A             Transition probability matrix
%                       ( Q x Q )
%       M             Mixing conditions (scalar)
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
%function [prior,A,[gammati,[epsilontij]]] =  = chmm_estimatemodel(b,prior,A,[method])

% Determine input parameters
if  (~exist('method','var'))
    method = 1;
end;

% Determine data parameters
N=length(O);
Q=size(A,1);

% Reserve memory space
epsilontij = zeros(N-1,Q,Q);
gammatik = zeros(N-1,Q,M);
b = zeros(N-1,Q);

% Fetch initial fwd-bwd coefficients
[~,alphati,betati] = chmm_probobs(b,prior,A,method);

% Iterate estimation procedure
reps=10;
for rep = (1:reps)
    
    % Calculate epsilon ( P(q_t=S_i, q_{t+1}=S_j | O, (prior, A, B)) )
    for t=(1:(N-1))
        for i=(1:Q)
            
            for k=(1:M)
                % Gamma parameter
                gammatik(t,i,k) = ((alphati(t,i)*betati(t,i))/sum(alphati(t,:).*betati(t,:)));
                gammatik(t,i,k) = gammatik(t,i,k) * ((c(i,k)*pdf('normal',O(t),mu(i,k),U(i,k)))/sum(c(i,:)*pdf('normal',O(t),mu(i,:),U(i,:))));
            end;
        end;
    end;
    % Estimate Gaussian mixing condition coefficients
    for i=(1:Q)
        for k=(1:M)
            % Mixing coefficients
            c(i,k) = sum(gammatik(:,i,k),1) / sum(sum(gammatik(:,i,:),1),3);
        end;
        % Mixing means
        mu(i,k) = sum(gammatik(:,i,k).*(O.')) / sum(gammatik(:,i,k));
        % Mixing variance
        U(i,k) = sum(gammatik(:,i,k).*((O.'-mu(i,k))*((O.'-mu(i,k)).'))) / sum(gammatik(:,i,k));
    end;
    
    for t=(1:(N-1))
        for i=(1:Q)
            for k=(1:M)
                b(t,i) = b(t,i) + c(i,k)*pdf('normal',O(t),mu(i,k),U(i,k));
            end;
            for j=(1:Q)
                epsilontij(t,i,j) = alphati(t,i)*A(i,j)*b(j,t+1)*betati(t+1,j);
            end;
        end;
        
        epsilontij(t,:,:) = epsilontij(t,:,:) / sum(sum(epsilontij(t,:,:)));
        %gammati(t,:) = sum(epsilontij(t,:,:),3);
    end;
    
    % Update parameters
    prior = gammati(1,:);
    for i=(1:Q)
        %sumgammati = sum(gammati(:,i));
        for j=(1:Q)
            A(i,j) = sum(epsilontij(:,i,j),1) / sum(sum(gammati(:,i,:),1),3);
        end;
    end;
    
    % Update display
    [logp,alphati,betati] = chmm_probobs(b,prior,A,method);
    disp([' Iteration ' num2str(rep) ' of ' num2str(reps) ' completed (log p=' num2str(logp) ').']);

end;

% Check output arguments
if (nargout<4)
    clear('epsilontij');
    if (nargout<3)
        clear('gammati');
    end;
end;
