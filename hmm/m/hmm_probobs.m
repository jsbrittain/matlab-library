function [logp,alphati,betati] = hmm_probobs(O,prior,A,B,method)
%function [logp,[alphati,[betati]]] = hmm_probobs(O,prior,A,B,[method])
%
% Probability of the observation sequence given the model
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
%       method        Estimation method for alpha / beta parameters
%                       0: Direct (can lead to long sequence instability)
%                       n: Apply scaling function every (n) steps
%
% Output argument
%       log p         Log observation probability
%
% Reference
%   Rabiner, K. R. (1989) A Tutorial on Hidden Markov Models and Selected
%   Applications in Speech Recognition. Proceedings of the IEEE, 77(2):
%   257-286.
%
%function [logp,[alphati,[betati]]] = hmm_probobs(O,prior,A,B,[method])

% Determine input parameters
if  (~exist('method','var'))
    method = 1;
end;

% Determine data parameters
N=length(O);
Q=size(A,1);

% Forward algorithm

% Reserve memory
alphati = zeros(N,Q);
c = ones(N,1);

% Initialise
alphati(1,:) = (prior.').*B(:,O(1));
if (method>0)
    c(1) = 1/sum(alphati(1,:));
    alphati(1,:) = alphati(1,:)*c(1);
end;

% Induction
for t=(2:N)
    % Compute forward coefficients
    for j=(1:Q)
        %switch (method)
            %case 0,     % Directly (leads to long sequence instability without scaling)
                alphati(t,j) = sum(alphati(t-1,:).*(A(:,j).')) * B(j,O(t));
            %otherwise,  % Via a scaling function (Rabiner (1989), Eq.92a,b)
            %    alphati(t,j) = sum(alphati(t-1,:).*A(j,:).*(B(:,O(t)).'));
        %end;
    end;
    % Scale coefficients
    if (method>0)
        if (mod(t-1,method)==0)
            % Apply scaling every `method' (n) steps (otherwise default c(t)=1)
            c(t) = 1/sum(alphati(t,:));
        end;
        alphati(t,:) = alphati(t,:) * c(t);
    end;
end;

% Termination
logp = -sum(log(c));

% Backward algorithm
if (nargout>2)
    
    % Reserve memory
    betati = zeros(N,Q);
    d = zeros(N,1);
    
    % Initialise
    betati(N,:) = c(N);
%     if (method>0)
%         d(1) = 1/sum(betati(1,:));
%         betati(1,:) = betati(1,:)*d(1);
%     end;
    
    % Induction
    for t=((N-1):(-1):1)
        for i=(1:Q)
            % Compute backward coefficients and apply rescaling (can be unity)
            betati(t,i) = c(t) * sum( A(i,:).*(B(:,O(t+1)).').*betati(t+1,:) );
        end;
%         % Scale coefficients
%         if (method>0)
%             if (mod(t-1,method)==0)
%                 % Apply scaling every `method' (n) steps (otherwise default d(t)=1)
%                 d(t) = 1/sum(betati(t,:));
%             end;
%             betati(t,:) = betati(t,:) * d(t);
%         end;
    end;
end;

% Check output parameters
if (nargout<1)
    clear('alphati');
end;
