function [O,q] = hmm_generate(prior,A,B,N)
%function [O,q] = hmm_generate(prior,A,B,N)
%
% Input parameters
%       prior         Initial state distribution
%                       ( Q x 1 )
%       A             Transition probability matrix
%                       ( Q x Q )
%       B             Emission probability matrix
%                       ( Q x M )
%       N             Sequence length (samples)
%
% Output argument
%       O             Observation sequence
%
% Reference
%   Rabiner, K. R. (1989) A Tutorial on Hidden Markov Models and Selected
%   Applications in Speech Recognition. Proceedings of the IEEE, 77(2):
%   257-286.
%
%function [O,q] = hmm_generate(prior,A,B,N)

% Reserve memory
q=ones(1,N);
O=ones(1,N);

% Recursively estimate time-series
for t=(1:N)
    % Determine state
    if (t==1)
        [pq,inq]=sort(prior);
    else
        [pq,inq]=sort(A(q(t-1),:));
    end;
    q(t)=inq(find(rand(1)<cumsum(pq),1,'first'));
    % Determine observation
    [po,ino]=sort(B(q(t),:));
    O(t)=ino(find(rand(1)<cumsum(po),1,'first'));
end;
