function D = hmm_dist(lambda1,lambda2,N)
%
%
%
%
%
%
%
%
%

% Observation sequence
if (~exist('N','var'))
    N = [];
end;

% Extract HMM model 1
prior1 = lambda1{1};
A1 = lambda1{2};
B1 = lambda1{3};
% Extract HMM model 2
prior2 = lambda2{1};
A2 = lambda2{2};
B2 = lambda2{3};

% Extract parameters
Q = size(A1,1);

% Set observation matrix to I should it be empty
if (isempty(B1))
    B1 = eye(Q);
end;
if (isempty(B2))
    B2 = eye(Q);
end;
if (isempty(N))
    N = 1000;
end;

% Generate observation sequence
[O2,dummy] = hmm_generate(prior2,A2,B2,N);

% P( O_2 | lambda_1 )
[dummy,logP1] = hmm_viterbi(O2,prior1,A1,B1);
% P( O_2 | lambda_2 )
[dummy,logP2] = hmm_viterbi(O2,prior2,A2,B2);

% Distance measure
D = (logP1 - logP2) / N;
