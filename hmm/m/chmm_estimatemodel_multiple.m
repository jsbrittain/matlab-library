function [prior,A,gammati,epsilontij] = chmm_estimatemodel_multiple(b,prior,A,method,varargin)
%function [prior,A,[gammati,[epsilontij]]] = chmm_estimatemodel_multiple(b,prior,A,[method,[varargin]])
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
%       method        Probability scaling step (samples)
%                       0: no scaling, 1: default
%
% Output arguments
%       prior         Initial state distribution
%       A             Transition probability matrix
%
% References
%
%   Rabiner, K. R. (1989) A Tutorial on Hidden Markov Models and Selected
%   Applications in Speech Recognition. Proceedings of the IEEE, 77(2):
%   257-286.
%
%   Rahimi, A. An Erratum for ``A Tutorial on Hidden Markov Models and
%   Selected Applications in Speech Recognition''
%   http://alumni.media.mit.edu/~rahimi/rabiner/rabiner-errata/rabiner-errata.html
%
%function [prior,A,[gammati,[epsilontij]]] = chmm_estimatemodel_multiple(b,prior,A,[method,[varargin]])

% Determine input parameters
if  (~exist('method','var'))
    method = 1;
end;

% Specify default parameters
maxreps = 500;          % Maximum repetition count
th = 1e-6;

% Interpret varargin
if (mod(length(varargin),2)~=0)
    error(' Optional arguments must occur in parameter-value pairs.');
end;
optargs='';
for n=(1:2:length(varargin))
    switch (lower(varargin{n}))
        case 'threshold'
            th = varargin{n+1};
        otherwise
            optargs = [optargs ',''' varargin{n} ''',' num2str(varargin{n+1})];
    end;
end;

% Determine data parameters
K=length(b);
Q=size(A,1);

% Reserve memory space
epsilontij = cell(K,1);
gammati = cell(K,1);

% Compute probability of measurement being in each state
logp0 = zeros(1,K); alphati = cell(K,1); betati = cell(K,1); Nk = 0;
for k=(1:K)
    [logp0(k),alphati{k},betati{k}] = chmm_probobs(b{k},prior,A,method);
    Nk = Nk + size(alphati,1);
end;
logp0=sum(logp0);

% Iterate estimation procedure
logp = 0; rep = 0;
while (((1-logp/logp0)>th) && (rep < maxreps))
    
    % Expectation-Maximisation step (Baum-Welch step)
    for k=(1:K)
        [epsilontij{k},gammati{k}] = chmm_EM(b{k},A,alphati{k},betati{k});
    end;
    
    % Update prior (JSB)
    prior = zeros(size(prior));
    for k=(1:K)
        prior = prior + gammati{k}(1,:)/K;
    end;
    
    % Update transition matrix (Rahimi Eq.13)
    Anum = zeros(Q); Aden = zeros(Q);
    for i=(1:Q)
        for j=(1:Q)
            for k=(1:K)
                Anum(i,j) = Anum(i,j) + squeeze(sum(epsilontij{k}(:,i,j),1));
                Aden(i,j) = Aden(i,j) + sum(gammati{k}(:,i),1);
            end;
        end;
    end;
    A = Anum./Aden;
    
    % Update display
    logp0 = logp; logp = 0;
    for k=(1:K)
        [lgp,alphati{k},betati{k}] = chmm_probobs(b{k},prior,A,method);
        logp = logp+lgp;
    end;
    
    % Display progress
    disp([' Iteration ' num2str(rep) ' (log p=' num2str(logp/K/Nk) ')']);
    rep = rep + 1;
    
end;

% Check output arguments
if (nargout<4)
    clear('epsilontij');
    if (nargout<3)
        clear('gammati');
    end;
end;
