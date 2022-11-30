function [prior,A,gammati,epsilontij] = chmm_estimatemodel(b,prior,A,method)
%function [prior,A,[gammati,[epsilontij]]] =  = chmm_estimatemodel(b,prior,A,[method])
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
%function [prior,A,[gammati,[epsilontij]]] =  = chmm_estimatemodel(b,prior,A,[method])

% Determine input parameters
if  (~exist('method','var'))
    method = 1;
end;

% Determine data parameters
N=size(b,2);
Q=size(A,1);

% Specify iteration threshold
th=1e-5;                % Percentage change cutoff

% Fetch initial fwd-bwd coefficients
[logp0,alphati,betati] = chmm_probobs(b,prior,A,method);

% Iterate estimation procedure
logp = 0; rep = 1;
while ((1-logp/logp0)>th)
    
    % Expectation-Maximisation step (Baum-Welch step)
    [epsilontij,gammati] = chmm_EM(b,A,alphati,betati);
    
    % Update parameters
    prior = gammati(1,:);
    for i=(1:Q)
        for j=(1:Q)
            A(i,j) = sum(epsilontij(:,i,j),1) / sum(gammati(:,i),1);
        end;
    end;
    
    % Update display
    logp0=logp;
    [logp,alphati,betati] = chmm_probobs(b,prior,A,method);
    disp(['  Iteration ' num2str(rep) ' (log p=' num2str(logp) ', ' num2str((1-logp/logp0)/(100*th)) ')']);
    rep = rep + 1;
    
end;

% Check output arguments
if (nargout<4)
    clear('epsilontij');
    if (nargout<3)
        clear('gammati');
    end;
end;
