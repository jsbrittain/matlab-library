function [epsilontij,gammati] = chmm_EM(b,A,alphati,betati)
%
% Expectation-Maximisation step
%
%
% References
%   Rahimi, A. An Erratum for ``A Tutorial on Hidden Markov Models and
%   Selected Applications in Speech Recognition''
%   http://alumni.media.mit.edu/~rahimi/rabiner/rabiner-errata/rabiner-errata.html
%

% Determine data parameters
N=size(b,2);
Q=size(A,1);

% Reserve variable space
epsilontij = zeros(N-1,Q,Q);
gammati = zeros(N-1,Q);

% Calculate epsilon ( P(q_t=S_i, q_{t+1}=S_j | O, (prior, A)) )
for t=(1:(N-1))
    for i=1:Q
        for j=1:Q
            % See Rabiner Erratum (by Rahimi)
            epsilontij(t,i,j) = alphati(t,i)*A(i,j)*b(j,t+1)*betati(t+1,j);
        end;
    end;
    gammati(t,:) = sum(epsilontij(t,:,:),3);
end;
