function [psi,psid,x] = psin( ph )
% Hyper-sphere PSI test
%
% An n-sphere corresponds to 1 radial length and (n-1) angles
%
% Takes a matrix [ N x M ] of phases. Phase differences are computed
% between all pairs of nodes ( m in M ) without duplicates.
%   N   time
%   M   nodes
%

% Convert phases to pairs of phase-differences
domega = zeros(size(ph,1),size(ph,2)*(size(ph,2)-1)/2);
kindex = 1;
for k1=(1:size(ph,2))
    for k2=((k1+1):size(ph,2))
        domega(:,kindex) = ph(:,k1)-ph(:,k2);
        kindex = kindex + 1;
    end
end
ph = domega;

% Number of observables
N = size(ph,1);
INdim = size(ph,2);
M = INdim + 1;              % M-sphere coordinates (input_dim - 1)

r = 1;
x = r*ones( N, M );
for n = (1:N)       % Each time point
    for m = (1:M)       % Each Cartesian coordinate
        for k = (1:m)       % Each product term
            if ( k == M )
                continue;   % No end cosine
            elseif ( k == m )
                % Last term for all ( except k = M )
                x(n,m) = x(n,m) * cos(ph(n,m));         % Cosine (k)
            else
                % Inner terms
                x(n,m) = x(n,m) * sin(ph(n,k));         % Sine (k)
            end
        end
    end
end

psid = mean( x, 1 );
psi = norm( psid );
