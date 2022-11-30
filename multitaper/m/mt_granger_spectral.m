function I = mt_granger_spectral( S )

% Determine parameters
M = size(S,1);              % Channel count

% Reverse memory
H = zeros(size(S));
Sigma = H;

% Spectral factorisation
for fn = (1:size(S,3))
    % Spectral factorisation
    switch (1)
        case 0,         % Schur decomposition
            [H(:,:,fn),Sigma(:,:,fn)]=schur(S(:,:,fn),'complex');
        case 1,         % SVD
            % Assumes first and third outputs are Hermitian symmetrical (tend to be)
            [H(:,:,fn),Sigma(:,:,fn),~]=svd(S(:,:,fn));
    end;
end;

% FFT-derived Granger causality
I = zeros(size(S));
for i1 = (1:M)
    for i2 = (1:M)
        if (i1==i2), continue; end;
        I(i2,i1,:) = log(S(i1,i1,:)./(S(i1,i1,:)-(Sigma(i2,i2,:)-(Sigma(i1,i2,:).^2)./Sigma(i1,i1,:)).*(abs(H(i1,i2,:)).^2)));
    end;
end;
