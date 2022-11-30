function S = varspectrum(A,freqs,rate,S2)
%function S = varspectrum(A,freqs,rate,S2)
%
%
%
%function S = varspectrum(A,freqs,rate,S2)

% Reserve memory
m=size(A,1);
p=size(A,3);
S=zeros(m,m,length(freqs),size(A,4));

% Recurse time-points
for t=(1:size(A,4))
    % Update displays
    disp(['Computing time point ' num2str(t) '/' num2str(size(A,4))]);
    
    % Determine spectra per frequency
    for n=(1:length(freqs))

        % Compute rational transfer function spectra
        H = eye(m);
        for k=(1:p)
            H = H + A(:,:,k,t)*exp(-1i*2*pi*freqs(n)/rate*k);
        end;
        S(:,:,n,t) = inv(H)*S2*(inv(H)');

        % Remove numerical discrepancy resulting in complex auto-spectra
        S(:,:,n,t) = S(:,:,n,t) - diag(1i*imag(diag(S(:,:,n,t))));

    end;

end;
