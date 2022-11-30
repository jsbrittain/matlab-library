function S = wvarspectrum(A,freqs,rate,S2,lambda)
%function S = wvarspectrum(A,freqs,rate,S2,lambda)
%
%
%
%function S = wvarspectrum(A,freqs,rate,S2,lambda)

% Reserve memory
m=size(A,1);
p=size(A,3);
S=zeros(m,m,length(freqs),size(A,4));

if (size(S2,3)==1)
    S2 = repmat(S2,[1 1 size(A,4)]);
end;

% Warped frequency vector
wf = inline(['atan2((1-' num2str(lambda) '^2)*sin(x),((1+' num2str(lambda) '^2)*cos(x)-2*' num2str(lambda) '));']);
wfreqs = wf(freqs*2*pi/rate);

% Recurse time-points
for t=(1:size(A,4))
    % Update displays
    disp(['Computing time point ' num2str(t) '/' num2str(size(A,4))]);
    
    % Determine spectra per frequency
    for n=(1:length(wfreqs))

        % Compute rational transfer function spectra
        H = eye(m);
        for k=(1:p)
            H = H + A(:,:,k,t)*exp(-1i*2*pi*wfreqs(n)/rate*k);
        end;
        S(:,:,n,t) = inv(H)*S2(:,:,t)*(inv(H)');

        % Remove numerical discrepancy resulting in complex auto-spectra
        S(:,:,n,t) = S(:,:,n,t) - diag(1i*imag(diag(S(:,:,n,t))));

    end;

end;

% Scale by frequency sampling
S = S .* repmat( permute([0 diff(wfreqs)],[1 3 2]), [ size(S,1) size(S,2) 1 size(S,4) ] );
