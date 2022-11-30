function [Iji,Iij] = vargranger(A,freqs,rate,C)
%function S = vargranger(A,freqs,rate,C)
%
%
%
%function S = vargranger(A,freqs,rate,C)

% Reserve memory
m=size(A,1);
p=size(A,3);
Iji=zeros(length(freqs),1);
Iij = Iji;

% Determine spectra per frequency
for n=(1:length(freqs))

    % Compute rational transfer function spectra
    H = eye(m);
    for k=(1:p)
        H = H + A(:,:,k)*exp(-1i*2*pi*freqs(n)/rate*k);
    end;
    S = inv(H)*C*(inv(H)');

    % Remove numerical discrepancy resulting in complex auto-spectra
    S = S - diag(1i*imag(diag(S)));
    
    Iji(n) = -log( 1 - (C(1,1)-(C(1,2)^2)/C(2,2))*(abs(H(2,1))^2)/S(2,2) );
    Iij(n) = -log( 1 - (C(2,2)-(C(2,1)^2)/C(1,1))*(abs(H(1,2))^2)/S(1,1) );

end;
