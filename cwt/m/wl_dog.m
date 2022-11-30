function psi=wl_dog(nu,m);
%function psi=wl_dog(nu,m);
%
% Wavelet: Derivative of Gaussian (DOG)
% Time-domain
%
% Common parameters
%   nu      Time
%
% Additional parameters
%   m       Derivative order
%

% Normalisation
norm=((-1)^(m+1))/sqrt(gamma(m+.5));

% Determine m'th derivative of gaussian
switch m
    case 0
        coeffs=  1;
    case 1
        coeffs= -nu;
    case 2
        coeffs= (nu.^2)-1;
    case 3
        coeffs=-(nu.^3)+3*nu;
    case 4
        coeffs= (nu.^4)-6*(nu.^2)+3;
    case 5
        coeffs=-(nu.^5)+10*(nu.^3)-15*nu;
    case 6
        coeffs= (nu.^6)-15*(nu.^4)+45*(nu.^2)-15;
    otherwise
        error([int2str(m) 'th derivative of Gaussian not yet supported.']);
end;
deriv=coeffs.*exp(-(nu.^2)/2);

% Normalised time-domain wavelet function
psi=norm*deriv;
