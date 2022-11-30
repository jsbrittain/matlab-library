function out=gmean(M,dim);
%function out=gmean(M,[dim]);
%
% Geometric mean
%
% Input parameters
%       M       Data matrix
%       dim     (opt) Take mean over dimension (default: 1)
%
%function out=gmean(M,[dim]);

if ((nargin<1) | (nargin>2))
    error(' Incorrect number of input arguments');
end;
if (~exist('dim'))
    dim=find(size(M)>1);        % Find non-singleton dimensions
    if (isempty('dim'))
        dim=1;
    else
        dim=dim(1);
    end;
end;

% Take geometric mean using an arithmetic mean in the log domain
if (isreal(M))
    out=exp(mean(log(M),dim));
    %out=prod(M,dim).^(1/size(M,dim));
else
    % Remove log(0) problem (for both real and imaginary parts)
    %M(real(M)<=exp(-100))=exp(-100)+i*imag(M(real(M)<=exp(-100)));
    %M(imag(M)<=exp(-100))=real(M(imag(M)<=exp(-100)))+i*exp(-100);
    % Determine geometric mean
    out=real(exp(mean(log(real(M)),dim)))+i*real(exp(mean(log(imag(M)),dim)));
end;
