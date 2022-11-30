function X = nanpad( varargin )
%
% Construct a ( Nmax x M ) matrix from M vector inputs of length Nm each
% with each column padded with nan's.
%
% This function expects column vector inputs and performs no error
% checking.
%

% Get maximum vector length
Nmax = 0;
M = length(varargin);
for m = (1:M)
    Nmax = max( Nmax, length(varargin{m}) );
end;

% Construct NaN matrix
X = nan( Nmax, M );

% Fill each column with its corresponding elements
for m = (1:M)
    X((1:length(varargin{m})),m) = varargin{m};
end;
