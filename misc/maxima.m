function out = maxima(x)
%function out = maxima(x)
%
% Find local maxima in the data
%
% Function written by David Sampson, available from
%   http://www-h.eng.cam.ac.uk/help/tpl/programs/Matlab/minmax.html
%

% Unwrap to vector
x = x(:);

% Identify whether signal is rising or falling
upordown = sign(diff(x));

% Find points where signal is rising before, falling after
maxflags = [upordown(1)<0; diff(upordown)<0; upordown(end)>0];
out      = find(maxflags);  % Maxima
