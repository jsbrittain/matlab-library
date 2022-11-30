function [h,y,v1] = plotn( X, varargin )
% Plot n-dimensional data projected into 2D (i.e. a monitor!)

N = size(X,1);
M = size(X,2);

% Setup 2-dimensional unit-vector axis planes
switch ( 1 )
    case 1,     % Evenly distributed around the unit circle
        dph = 2*pi/M;
        ph = (0:(M-1))*dph + pi/2;  % (-pi) ensures one axis points up!
        v1 = exp(1i*ph);            % Unit vectors in 2d using complex notation
end;

% Transform n-dimensional matrix into a 2-dimensional vector
y = sum( X.*repmat(v1,N,1), 2 );        % Scale unit vectors and sum

% Plot in 2D
h = plot( real(y), imag(y), varargin{:} );

% Plot axes
if ( true )
    hold('on');
    plot( [zeros(1,length(v1)); real(v1)], [zeros(1,length(v1)); imag(v1)], 'k' );
    axis image; axis off
end;
