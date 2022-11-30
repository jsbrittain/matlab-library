function [x, P, K, Pp, xp] = kf_filtr(newz, Q, R, theta, H, Pp, xp);
%function [x, P, K, Pp, xp] = kf_filtr(newz, Q, R, theta, H, Pp, xp);
%
%  Multivariate Kalman Filter
%  Recursive form
%
%  Random process to be estimated
%    x(k+1) = theta(k)*x(k) + w(k)        // Process
%    z(k)   = H(k)*x(k) + v(k)            // Measurement
%
%  Function parameters
%    newz     observations (measurements)
%    Q        process noise covariance matrix
%    R        measurement noise covariance matrix
%    theta    state transition matrix
%    H        measurement-to-state matrix
%    Pp       'a priori' error covariance matrix
%    xp       'a priori' estimate sequence
%
%function [x, P, K, Pp, xp] = kf_filtr(newz, Q, R, theta, H, Pp, xp);
%

% Initial conditions
channels=size(newz,1);
I=ones(channels,1);

% Kalman gain
K=Pp.*H ./ ( H.*Pp.*H + R );

% Update estimate with measurement
x=xp + K.*(newz - H .* xp);

% Compute error covariance for updated estimate
P=(I-K.*H).*Pp;

% Project ahead
xp=theta.*x;
Pp=theta.*P.*theta + Q;
