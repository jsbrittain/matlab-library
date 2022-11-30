function [x, P, K, Pp, xp] = ekf_filtr(newz, Q, R, f, h, Pp, xp);
%function [x, P, K, Pp, xp] = ekf_filtr(newz, Q, R, f, h, Pp, xp);
%
%  Extended Kalman Filter
%  Recursive form
%
%  Random process to be estimated
%    x(k) = f( x(k-1), w(k-1) )       // Process (No control inputs)
%    z(k) = h( x(k), v(k) )           // Measurement
%
%  Function parameters
%    newz     observations (measurements)
%    Q        process noise covariance matrix
%    R        measurement noise covariance matrix
%    f        state transition function (handle)
%    h        measurement-to-state function (handle)
%    Pp       'a priori' error covariance matrix
%    xp       'a priori' estimate sequence
%
%function [x, P, K, Pp, xp] = ekf_filtr(newz, Q, R, f, h, Pp, xp);
%

% Determine function handles
if (~isa(f,'function_handle'))
    error('EKF: f must be a function handle');
end;
if (~isa(h,'function_handle'))
    error('EKF: h must be a function handle');
end;

% Initial conditions
channels=size(xp,1);
I=eye(channels);
bigO=zeros(size(xp));

% Construct Jacobian for measurement update (correction-step)
H=jacobian(h,{xp,bigO},1);         % part deriv h wrt x
V=jacobian(h,{xp,big0},2);         % part deriv h wrt v

% Kalman gain
K=Pp*H' * inv( H*Pp*H' + V*R*V' );

% Update estimate with measurement
x=xp + K*(newz - feval(h,xp,0));

% Compute error covariance for updated estimate
P=(I-K*H)*Pp;

% Compute Jacobian for time update (prediction-step)
A=jacobian(f,{x,big0},1);      % part deriv f wrt x
W=jacobian(f,{x,big0},2);      % part deriv f wrt w

% Project ahead
xp=feval(f,x,0);
Pp=A*P*A' + W*Q*W';
