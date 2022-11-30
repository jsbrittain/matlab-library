function [x, P, K, Pp, xp] = kalm_filt(z, Q, R, theta, H)
%function [x, P, K, Pp, xp] = kalm_filt(z, Q, R, theta, H)
%
%  Multivariate Kalman filter
%    Specialised to same state/measurement order matrices
%
%  Random process to be estimated
%    x(k+1) = theta(k)*x(k) + w(k)        // Process
%    z(k)   = H(k)*x(k) + v(k)            // Measurement
%
%  Function parameters
%    z             Observations (measurements)
%    Q             Process noise covariance matrix (associated with w(k))
%    R             Measurement noise covariance matrix (associated with v(k))
%    H             Measurement-to-state matrix (Optional, default: I)
%    theta         State transition matrix (Optinal, default: I)
%
%function [x, P, K, Pp, xp] = kalm_filt(z, Q, R, theta, H)

len=size(z,2);
channels=size(z,1);

if ( channels > len )
    warning('Expecting time-series to be rows - is this data transpose?');
end

% Setup structures for faster computation
x=zeros(channels, len);
xp=zeros(channels, len+1);
P=zeros(channels, channels, len);
Pp=zeros(channels, channels, len+1);
K=zeros(channels, channels, len);


% Initial conditions
xp(:,1)=z(:,1);
Pp(:,:,1)=eye(channels, channels);
I=eye(channels);

% Set system parameters
if nargin<5
    if nargin<4
        theta=eye(channels);  % State transition matrix
    end
    H=eye(channels);          % Measurment-to-state matrix
end

% Perform prediction-correction loop
for k=1:len
    % Kalman gain
    K(:,:,k)=Pp(:,:,k)*H' / (H*Pp(:,:,k)*H'+R);
    
    % Update estimate with measurement
    x(:,k)=xp(:,k) + K(:,:,k)*(z(:,k) - H * xp(:,k));
    
    % Compute error covariance for updated estimate
    P(:,:,k)=(I-K(:,:,k)*H)*Pp(:,:,k);
    
    % Project ahead
    xp(:,k+1)=theta*x(:,k);
    Pp(:,:,k+1)=theta*P(:,:,k)*theta' + Q;
end
