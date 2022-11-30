function [x, P, K, L] = minimax_filt(z, theta, H, gamma, U, V, W);
%function [x, P, K, L] = minimax_filt(z, theta, H, gamma, U, V, W);
%
%  Multivariate Minimax filter
%
%  Random process to be estimated
%    x(k+1) = theta(k)*x(k) + w(k)        // Process
%    z(k)   = H(k)*x(k) + v(k)            // Measurement
%
%  Function parameters
%    z             Observations (measurements)
%    Q             Process noise covariance matrix (associated with w(k))
%                   (CONSTANT)
%    R             Measurement noise covariance matrix (associated with v(k))
%                   (CONSTANT)
%    H             Measurement-to-state matrix (Optional, default: I)
%                   (CONSTANT)
%    theta         State transition matrix (Optinal, default: I)
%                   (CONSTANT)
%
%function [x, P, K, L] = minimax_filt(z, theta, H, gamma, U, V, W);

len=size(z,2);
channels=size(z,1);

% Setup structures for faster computation
x=zeros(channels, len);
L=zeros(channels, channels, len);
K=zeros(channels, channels, len);
P=zeros(channels, channels, len);

% Initial conditions
x(:,1)=z(:,1);
P(:,:,1)=zeros(channels, channels);
I=eye(channels);

% Set system parameters
if nargin<5
    if nargin<4
        theta=eye(channels);  % State transition matrix
    end;
    H=eye(channels);          % Measurment-to-state matrix
end;

% Perform prediction-correction loop
for k=1:len

    % Pre-gain
    L(:,:,k)=inv( I - gamma*U*P(:,:,k) + H'*inv(V)*H*P(:,:,k) );
    
	% Minimax gain
	K(:,:,k)=theta*P(:,:,k)*L(:,:,k)*H'*inv(V);
	
    % State estimator
    x(:,k+1)=theta*x(:,k) + K(:,:,k)*(z(:,k) - H*x(:,k));
    
    % Error covariance matrix
    P(:,:,k+1)=theta*P(:,:,k)*L(:,:,k)*theta' + W;
    
    % Ensure eigenvalues of P are < 1 in magnitude
    lambda=eig(P(:,:,k+1));
    if (sum(abs(lambda) >= 1) > 0)
        disp('gamma is too large');
        return;
    end;
end;
