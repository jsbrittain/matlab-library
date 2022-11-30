function [xN, PN] = minimax_fb(z, theta, H, gamma, U, V, W);
%function [xN, PN] = minimax_fb(z, theta, H, gamma, U, V, W);
%
%  Multivariate Minimax Forward Backward Filter
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
%function [xN, PN] = minimax_fb(z, theta, H, gamma, U, V, W);
%

channels=size(z,1);
len=size(z,2);

% Setup structures for faster computation
xN=zeros(channels, len);
PN=zeros(channels, channels, len);

% Initial conditions
xp(:,1)=z(:,1);
Pp(:,:,1)=zeros(channels, channels);
I=eye(channels);

% Perform forward and backward filter sweeps
[x1, P1, K1, L1] = minimax_filt(z, theta, H, gamma, U, V, W);
[x2, P2, K2, L2] = minimax_filt(fliplr(z), theta, H, gamma, U, V, W);

% Combine forward and backward filters to produce a smoothed estimate
for t=1:len
    PN(:,:,t)=inv( inv(P1(:,:,t)) + inv(P2(:,:,len-t+1)) );
    xN(:,t)=PN(:,:,t)*[inv(P1(:,:,t))*x1(:,t) + inv(P2(:,:,len-t+1))*x2(:,len-t+1)];
end;
