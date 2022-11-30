function [xN, PN] = kalm_fb(z, Q, R, theta, H);
%function [xN, PN] = kalm_fb(z, Q, R, theta, H);
%
%  Multivariate Kalman Forward Backward Smoother
%  (Batch processing version)
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
%function [xN, PN] = kalm_fb(z, Q, R, theta, H);

len=size(z,2);
channels=size(z,1);

if nargin<4
    theta=eye(channels);
    H=eye(channels);
end;

% Setup structures for faster computation
xN=zeros(channels, len);
PN=zeros(channels, channels, len);

% Initial conditions
xp(:,1)=z(:,1);
Pp(:,:,1)=zeros(channels, channels);
I=eye(channels);

% Perform forward and backward filter sweeps
[x1, P1, K1, Pp1, xp1]=kalm_filt(z, Q, R, theta, H);
[x2, P2, K2, Pp2, xp2]=kalm_filt(fliplr(z), Q, R, theta, H);

% Combine forward and backward filters to produce a smoothed estimate
for t=1:len
    PN(:,:,t)=inv( inv(P1(:,:,t)) + inv(Pp2(:,:,len-t+1)) );
    xN(:,t)=PN(:,:,t)*[inv(P1(:,:,t))*x1(:,t) + inv(Pp2(:,:,len-t+1))*xp2(:,len-t+1)];
end;