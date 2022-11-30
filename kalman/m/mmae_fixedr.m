function [condx, condP, Pp, xp, p]=mmae_fixedr(newz, Q, R, theta, H, Pp, xp, p, channels, filters);
%function [condx, condP, Pp, xp, p]=mmae_fixedr(newz, Q, R, theta, H, Pp, xp, p, channels, filters);
%
% Multi-Modal Adaptive Estimator (MMAE)
% Implementation: Fixed Model ; Multivariate ; Recursive
%
% Outputs
%   condx       Conditioned estimate of state x.
%   condP       Conditioned error covariance matrix.
%   Pp          'a posteriori' error covariance matrix for all filters of the
%               filter bank.
%   xp          'a posteriori' state estimate for all filters of the filter
%               bank.
%   p           Filter bank probabilities
%
% Parameters
%   newz        New measurement samples (for time t).
%   Q           Process noise covariance.
%   R           Measurement noise covariance.
%   theta       State transition matrix.
%   H           State-to-measurement transition matrix.
%   Pp          'a priori' error covariance matrix.
%   xp          'a priori' state estimates.
%   p           Model probability
%   channels    Number data channels
%   filters     No. of filter to implement in the MMAE algorithm.
%
% function [condx, condP, Pp, xp, p]=mmae_fixedr(newz, Q, R, theta, H, Pp, xp, p, channels, filters);

% Allocate memory
x=zeros(channels,filters);
P=zeros(channels,channels,filters);
f=zeros(channels,channels,filters);

% Calculate conditional probability (Welch MMAE (a))
for j=1:filters
	C=H*Pp(:,:,j)*H' + R(:,:,j);
	f(:,:,j)=mmae_dens(newz, H*xp(:,j), C);
end;

% Calculate filter probability (Welch MMAE (b))
fp_sum=0;
for j=1:filters
    fp_sum=fp_sum + f(:,:,j)*p(:,:,j);
end;
for j=1:filters
    p(:,:,j)=f(:,:,j)*p(:,:,j)/fp_sum;
end;

% Model conditioned estimate (Welch MMAE (c)) - includes 'a priori' estimates for the next time loop
for j=1:filters
    [x(:,j), P(:,:,j), K(:,:,j), Pp(:,:,j), xp(:,j)] = kalm_filtr(newz, Q(:,:,j), R(:,:,j), theta, H, Pp(:,:,j), xp(:,j));
end;

% Compute model-conditioned state estimate (Welch MMAE (d))
condx=0;
for j=1:filters
    condx = condx + p(:,:,j)*x(:,j);
end;

% Compute model-conditioned error covariance matrix(Welch MMAE (e))
condP=0;
for j=1:filters
    e=condx - x(:,j);
    condP = condP + p(:,:,j)*(P(:,:,j) + e*e');
end;