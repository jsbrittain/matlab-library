function [condx, condP, Pp, xp, p, condxp, condPp]=mmae_dynamicr(newz, Q, R, theta, H, Pp, xp, channels, filters);
%function [condx, condP, Pp, xp, p, condxp, condPp]=mmae_dynamicr(newz, Q, R, theta, H, Pp, xp, channels, filters);
%
% Multi-Modal Adaptive Estimator (MMAE)
% Implementation: Dynamic Model ; Multivariate ; Recursive
%
% Outputs
%   condx       Conditioned estimate of state x.
%   condP       Conditioned error covariance matrix.
%   Pp          'a posteriori' error covariance matrix for all filters of the
%               filter bank.
%   xp          'a posteriori' state estimate for all filters of the filter
%               bank.
%   p           Filter bank probabilities
%   condxp      Conditioned estimate of 'a priori' x for (t+1) (Optional)
%   condPp      Conditioned 'a priori' error cov. matrix for (t+1) (Optional)
%
% Parameters
%   newz        New measurement samples (for time t).
%   Q           Process noise covariance.
%   R           Measurement noise covariance.
%   theta       State transition matrix.
%   H           State-to-measurement transition matrix.
%   Pp          'a priori' error covariance matrix.
%   xp          'a priori' state estimates.
%   channels    Number data channels
%   filters     No. of filter to implement in the MMAE algorithm.
%
% function [condx, condP, Pp, xp, p, condxp, condPp]=mmae_dynamicr(newz, Q, R, theta, H, Pp, xp, channels, filters);

% Allocate memory
x=zeros(channels,filters);
P=zeros(channels,channels,filters);
f=zeros(channels,channels,filters);

% Calculate conditional probability (Welch MMAE (a))
for chan=1:channels
	for j=1:filters
		C=H(chan,chan)*Pp(chan,chan,j)*H(chan,chan)' + R(chan,chan,j);
		f(chan,chan,j)=mmae_dens(newz(chan), H(chan,chan)*xp(chan,j), C);
	end;
end;

% Calculate filter probability (Welch MMAE (b))
f_sum=sum(f,3);
for j=1:filters
    p(:,:,j)=f(:,:,j)/f_sum;
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

% Produce conditioned 'a priori' estimates
if nargout>5
    condxp=0;
    condPp=0;
	for j=1:filters
        condxp = condxp + p(:,:,j)*xp(:,j);
        condPp = condPp + p(:,:,j)*Pp(:,:,j);
	end;
end;
