function [x,P]=kf_coh(f,cl,f_max,Q);
%function [x,P]=kf_coh(f,cl,f_max);
%
% Hybrid Kalman filter and Fourier based analysis
% Kalman filtering of Fourier-based coherence estimates over trials.
%
% Arguments:
%   f_max   (opt) Maximum frequency 
%
%function [x,P]=kf_coh(f,cl,f_max);

% Determine input parameters
M=size(f,3);
if (nargin<3)
    f_max=size(f,1);
else
    f_max=find(f(:,1,1)>f_max);
    f_max=f_max(1)-1;
end;
if (nargin<4)
    Q=1e-5*eye(f_max);  % Default noise covariance
else
    if (size(Q,1)==1)
        Q=Q*eye(f_max);
    end;
end;

% Initialise variable space
x=zeros(f_max,M);
P=zeros(f_max,f_max,M);

% Kalman fitering parameters and initial conditions
xp=reshape(f(1:f_max,4,1),f_max,1);              % Initial estimates
R=diag((2/cl(1).seg_tot)*((abs(xp).^2).*((1-(abs(xp).^2)).^2)),0);   % PBMB (6.4)

Pp=R;                                            % Error covariance matrix
theta=eye(f_max);                                % State-transition matrix
H=eye(f_max);                                    % Process-to-measurement matrix

x(:,1)=xp;
P(:,:,1)=Pp;
for ind=2:M
    % Construct state vector
    newz=reshape(f(1:f_max,4,ind),f_max,1);
    R=diag((2/cl(ind).seg_tot)*((abs(newz).^2).*((1-(abs(newz).^2)).^2)),0);
    
    % Perform Kalman filter
    disp(['Running Kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    [x(:,ind), P(:,:,ind), K, Pp, xp] = kalm_filtr(newz, Q, R, theta, H, Pp, xp);
end;
