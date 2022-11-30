function [sp,Psp]=kf_spm(sp,seg_count,Q);
%function [sp,Psp]=kf_spm(sp,seg_count,Q);
%
% Hybrid Kalman filter and Fourier based analysis
% Kalman filtering of Fourier-based coherence estimates over trials.
% Takes sp matrix as input.
%
% Reduced vector implementation returning spectral coefficients.
%
%function [sp,Psp]=kf_spm(sp,f_max,Q);

% Check sp columns
if (size(sp,2)>3)
    % Cross-spectra defined as seperate real and complex components
    sp(:,3,:)=sp(:,3,:)+i*sp(:,4,:);
    sp=sp(:,1:3,:);
end;

% Determine input parameters
M=size(sp,3);
f_max=size(sp,1);
if (nargin<3)                           % Determine Q
    Q=1e-3*ones(3*f_max,1);             % Default noise covariance
else
    if (size(Q,1)==1)
        Q=Q*ones(3*f_max,1);
    end;
end;

% Initialise variable space
x=zeros(3*f_max,M);
P=zeros(3*f_max,M);

% Kalman fitering parameters and initial conditions
x(:,1)=reshape(sp(1:f_max,:,1),3*f_max,1);         % Initial estimates
Rspect=ones(f_max,1)*(0.4343^2)/seg_count;
R=[Rspect; Rspect; Rspect];

P(:,1)=R;                                          % State error-covariance matrix
theta=ones(3*f_max,1);                             % State-transition matrix
H=ones(3*f_max,1);                                 % Process-to-measurement matrix

disp(['Running Kalman filter. Groups=' int2str(M)]);
xp=x(:,1);
Pp=P(:,1)+Q;
for ind=2:M
    % Construct state vector
    newz=reshape(sp(1:f_max,:,ind),3*f_max,1);
    
    % Perform Kalman filter
    [x(:,ind), P(:,ind), K, Pp, xp] = kf_filtr(newz, Q, R, theta, H, Pp, xp);
end;

% Formulate frequency-domain matrix output
sp=reshape(x(1:(f_max*3),:),f_max,3,M);
Psp=reshape(P(1:(f_max*3),:),f_max,3,M);
