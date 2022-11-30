function [f,P]=kf_spm(f,t,cl,f_max,Q);
%function [f,P]=kf_spm(f,t,cl,f_max,Q);
%
% Hybrid Kalman filter and Fourier based analysis
% Kalman filtering of Fourier-based coherence estimates over trials.
%
%function [f,P]=kf_spm(f,t,cl,f_max,Q);

% Determine input parameters
M=size(f,3);
if (nargin<4)
    f_max=size(f,1);
else
    f_max=find(f(:,1,1)>f_max);
    f_max=f_max(1)-1;
end;
if (nargin<5)
    Q=1e-5*ones(4*f_max,1);  % Default noise covariance
else
    if (size(Q,1)==1)
        Q=Q*ones(4*f_max,1);
    end;
end;

% Initialise variable space
x=zeros(4*f_max,M);
P1=zeros(4*f_max,M);

% Kalman fitering parameters and initial conditions
coh=reshape(f(1:f_max,4,1),f_max,1);
xp=reshape(f(1:f_max,2:5,1),4*f_max,1);          % Initial estimates
Rspect=ones(f_max,1)*(0.4343^2)/cl(1).seg_tot;
Rcoh=(2/cl(1).seg_tot)*((abs(coh).^2).*((1-(abs(coh).^2)).^2));   % PBMB (6.4)
Rphase=(1/(2*cl(1).seg_tot))*((1./(abs(coh).^2))-1);
R=[Rspect; Rspect; Rcoh; Rphase];

Pp=R;                                              % State error-covariance matrix
theta=ones(4*f_max,1);                             % State-transition matrix
H=ones(4*f_max,1);                                 % Process-to-measurement matrix

x(:,1)=xp;
P1(:,1)=Pp;
for ind=2:M
    % Construct state vector
    coh=reshape(f(1:f_max,4,ind),f_max,1);
    newz=reshape(f(1:f_max,2:5,ind),4*f_max,1);
    
    % Construct measurement noise covariance matrix
    Rcoh=(2/cl(ind).seg_tot)*((abs(coh).^2).*((1-(abs(coh).^2)).^2));
    Rphase=(1/(2*cl(ind).seg_tot))*((1./(abs(coh).^2))-1);
    R=[Rspect; Rspect; Rcoh; Rphase];
    
    % Perform Kalman filter
    disp(['Running Kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    [x(:,ind), P1(:,ind), K, Pp, xp] = kf_filtr(newz, Q, R, theta, H, Pp, xp);
end;

% Formulate output
f=f(1:f_max,1,:);
f(:,2:5,:)=reshape(x,f_max,4,M);
P=zeros(f_max,1,M);
P(:,2:5,:)=reshape(P1,f_max,4,M);
