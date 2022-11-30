function [f,Pf,t,Pt]=kf_spm(f,t,cl,f_max,Q);
%function [f,Pf,t,Pt]=kf_spm(f,t,cl,f_max,Q);
%
% Hybrid Kalman filter and Fourier based analysis
% Kalman filtering of Fourier-based coherence estimates over trials.
%
% Reduced vector implementation returning frequency-domain parameters,
% associated error covariance parameters, and optionally the set of
% time-domain parameters and associated error-covariance parameters.
%
%function [f,Pf,t,Pt]=kf_spm(f,t,cl,f_max,Q);

% Check parameter validity
if (nargin<3)
    error('You must provide at least the first three input arguments.');
end;
if ((nargout~=2) & (nargout~=4))
    error('Number of output arguments must be either 2 or 4.');
end;
if (nargout==4)
    tout=1;     % Assess time-domain parameters
else
    tout=0;
end;

% Determine input parameters
M=size(f,3);
t_max=size(t,1);
if ((nargin<4) | (f_max==0))        % Determine f_max
    f_max=size(f,1);
else
    f_first=find(f(:,1,1)>f_max);
    if (isempty(f_first))
        f_max=size(f,1);
    else
        f_max=f_first(1)-1;
    end;
    if (f_max==0)
        error('Specified f_max too small.');
    end;
end;
if (nargin<5)               % Determine Q
    Q=1e-5*ones(4*f_max+tout*t_max,1);  % Default noise covariance
else
    if (size(Q,1)==1)
        Q=Q*ones(4*f_max+tout*t_max,1);
    end;
end;

% Initialise variable space
x=zeros(4*f_max+tout*t_max,M);
P=zeros(4*f_max+tout*t_max,M);

% Kalman fitering parameters and initial conditions
coh=reshape(f(1:f_max,4,1),f_max,1);
xp=reshape(f(1:f_max,2:5,1),4*f_max,1);          % Initial estimates
Rspect=ones(f_max,1)*(0.4343^2)/cl(1).seg_tot;
Rcoh=(2/cl(1).seg_tot)*((abs(coh).^2).*((1-(abs(coh).^2)).^2));   % PBMB (6.4)
Rphase=(1/(2*cl(1).seg_tot))*((1./(abs(coh).^2))-1);
if (tout)
    xp=[xp; t(:,2,1)];
    Rcum=ones(t_max,1)*((cl(1).q_c95/1.96)^2);
else
    Rcum=[];
end;
R=[Rspect; Rspect; Rcoh; Rphase; Rcum];

Pp=R;                                              % State error-covariance matrix
theta=ones(4*f_max+tout*t_max,1);                  % State-transition matrix
H=ones(4*f_max+tout*t_max,1);                      % Process-to-measurement matrix

x(:,1)=xp;
P(:,1)=Pp;
for ind=2:M
    % Construct state vector
    coh=reshape(f(1:f_max,4,ind),f_max,1);
    newz=reshape(f(1:f_max,2:5,ind),4*f_max,1);
    
    % Construct measurement noise covariance matrix
    Rcoh=(2/cl(ind).seg_tot)*((abs(coh).^2).*((1-(abs(coh).^2)).^2));
    Rphase=(1/(2*cl(ind).seg_tot))*((1./(abs(coh).^2))-1);
    if (tout)
        newz=[newz; t(:,2,ind)];
        Rcum=ones(t_max,1)*((cl(ind).q_c95/1.96)^2);
    end;
    R=[Rspect; Rspect; Rcoh; Rphase; Rcum];
    
    % Perform Kalman filter
    disp(['Running Kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    [x(:,ind), P(:,ind), K, Pp, xp] = kf_filtr(newz, Q, R, theta, H, Pp, xp);
end;

% Formulate frequency-domain matrix output
f=f(1:f_max,1,:);
f(:,2:5,:)=reshape(x(1:(f_max*4),:),f_max,4,M);
Pf=zeros(f_max,1,M);
Pf(:,2:5,:)=reshape(P(1:(f_max*4),:),f_max,4,M);

% Formulate time-domain matrix output
if (tout)
    t=t(1:t_max,1,:);
    t(:,2,:)=reshape(x((f_max*4+1):size(x,1),:),t_max,1,M);
    Pt=zeros(t_max,1,M);
    Pt(:,2,:)=reshape(P((f_max*4+1):size(P,1),:),t_max,1,M);
end;
