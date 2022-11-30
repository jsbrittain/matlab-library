function [f,Pf,t,Pt]=kf_spm(f,t,cl,f_max,Q);
%function [f,Pf,t,Pt]=kf_spm(f,t,cl,f_max,Q);
%
% Hybrid Kalman filter and Fourier based analysis
% Kalman filtering of Fourier-based coherence estimates over trials.
%
% This version permits the evaluation of time-domain parameters
% via a second Kalman filter (done this way for efficiency)
% NB: Current implementation still exceptionally slow for time-domain.
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
    tout=1;     % Calculate time-domain parameters
else
    tout=0;
end;

% Determine input parameters
M=size(f,3);
t_max=size(t,1);
if (nargin<4)
    f_max=size(f,1);
else
    f_max=find(f(:,1,1)>f_max);
    f_max=f_max(1)-1;
end;
if (nargin<5)
    Qx=1e-5*eye(4*f_max);  % Default noise covariance
    if (tout), Qy=1e-5*eye(t_max); end;
else
    if (size(Q,1)==1)
        Qx=Q*eye(f_max);
        Qy=Q*eye(f_max);
    else
        Qx=Q;
        Qy=Q;
    end;
end;

% Initialise variable space
x=zeros(4*f_max,M);
Px=zeros(4*f_max,M);
if (tout)
    y=zeros(t_max,M);
    Py=zeros(t_max,M);
end;

% Kalman fitering parameters and initial conditions
coh=reshape(f(1:f_max,4,1),f_max,1);
xp=[reshape(f(1:f_max,2:5,1),4*f_max,1)];   % Initial estimates

Rspect=ones(f_max,1)*(0.4343^2)/cl(1).seg_tot;                  % PBMB (6.1)
Rcoh=(2/cl(1).seg_tot)*((abs(coh).^2).*((1-(abs(coh).^2)).^2)); % PBMB (6.4)
Rphase=(1/(2*cl(1).seg_tot))*((1./(abs(coh).^2))-1);            % PBMB (6.7)
Rx=diag([Rspect; Rspect; Rcoh; Rphase],0);

Ppx=Rx;                                      % State error-covariance matrix
thetax=eye(4*f_max);                         % State-transition matrix
Hx=eye(4*f_max);                             % Process-to-measurement matrix

x(:,1)=xp;                                   % Assign a priori estimate to
Px(:,1)=diag(Ppx);                           % first output.

% Time-domain filtering
if (tout)
	yp=t(:,2,1);                                 % Initial state estimate
	Ry=diag(ones(t_max,1)*((cl(1).q_c95/1.96)^2),0); % PBMB (6.10)
	Ppy=Ry;                                      % State error-covariance matrix
	thetay=eye(t_max);                           % State-transition matrix
	Hy=eye(t_max);                               % Process-to-measurement matrix
    y(:,1)=yp;
    Py(:,1)=diag(Ppy);
end;

for ind=2:M
    disp(['Running Kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    
    % Frequency-domain analysis
    %===========================
    
    % Construct state vector
    coh=reshape(f(1:f_max,4,ind),f_max,1);
    newz=reshape(f(1:f_max,2:5,ind),4*f_max,1);
    
    % Construct measurement noise covariance matrix
    Rcoh=(2/cl(ind).seg_tot)*((abs(coh).^2).*((1-(abs(coh).^2)).^2));
    Rphase=(1/(2*cl(ind).seg_tot))*((1./(abs(coh).^2))-1);
    Rx=diag([Rspect; Rspect; Rcoh; Rphase],0);
    
    % Perform Kalman filter
    [x(:,ind), Pout, K, Ppx, xp] = kalm_filtr(newz, Qx, Rx, thetax, Hx, Ppx, xp);
    Px(:,ind)=diag(Pout);
    
    % Time-domain analysis
    %======================
    
    if (tout)
        % Construct state vector
        newz=reshape(t(:,2,ind),t_max,1);
        
        % Construct measurement noise covariance matrix
        Ry=diag(ones(t_max,1)*((cl(ind).q_c95/1.96)^2),0);
        
        % Perform Kalman filter
        [y(:,ind), Pout, K, Ppy, yp] = kalm_filtr(newz, Qy, Ry, thetay, Hy, Ppy, yp);
        Py(:,ind)=diag(Pout);
    end;
end;

% Formulate output (spectral coefficients)
f=f(1:f_max,1,:);
f(:,2:5,:)=reshape(x,f_max,4,M);
Pf=zeros(f_max,1,M);
Pf(:,2:5,:)=reshape(Px,f_max,4,M);

% Formulate output (time-domain coefficients)
if (tout)
	t=t(:,1,:);
	t(:,2,:)=reshape(y,t_max,1);
	Pt=zeros(t_max,1,M);
	Pt(:,2,:)=reshape(Py,t_max,1);
end;
