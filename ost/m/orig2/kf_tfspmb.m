function [f,Pf,t,Pt,sp,Psp]=kf_tfspm(f,arg1,arg2,arg3,arg4,arg5);
%function [f,Pf,t,Pt,sp,Psp]=kf_tfspm(f,[t],cl,[sp],[f_max],[Q]);
%
% Hybrid Kalman filter and Fourier based analysis
% Kalman filtering of Fourier-based coherence estimates over trials.
% Takes f and t matrices as input.
%
% Reduced vector implementation returning frequency-domain parameters,
% associated error covariance parameters, the set of time-domain parameters
% and associated error-covariance matrices along with all spectral
% coefficients from the sp matrix 
%
% NB: No optional parameters may be skipped if subsequent parameters are
%     provided.
%
%function [f,Pf,t,Pt,sp,Psp]=kf_tfspm(f,[t],cl,[sp],[f_max],[Q]);

% Design time options
stab_coh=1;                 % Variance stabilizing transform for coherence

% Check parameter validity
tout=0;
spout=0;
if (nargin<4)
    error('You must provide at least the first three input arguments.');
end;
if ((mod(nargout,2)~=0))
    error('Output arguments must be returned in pairs.');
end;
if (nargout>2)
    tout=1;     % Assess time-domain parameters
end;
if (nargout>4)
    spout=1;    % Assess spectral coefficients
end;

% Determine input parameters
f_max=-1;
Q=-1;
switch (nargin)
    case 2
        cl=arg1;
    case 3
        if (size(arg1,3)>1)     % if a t matrix
            t=arg1;
            cl=arg2;
        else
            cl=arg1;
            f_max=arg2;
        end;
    case 4
        if (size(arg1,3)>1)     % if a t matrix
            t=arg1;
            cl=arg2;
            if (size(arg3,2)>1) % if a sp matrix
                sp=arg3;
            else
                f_max=arg3;
            end;
        else
            cl=arg1;
            f_max=arg2;
            Q=arg3;
        end;
    case 5
        if (size(arg3,1)>1)     % if arg3 is an sp matrix
            t=arg1;
            cl=arg2;
            sp=arg3;
            f_max=arg4;
        else
            t=arg1;
            cl=arg2;
            f_max=arg3;
            Q=arg4;
        end;
    case 6
        t=arg1;
        cl=arg2;
        sp=arg3;
        f_max=arg4;
        Q=arg5;
    otherwise
        error('Incorrect number of input arguments');
end;

% Determine data parameters
M=size(f,3);
t_max=size(t,1);

% Determine f_max in terms of samples
if (isempty(f_max))
    f_max=-1;
end;
if (f_max~=-1)
    f_first=find(f(:,1,1)>f_max);
    if (isempty(f_first))
        f_max=size(f,1);
    else
        f_max=f_first(1)-1;
    end;
    if (f_max==0)
        error('Specified f_max too small.');
    end;
else
    f_max=size(f,1);
end;

% Determine state vector length
vlen=4*f_max+tout*t_max+spout*3*f_max;

% Determine Q
if (Q==-1)
    Q=1e-5*ones(vlen,1);    % Default noise covariance
else
    if (size(Q,1)==1)
        Q=Q*ones(vlen,1);
    end;
end;

% Initialise variable space
x=zeros(vlen,M);
P=zeros(vlen,M);

% Kalman fitering parameters and initial conditions
coh=reshape(f(1:f_max,4,1),f_max,1);
Rspect=ones(f_max,1)*(0.4343^2)/cl(1).seg_tot;
if (stab_coh)   % Variance stabilising transform, propagate atanh mag. coherency
    xp=[reshape(f(1:f_max,2:3,1),2*f_max,1); atanh(sqrt(coh)); reshape(f(1:f_max,5,1),f_max,1)];
    Rcoh=ones(f_max,1)/(2*cl(1).seg_tot);       % PBMB (6.5)
else
    xp=reshape(f(1:f_max,2:5,1),4*f_max,1);          % Initial estimates
    Rcoh=(2/cl(1).seg_tot)*(coh.*((1-coh).^2));   % PBMB (6.4)
end;
Rphase=(1/(2*cl(1).seg_tot))*((1./coh)-1);
if (tout)
    xp=[xp; t(:,2,1)];
    Rcum=ones(t_max,1)*((cl(1).q_c95/1.96)^2);
else
    Rcum=[];
end;
if (spout)
    xp=[xp; reshape(sp(1:f_max,:,1),3*f_max,1)];
    Rauto1=(sp(1:f_max,1,1).^2)/cl(1).seg_tot;
    Rauto2=(sp(1:f_max,2,1).^2)/cl(1).seg_tot;
    Rcross=(sp(1:f_max,3,1).^2)/cl(1).seg_tot;
else
    Rauto1=[];
    Rauto2=[];
    Rcross=[];
end;
R=[Rspect; Rspect; Rcoh; Rphase; Rcum; Rauto1; Rauto2; Rcross];

Pp=R;                                              % State error-covariance matrix
theta=ones(vlen,1);                                % State-transition matrix
H=ones(vlen,1);                                    % Process-to-measurement matrix

x(:,1)=xp;
P(:,1)=Pp;
for ind=2:M
    % Construct state vector
    coh=reshape(f(1:f_max,4,ind),f_max,1);
    if (stab_coh)
        newz=[reshape(f(1:f_max,2:3,ind),2*f_max,1); atanh(sqrt(coh)); reshape(f(1:f_max,5,ind),f_max,1)];
    else
        newz=reshape(f(1:f_max,2:5,ind),4*f_max,1);
        Rcoh=(2/cl(ind).seg_tot)*(coh.*((1-coh).^2));
    end;
    
    % Construct measurement noise covariance matrix
    Rphase=(1/(2*cl(ind).seg_tot))*((1./coh)-1);
    if (tout)
        newz=[newz; t(:,2,ind)];
        Rcum=ones(t_max,1)*((cl(ind).q_c95/1.96)^2);
    end;
    if (spout)
        newz=[newz; reshape(sp(1:f_max,:,ind),3*f_max,1)];
        Rauto1=(sp(1:f_max,1,ind).^2)/cl(1).seg_tot;
        Rauto2=(sp(1:f_max,2,ind).^2)/cl(1).seg_tot;
        Rcross=(sp(1:f_max,3,ind).^2)/cl(1).seg_tot;
    end;
    R=[Rspect; Rspect; Rcoh; Rphase; Rcum; Rauto1; Rauto2; Rcross];
    
    % Perform Kalman filter
    disp(['Running Kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    [x(:,ind), P(:,ind), K, Pp, xp] = kf_filtr(newz, Q, R, theta, H, Pp, xp);
end;

% Inverse variance stabalizing transform
if (stab_coh)
    start=4*f_max+1;
    stop=5*f_max;
    x(start:stop,:)=tanh(x(start:stop,:)).^2;
%    P(start:stop,:)=P(start:stop,:);           % Inverse transform variance
end;

% Formulate frequency-domain matrix output
f=f(1:f_max,1,:);
f(:,2:5,:)=reshape(x(1:(f_max*4),:),f_max,4,M);
Pf=zeros(f_max,1,M);
Pf(:,2:5,:)=reshape(P(1:(f_max*4),:),f_max,4,M);

% Formulate time-domain matrix output
if (tout)
    t=t(1:t_max,1,:);
    start=f_max*4+1;
    stop=start+t_max-1;
    t(:,2,:)=reshape(x(start:stop,:),t_max,1,M);
    Pt=zeros(t_max,1,M);
    Pt(:,2,:)=reshape(P(start:stop,:),t_max,1,M);
end;

% Formulate spectral coefficients
if (spout)
    start=stop+1;
    stop=start+3*f_max-1;
    sp=reshape(x(start:stop,:),f_max,3,M);
    Psp=reshape(P(start:stop,:),f_max,3,M);
end;
