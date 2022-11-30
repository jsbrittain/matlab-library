function [f,Pf,t,Pt]=kf_tfmfb(f,t,cl,f_max,Q);
%function [f,Pf,t,Pt]=kf_tfmfb(f,t,cl,f_max,Q);
%
% Hybrid Kalman filter and Fourier based analysis
% Kalman filtering of Fourier-based coherence estimates over trials
% with forward-backward smoothing.
% Takes f and t matrices as input.
%
% Reduced vector implementation returning frequency-domain parameters,
% associated error covariance parameters, and optionally the set of
% time-domain parameters and associated error-covariance parameters.
%
%function [f,Pf,t,Pt]=kf_tfmfb(f,t,cl,f_max,Q);

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

% Perform forward filtering
[f1,Pf1,t1,Pt1]=kf_tfm(f,t,cl,f_max,Q);

% Formulate last state estimate for backward recursion
f_max=size(f1,1);
t_max=size(t1,1);
M=size(f,3);
x0=zeros(4*f_max+t_max,1);
P0=zeros(4*f_max+t_max,1);
x0(1:(4*f_max))=reshape(f1(:,2:5,M),4*f_max,1);
x0((4*f_max+1):(4*f_max+t_max))=reshape(t1(:,2,M),t_max,1);
P0(1:(4*f_max))=reshape(Pf1(:,2:5,M),4*f_max,1);
P0((4*f_max+1):(4*f_max+t_max))=reshape(Pt1(:,2,M),t_max,1);

% Perform backward recursion
[f2,Pf2,t2,Pt2]=kf_tfm(flipdim(f,3),flipdim(t,3),flipdim(cl,1),f_max,Q,x0,P0);

% Flip backward estimates to correspond to forward estimates
f2=flipdim(f2,3);
Pf2=flipdim(Pf2,3);
t2=flipdim(t2,3);
Pt2=flipdim(Pt2,3);

% Allocate variable space (and overwrite input parameters)
f=zeros(size(f1));
Pf=zeros(size(Pf1));
t=zeros(size(t1));
Pt=zeros(size(Pt1));

% Fill first column of f and t with freq/time data
f(:,1,:)=f1(:,1,:);
t(:,1,:)=t1(:,1,:);

% Ensure no zero valued data
vsmall=eps;
f1(f1==0)=vsmall; f2(f2==0)=vsmall;
Pf1(Pf1==0)=vsmall; Pf2(Pf2==0)=vsmall;
t1(t1==0)=vsmall; t2(t2==0)=vsmall;
Pt1(Pt1==0)=vsmall; Pt2(Pt2==0)=vsmall;

% Combine into a single smoothed estimate
for ind=1:size(f1,3)
    % Frequency domain parameters
    Pf(:,2:5,ind)=1./((1./Pf1(:,2:5,ind)) + 1./conj(Pf2(:,2:5,ind)));
    f(:,2:5,ind)=Pf(:,2:5,ind).*((1./Pf1(:,2:5,ind)).*f1(:,2:5,ind) + (1./conj(Pf2(:,2:5,ind))).*conj(f2(:,2:5,ind)));
    
    % Time domain parameters
    Pt(:,2,ind)=1./((1./Pt1(:,2,ind)) + 1./conj(Pt2(:,2,ind)));
    t(:,2,ind)=Pt(:,2,ind).*((1./Pt1(:,2,ind)).*t1(:,2,ind) + (1./conj(Pt2(:,2,ind))).*conj(t2(:,2,ind)));
end;
