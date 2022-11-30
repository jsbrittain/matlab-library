function [sp,Psp]=kf_tfspm2w_mmd(sp,cl,Q,smooth);
%
% Hybrid Kalman-Fourier analysis
% Utilises a Dynamic Multi-Modal Model
%
% Performs periodogram filtering to produce spectral statistics across
% non-stationary (or assumed stationary) trials.  Includes self-startup by
% estimating the 'true'-local spectrum as the previous state estimate,
% beginning with the first periogram (with coherence=1), then continually
% smoothing and updating the spectra producing a running coherence
% estimate (assumes a slowly varying inter-trial spectrum).
%
% This version addresses the following issues:
%       startup         Periodogram estimate used for first variance
%       variance        Running spectral estimate used for successive
%                       variance estimates (removes disjoint analysis)
% Issues
%       Qpercentage     Need to justify this measure and not a scalar Q
%       distribution    Periodogram distribution not gaussian
%
% Input parameters
%       sp              Spectral matrix of Fourier coefficients.
%       cl              Confidence limits structure for Fourier analysis.
%       Qpercentage     (opt) Process noise variance (as a percentage of R) (default 0.1).
%       smooth          (opt) Smooth spectra (Backward filter spectra after forward filtering, 1=Yes 0=No).
%
% Output parameters
%       sp              Spectral matrix, cols: (Dim 3 represent group entities)
%                           1       Auto-spectra ch.1
%                           2       Auto-spectra ch.2
%                           3       Cross-spectra ch.1,2
%                           4       Coherence
%       Psp             Estimated error covariance matrix corresponding to sp.
%

% Determine data parameters
M=size(sp,3);       % Group count
fpts=size(sp,1);   % Number of frequency points
vlen=3*fpts;       % Length of state vector

% Determine process noise covariance (Q)
if (~exist('Q'))
    Q=0.1;        % Default value
end;
if (~exist('smooth'))
    smooth=0;
end;
filters=5;

% Initialise variable space
x=zeros(vlen,M);
P=zeros(vlen,M);
xp=zeros(vlen,filters);
Pp=zeros(vlen,filters);

% Initial conditions
xp(:,1)=reshape(sp(:,:,1),vlen,1);      % Initial state vector
R=(xp(:,1).^2)/cl(1).seg_tot;           % Measurement error-covariance matrix
R(2*fpts+1:vlen)=xp(1:fpts,1).*xp(fpts+1:2*fpts,1)/cl(1).seg_tot;
Pp(:,1)=R;

% Setup differences in filter bank
for ind=1:filters
    xp(:,ind)=xp(:,1);
    Pp(:,ind)=Pp(:,1);
end;

theta=ones(vlen,1);                     % State-transition matrix
H=ones(vlen,1);                         % Process-to-measurement matrix

for ind=1:M
    % Construct state vector
    newz=reshape(sp(:,:,ind),vlen,1);
    % Construct measurement noise covariance (from updated spectrum)
    if (ind==1)
        spec=newz;
    else
        spec=x(:,ind-1);
    end;
    % Cross-spectral distribution - Product of power spectra --- check ---
    R=(spec.^2)/cl(ind).seg_tot;
    R(2*fpts+1:vlen)=spec(1:fpts).*spec(fpts+1:2*fpts)/cl(ind).seg_tot;
    for j=1:filters
        R2(:,j)=R.*(10^(j-(filters+1)/2));
    end;
    % Kalman filter
    disp(['Running kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    [x(:,ind),P(:,ind),Pp,xp,p]=kf_mmdrv(newz,Q*R2,R2,theta,H,Pp,xp,vlen,filters);
end;

% Backward sweep
if (smooth==1)
	xN=zeros(vlen,M);
	PN=zeros(vlen,M);
    xN(:,M)=x(:,M);
    PN(:,M)=P(:,M);
	for ind=M-1:-1:1
        disp(['Running kalman filter for group ' int2str(ind) ' of ' int2str(M) ' (B)']);
        A=P(:,ind)./Pp(:,ind+1);
        xN(:,ind)=x(:,ind)+A.*(xN(:,ind+1)-xp(:,ind+1));
        PN(:,ind)=P(:,ind)+A.*(PN(:,ind+1)-Pp(:,ind+1)).*A;
	end;
	x=xN;
	P=PN;
end;

% Formulate frequency domain output
sp=reshape(x,fpts,3,M);
Psp=reshape(P,fpts,3,M);

% Calculate coherence
for ind=1:M
    sp(:,4,ind)=(abs(sp(:,3,ind)).^2)./((sp(:,1,ind)).*(sp(:,2,ind)));
end;
