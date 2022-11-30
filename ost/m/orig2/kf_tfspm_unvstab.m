function [sp,Psp]=kf_tfspm_unvstab(sp,cl,power,Q,smooth);
%
% (REWRITE) Kalman-Fourier analysis, un-variance stabilised (named to avoid
% clashing with kf_tfspm), this version does not variance stabilise as this
% would lead to logarithmic averaging, instead it relies on a disjoint
% sections analysis to provide an estimate of the true power spectrum which
% is subsequently used to provide variance estimates for the individual
% periodograms.
%
% NB: This methodology throws up some contradictions, for
% example, assuming a non-stationary spectrum then using the stationary
% spectrum as a parameter in the non-stationary analysis!
%
% Returns the coherence calcualted from the sp matrix as a 4th column.
%
% Input parameters
%       sp              Spectral matrix of Fourier coefficients.
%       cl              Confidence limits structure for Fourier analysis.
%       power           Spectral power estimated over all segments (disjoint analysis).
%       Qpercentage     (opt) Process noise variance (as a percentage of R) (default 0.1).
%       smooth          (opt) Smooth spectra (Backward filter spectra after forward filtering, 1=Yes 0=No).
%
% Output parameters
%       sp              Spectral matrix, cols; (Dim 3 represent group entities)
%                           1       Auto-spectra ch.1
%                           2       Auto-spectra ch.2
%                           3       Cross-spectra ch.1,2
%                           4       Coherence
%       Psp             Covariance matrix corresponding to sp matrix.
%

% Determine data parameters
M=size(sp,3);       % Group count
f_max=size(sp,1);   % Number of frequency points
vlen=3*f_max;       % Length of state vector

% Determine process noise covariance (Q)
if (exist('Q'))
    if (length(Q)==1), Q=Q.*ones(vlen,1); end;
else
    Q=0.1*ones(vlen,1);        % Default value
end;
if (~exist('smooth'))
    smooth=0;
end;

% Initialise variable space
x=zeros(vlen,M);
P=zeros(vlen,M);
xp=zeros(vlen,M+1);
Pp=zeros(vlen,M+1);

% Check if power a single estimate
if (size(power,3)==1)
    power(:,:,1:M)=power;
end;

% Initial conditions
xp(:,2)=reshape(sp(:,:,1),vlen,1);              % Initial state vector
R=reshape(power(:,:,1).^2,vlen,1)./cl(1).seg_tot;      % Measurement error-covariance matrix

Pp(:,2)=R;                  % 'a priori' error-covar matrix
theta=ones(vlen,1);         % State-transition matrix
H=ones(vlen,1);             % Process-to-measurement matrix

x(:,1)=xp(:,2);
P(:,1)=R;
for ind=2:M
    % Construct state vector
    newz=reshape(sp(:,:,ind),vlen,1);
    % Construct measurement noise covariance
    R=reshape(power(:,:,ind).^2,vlen,1)./cl(ind).seg_tot;
    % Kalman filter
    disp(['Running kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    [x(:,ind),P(:,ind),K,Pp(:,ind+1),xp(:,ind+1)]=kf_filtrv(newz,Q.*R,R,theta,H,Pp(:,ind),xp(:,ind));
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
sp=reshape(x,f_max,3,M);
Psp=reshape(P,f_max,3,M);

% Calculate coherence
for ind=1:M
    sp(:,4,ind)=(abs(sp(:,3,ind)).^2)./((sp(:,1,ind)).*(sp(:,2,ind)));
end;
