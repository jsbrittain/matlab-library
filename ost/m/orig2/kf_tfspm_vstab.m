function [sp,Psp]=kf_tfspm_vstab(sp,cl,Q);
%
% (REWRITE) Kalman-Fourier analysis for variance stabilised estimates
%
% Performs a log10 transform on sp matrix (variance stabilising).
% Returns the coherence calcualted from the sp matrix as a 4th column.
%
% Input parameters
%       sp      Spectral matrix of Fourier coefficients.
%       cl      Confidence limits structure for Fourier analysis.
%       Q       (opt) Process noise variance (default 1e-5).
%
% Output parameters
%       sp      Spectral matrix, cols; (Dim 3 represent group entities)
%                   1       Auto-spectra ch.1
%                   2       Auto-spectra ch.2
%                   3       Cross-spectra ch.1,2
%                   4       Coherence
%       Psp     Covariance matrix corresponding to sp matrix.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% Smoothing in log10 space produces biased estimates, even if taking %
% the sample mean!                                                   %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine data parameters
M=size(sp,3);       % Group count
f_max=size(sp,1);   % Number of frequency points
vlen=3*f_max;       % Length of state vector
sp=log10(sp);       % Variance stabilising transform

% Determine process noise covariance (Q)
if (exist('Q'))
    if (length(Q)==1), Q=Q.*ones(vlen,1); end;
else
    Q=1e-5*ones(vlen,1);        % Default value
end;

% Initialise variable space
x=zeros(vlen,M);
P=zeros(vlen,M);

% Initial conditions
xp=reshape(sp(:,:,1),vlen,1);      % Initial state vector
R_const=ones(vlen,1)*(log10(exp(1))^2);
R=R_const./cl(1).seg_tot;

Pp=R;                       % 'a priori' error-covar matrix
theta=ones(vlen,1);         % State-transition matrix
H=ones(vlen,1);             % Process-to-measurement matrix

x(:,1)=xp;
P(:,1)=Pp;
for ind=2:M
    % Construct state vector
    newz=reshape(sp(:,:,ind),vlen,1);
    % Construct measurement noise covariance
    R=R_const./cl(ind).seg_tot;
    % Kalman filter
    disp(['Running kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    [x(:,ind),P(:,ind),K,Pp,xp]=kf_filtr(newz,Q,R,theta,H,Pp,xp);
end;

	% Create a running mean to check function
	ss11=sp(:,1,1);
	ss22=sp(:,2,1);
	ss12=sp(:,3,1);
	for ind=2:M
        disp(['Running mean filter for group ' int2str(ind) ' of ' int2str(M)]);
        ss11=(ind/(ind+1))*(ss11+(1/ind)*sp(:,1,ind));
        ss22=(ind/(ind+1))*(ss22+(1/ind)*sp(:,2,ind));
        ss12=(ind/(ind+1))*(ss12+(1/ind)*sp(:,3,ind));
        x(:,ind)=[ss11; ss22; ss12];
	end;

% Formulate frequency domain output
sp=reshape(x,f_max,3,M);
Psp=reshape(P,f_max,3,M);

% Calculate coherence
for ind=1:M
    sp(:,4,ind)=(abs(10.^sp(:,3,ind)).^2)./((10.^sp(:,1,ind)).*(10.^sp(:,2,ind)));
end;
