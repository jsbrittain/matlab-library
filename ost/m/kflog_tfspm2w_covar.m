function [sp,params]=kflog_tfspm2w_covar(sp,cl,Q,smooth,jackknife,jazalph)
%function [sp,params]=kflog_tfspm2w_covar(sp,cl,Q,smooth,jackknife,jazalph)
%
% Hybrid Kalman-Fourier analysis
%
% Performs kalman-filtering in the log-domain (For non-adaptive case this
% gives a running geometric mean).  Reduced vector implementation
% (without loss of generalisation) for speed.  Includes weighting on
% measurement variance to de-emphasise short segments.
%
% Kalman filtering of auto- and cross-spectra.
%
% This version optimised to minimise memory constraints (single precision
% vectors used where possible). Also, this version implements a
% frequency-selective Jazwinski algorithm (kf_tfspm2w_2 provides average
% error correction). Cross-spectra (necessarily) require reconstruction
% outside of the usual Kalman implementation since process noise variance
% becomes data-dependent.
%
% THIS VERSION WILL IMPLEMENT FULL COVARIANCE STRUCTURE IN PROCESS NOISE,
% WITH THE OPTION OF EXPANDING SINGLE-SEGMENT ESTIMATORS TO SIMILAR
% STRUCTURE (USEFUL FOR MULTITAPER ESTIMATES) - Q: JAZWINSKI
% IMPLEMENTATION?
%
% !! THIS IS A SINGLE CHANNEL SPECTRAL SMOOTHER APPROACH ONLY - SECOND CHANNEL IGNORED !!
%
% PS - I HAVEN'T DONE ANY OF THIS YET!
%
% Input parameters
%       sp              Spectral matrix of Fourier coefficients
%       cl              Confidence limits structure for Fourier analysis
%       Q               (opt) Process noise variance (default 1e-3)
%       smooth          (opt) Smooth spectra (1=Yes 0=No(default))
%       jackknife       (opt) Jackknife variance (1=Yes 0=No(default))
%       jazalph         (opt) Alpha-rate for Jazwinski adaptive noise algorithm
%
% Output parameters
%       sp              Spectral matrix, cols: (Dim 3 represent group entities)
%                           1       Auto-spectra ch.1
%       Psp             Estimated error covariance matrix corresponding to sp.
%
%function [sp,params]=kflog_tfspm2w_covar(sp,cl,Q,smooth,jackknife,jazalph)

% Determine data parameters
M=size(sp,3);           % Group count
fpts=size(sp,1);        % Number of frequency points
vlen=fpts;              % Length of state vector

% Determine process noise covariance (Q)
if (~exist('Q','var'))
    Q=0.1/100;          % Default value
end;
if (~exist('smooth','var'))
    smooth=false;
end;
if (~exist('jackknife','var'))
    jackknife=false;
end;
if (~exist('jazalph','var'))
    jazalph = [];
end;
if (isempty(jazalph))
    jazalph = 1;                % No process variance correction
end;

% Working class precision
classstr='single';

% Initialise variable space
x=zeros(vlen,M,classstr);
P=zeros(vlen,vlen,M,classstr);
xp=zeros(vlen,M+1,classstr);
Pp=zeros(vlen,vlen,M+1,classstr);
L=zeros(vlen,M,classstr);
w=zeros(vlen,M,classstr);
params.w=zeros(M,vlen,M,classstr);
evidence=zeros(vlen,M,classstr);
if (jackknife)
    jk=zeros(fpts,5,M,classstr);
    jkmean=zeros(fpts,5,M,classstr);
end;

% Transform spectra to log domain
warning('off','MATLAB:log:logOfZero');                              % Prevent warning on log(0) - dealt with below
sp(:,1,:)=log(sp(:,1,:));
warning('on','MATLAB:log:logOfZero');
varstab=1;%log(exp(1))^2;                                           % Variance stabilising constant (for various log bases)

% Transition matrices
theta=1;                                                            % State-transition matrix
H=1;                                                                % Process-to-measurement matrix

% Initial conditions (largely redundant due to first iteration)
Q=Q*eye(vlen,vlen,classstr);                                        % State error-covariance matrix (full)
R=eye(vlen,classstr)*varstab*psi(1,cl(1).seg_tot);                  % Measurement error-covariance matrix (diagonal)
xp(:,1)=sp(:,1,1);                                                  % Initial 'a priori' state vector
Pp(:,:,1)=R;                                                        % Initial 'a priori' error cov

% First iteration (ensures correct variance/L returned)
x(:,1)=xp(:,1);
xp(:,2)=x(:,1);
P(:,:,1)=R;
Pp(:,:,2)=P(:,:,1) + Q;
L(:,1)=ones(vlen,1); w(:,1)=ones(vlen,1);
if (jackknife)
    jk(:,:,1)=zeros(fpts,5,1);      % Zero variance for single point
end;
params.w(1,:,:)=w;

% Kalman filter - forward sweep
progress=0;
for ind=(2:M)
    % Display progress
    if (floor(ind/M*100)>progress)
        progress=floor(ind/M*100);
        disp(['Running kalman filter: ' int2str(progress) '%']);
    end;
    % Construct state vector
    newz=sp(:,1,ind);
    % Construct measurement noise covariance (from updated spectrum)
    R=eye(vlen,classstr)*varstab*psi(1,cl(ind).seg_tot);
    % Evidence
    evidence(:,ind) = exp(-((newz-xp(:,ind)).^2)./(2*diag(R)));
    % Jazwinski algorithm (update Q)
    enewz = (newz-xp(:,ind));
    qj = ( (enewz*(enewz.')) - (R+P(:,:,ind-1)) ); qj(qj<0)=0;
    
    for n=(1:size(qj,1))
        if (qj(n,n)<0)
            qj(:,n)=0; qj(n,:)=0;
        end;
    end;
    
    Q = jazalph*Q + (1-jazalph)*qj;
    Pp(:,:,ind) = P(:,:,ind-1) + Q;             % Update `a priori' state covar
    
    % Kalman filter auto-spectra
    [x(:,ind),P(:,:,ind),K,Pp(:,:,ind+1),xp(:,ind+1)] = kalm_filtr(newz,Q,R,theta,H,Pp(:,:,ind),xp(:,ind));
    
    % Construct weight matrix to determine effective no. segments
    Ka=diag(K); Kb=diag(1-K);                   % Appended is cross-spec weighting
    
    %%% REALLY UGLY SOLUTION %%%
    
    w(:,1:ind)=[Kb(:,ones(1,ind-1)).*w(:,(1:(ind-1))) Ka];
    if ((~smooth) || (ind==M))
        % Force cross-spectral weights to sum=1 (as not kalman filtered)
        w2=w(:,1:ind);
        % Determine effective no. segments (as gains sums to 1 at each step)
        L(:,ind)=1./sum(w2.^2,2);
        % Construct jackknife error estimates
        if (jackknife)
            jkrange = (1:ind);
            [jk(:,:,ind),jkmean(:,:,ind)]=kflog_jk(sp(:,:,jkrange),w(jkrange));
        end;
    end;
    % Keep record of weighting scheme
    params.w(ind,:,:)=w;
end;

% Kalman smoother - backward sweep
if (smooth)
    % Allocate memory
	xN=zeros(vlen,M,classstr);
	PN=zeros(vlen,vlen,M,classstr);
    % Initial conditions of smoother
    xN(:,M)=x(:,M);
    PN(:,:,M)=P(:,:,M);
    % Init backward weight matrix
    w2=w;
    % Recurse backward sweep
	for ind=(M-1:-1:1)
        % Display progress
        if (floor(ind/M*100)<progress)
            progress=floor(ind/M*100);
            disp(['Running kalman filter: ' int2str(100-progress) '% (B)']);
        end;
        % Determine backward gain and combine forward/backward state estimates
        A=P(:,:,ind)*inv(Pp(:,:,ind+1));
        xN(:,ind)=x(:,ind)+A*(xN(:,ind+1)-xp(:,ind+1));
        PN(:,:,ind)=P(:,:,ind)+A*(PN(:,:,ind+1)-Pp(:,:,ind+1))*(A.');
        
        % Adjust effective no. segments for backward sweep
        Aa=A; Ab=(1-A);
        % Re-normalise forward weights
        wf=[w(:,1:ind) zeros(vlen,M-ind)];
        wsum=sum(wf,2);
        wf=wf./wsum(:,ones(1,M));
        % Normalise conditioned cross-spectral weights
        wb=w2;
        % Determine effective no. segments (normalised bwd coeffs)
        nw=wb.*Aa(:,ones(1,size(wb,2)))+wf.*Ab(:,ones(1,size(wf,2)));
        L(:,ind)=1./sum(nw.^2,2);
        % Update backward sweep coefficients (un-normalised bwd coeffs)
        w2=w2.*Aa(:,ones(1,size(wb,2)))+wf.*Ab(:,ones(1,size(wf,2)));
        % Construct jackknife error estimates
        if (jackknife)
            jk(:,:,ind)=kflog_jk(sp,nw);
        end;
        % Keep record of weighting scheme
        params.w(ind,:,:)=nw;
	end;
	x=xN; P=PN; clear('xN','PN','x12N');
end;
clear('w','w2');

% Formulate frequency domain output
sp=reshape(exp(x),fpts,1,M);
params.Psp=P;
params.evidence=evidence;
clear('x','P','xp','Pp','x12');

% Generate parameters structure
params.L=[L; L; L];
if (jackknife)
    params.jk=jk; clear('jk');
    params.jkmean=jkmean; clear('jkmean');
end;

% Calculate coherence & phase
warning off MATLAB:divideByZero
%sp(:,4,:)=abs(sp(:,3,:)).^2./((sp(:,1,:)).*(sp(:,2,:)));
warning on MATLAB:divideByZero
%sp(:,5,:)=angle(sp(:,3,:));
