function [sp,params]=kf_tfspm2w_2(sp,cl,Q,smooth,jackknife,jazalph)
%function [sp,params]=kf_tfspm2w_2(sp,cl,Q,smooth,jackknife,jazalph)
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
% vectors used where possible).
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
%                           2       Auto-spectra ch.2
%                           3       Cross-spectra ch.1,2
%                           4       Coherence
%       Psp             Estimated error covariance matrix corresponding to sp.
%
%function [sp,params]=kf_tfspm2w_2(sp,cl,Q,smooth,jackknife,jazalph)

% Determine data parameters
M=size(sp,3);           % Group count
fpts=size(sp,1);        % Number of frequency points
vlen=4*fpts;            % Length of state vector
N=cl(1).pad_len;        % Max segment length

% Determine process noise covariance (Q)
if (~exist('Q'))
    Q=0.1/100;          % Default value
end;
if (~exist('smooth'))
    smooth=false;
end;
if (~exist('jackknife'))
    jackknife=false;
end;
if (~exist('jazalph'))
    jazalph = 1;        % No process variance correction
end;

% Initialise variable space
x=zeros(vlen,M);
P=zeros(1,M);
xp=zeros(vlen,M+1);
Pp=zeros(1,M+1);
L=zeros(1,M);
w=zeros(1,M);
evidence=zeros(vlen,M);
if (jackknife)
    jk=zeros(fpts,5,M);
    jkmean=zeros(fpts,5,M);
end;

% Transform spectra to log domain
warning off     % Prevent warning on log(0) - dealt with below
sp(:,1,:)=single(log(double(sp(:,1,:))));
sp(:,2,:)=single(log(double(sp(:,2,:))));
sp(:,4,:)=single(log(imag(double(sp(:,3,:)))));                     % NB: real/imag cross-spectra can be -ve leading
sp(:,3,:)=single(log(real(double(sp(:,3,:)))));                     %  |  to complex logarithms (filter seperately).
warning on
varstab=1;%log(exp(1))^2;                                           % Variance stabilising constant (for various log bases)

% Set log(0) values to -100 (log(0)=-Inf therefore test real part only. Imag part only ever {0,pi})
sp(real(sp)<-100)=-100+i*imag(sp(real(sp)<-100));

% Transition matrices
theta=1;                                                            % State-transition matrix
H=1;                                                                % Process-to-measurement matrix

% Initial conditions (largely redundant due to first iteration)
xp(:,1)=reshape(double(sp(:,:,1)),vlen,1);                          % Initial 'a priori' state vector
R=varstab*psi(1,cl(1).seg_tot);%*(N/cl(1).seg_len);                   % Measurement error-covariance matrix
Pp(:,1)=R;                                                          % Initial 'a priori' error cov

% First iteration (ensures correct variance,L returned)
x(:,1)=xp(:,1); P(:,1)=R;
xp(:,2)=x(:,1); Pp(:,2)=P(:,1)+Q;
L(1)=1; w(1)=1;
if (jackknife)
    % jk(:,:,1)=zeros(fpts,5,1);      % Zero variance for single point
    %jkmean(:,:,1)=;
end;
params.w(1,:)=single(w);

% Kalman filter - forward sweep
progress=0;
for ind=2:M
    % Display progress
    if (floor(ind/M*100)>progress)
        progress=floor(ind/M*100);
        disp(['Running kalman filter: ' int2str(progress) '%']);
    end;
    % Construct state vector
    newz=reshape(double(sp(:,:,ind)),vlen,1);
    % Construct measurement noise covariance (from updated spectrum)
    R=varstab*psi(1,cl(ind).seg_tot);%*(N/cl(ind).seg_len);
    % Evidence
    evidence(:,ind) = exp(-((newz-xp(:,ind)).^2)./(2*R));
    % Jazwinski algorithm (update Q)
    enewz = (newz(1:end/2)-xp(1:end/2,ind));
    qj = ( mean(enewz).^2 - (R+P(1,ind-1)) ); qj(qj<0)=0;
    Q = jazalph*Q + (1-jazalph)*qj;
    Pp(:,ind) = P(:,ind-1)+Q;                                               % Update `a priori' state covar
    % Kalman filter auto-spectra
    [x(:,ind),kfP,K,kfPp,xp(:,ind+1)]=kalm_filtrv(newz,Q,R,theta,H,Pp(:,ind),xp(:,ind));
    K=K(1); P(:,ind)=kfP(1); Pp(:,ind+1)=kfPp(1);                           % Propagate scalars (constant K,P,Pp over freqs)
    % Construct weight matrix to determine effective no. segments
    Ka=K; Kb=(1-K);
    w(:,1:ind)=[Kb(1:end*(~isempty(w)),ones(1,ind-1)).*w(:,1:ind-1) Ka];    % Weights maintained to sum=1
    if ((~smooth) | (ind==M))
        L(ind)=1./sum(w.^2,2);
        % Construct jackknife error estimates
        if (jackknife)
            jkrange = (1:ind);
            [jk(:,:,ind),jkmean(:,:,ind)]=kflog_jk(double(sp(:,:,jkrange)),w(jkrange));
        end;
    end;
    % Keep record of weighting scheme
    params.w(ind,:)=single(w);
end;

% Kalman smoother - backward sweep
if (smooth)
    % Allocate memory
	xN=zeros(vlen,M);
	PN=zeros(vlen,M);
    % Initial conditions of smoother
    xN(:,M)=x(:,M);
    PN(:,M)=P(:,M);
    % Init backward weight matrix
    w2=w;
    % Recurse backward sweep
    progress=0;
	for ind=M-1:-1:1
        % Display progress
        if (floor(ind/M*100)<progress)
            progress=floor(ind/M*100);
            disp(['Running kalman filter: ' int2str(100-progress) '% (B)']);
        end;
        % Determine backward gain and combine forward/backward state estimates
        A=P(:,ind)./Pp(:,ind+1);
        xN(:,ind)=x(:,ind)+A.*(xN(:,ind+1)-xp(:,ind+1));
        PN(:,ind)=P(:,ind)+A.*(PN(:,ind+1)-Pp(:,ind+1)).*A;
        
        % Adjust effective no. segments for backward sweep
        Aa=A; Ab=(1-A);
        % Re-normalise forward weights
        wf=[w(:,1:ind) zeros(1,M-ind)];
        wsum=sum(wf,2);
        wf=wf./wsum(:,ones(1,M));
        % Normalise backward weights
        wb=w2; wsum=sum(wb,2);
        wb=wb./wsum(:,ones(1,M));
        % Determine effective no. segments (normalised bwd coeffs)
        nw=wb.*Aa(:,ones(1,size(wb,2)))+wf.*Ab(:,ones(1,size(wf,2)));
        L(:,ind)=1./sum(nw.^2,2);
        % Update backward sweep coefficients (un-normalised bwd coeffs)
        w2=w2.*Aa(:,ones(1,size(wb,2)))+wf.*Ab(:,ones(1,size(wf,2)));
        % Construct jackknife error estimates
        if (jackknife)
            jk(:,:,ind)=kflog_jk(double(sp),nw);
        end;
        % Keep record of weighting scheme
        params.w(ind,:)=single(nw);
    end;
	x=xN; P=PN; clear('xN','PN');
end;
clear('w','w2');
L=single(L(ones(3*fpts,1),:));

% Formulate frequency domain output
x12=real(exp(x(2*fpts+1:3*fpts,:)))+i*real(exp(x(3*fpts+1:end,:)));    % (complex) cross-spectrum (time-domain)
sp=reshape([exp(x(1:2*fpts,:)); x12],fpts,3,M);
P=P(ones(3*fpts,1),:);
params.Psp=single(reshape(P,fpts,3,M));
params.evidence=evidence;
clear('x','P','xp','Pp','x12');

% Generate parameters structure
params.L=single(L);
if (jackknife)
    params.jk=jk; clear('jk');
    params.jkmean=jkmean; clear('jkmean');
end;

% Calculate coherence & phase
warning off MATLAB:divideByZero
sp(:,4,:)=abs(sp(:,3,:)).^2./((sp(:,1,:)).*(sp(:,2,:)));
warning on MATLAB:divideByZero
sp(:,5,:)=angle(sp(:,3,:));
