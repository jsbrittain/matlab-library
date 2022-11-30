function [sp,params]=kflog_tfspm2w_3(sp,cl,Q,smooth,jackknife,jazalph)
%function [sp,params]=kflog_tfspm2w_3(sp,cl,Q,smooth,jackknife,jazalph)
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
%function [sp,params]=kflog_tfspm2w_3(sp,cl,Q,smooth,jackknife,jazalph)

% Determine data parameters
M=size(sp,3);           % Group count
fpts=size(sp,1);        % Number of frequency points
vlen=2*fpts;            % Length of state vector
%N=cl(1).pad_len;        % Max segment length

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
P=zeros(vlen,M,classstr);
xp=zeros(vlen,M+1,classstr);
x12=zeros(2*fpts,M,classstr);
xp12=zeros(2*fpts,M+1,classstr);
Pp=zeros(vlen,M+1,classstr);
L=zeros(3*fpts,M,classstr);
w=zeros(3*fpts,M,classstr);
params.w=zeros(M,3*fpts,M,classstr);
evidence=zeros(vlen,M,classstr);
if (jackknife)
    jk=zeros(fpts,5,M,classstr);
    jkmean=zeros(fpts,5,M,classstr);
end;

% Transform spectra to log domain
warning('off','MATLAB:log:logOfZero');                              % Prevent warning on log(0) - dealt with below
sp(:,1,:)=log(sp(:,1,:));
sp(:,2,:)=log(sp(:,2,:));
sp(:,4,:)=log(imag(sp(:,3,:)));                                     % NB: real/imag cross-spectra can be -ve leading
sp(:,3,:)=log(real(sp(:,3,:)));                                     %  |  to complex logarithms (filter seperately).
warning('on','MATLAB:log:logOfZero');
varstab=1;%log(exp(1))^2;                                           % Variance stabilising constant (for various log bases)

% Set log(0) values to -100 (log(0)=-Inf therefore test real part only. Imag part only ever {0,pi})
sp(real(sp)<-100)=-100+i*imag(sp(real(sp)<-100));

% Transition matrices
theta=1;                                                            % State-transition matrix
H=1;                                                                % Process-to-measurement matrix

% Initial conditions (largely redundant due to first iteration)
xp(:,1)=[sp(:,1,1); sp(:,2,1)];                                     % Initial 'a priori' state vector
xp12(:,1)=[sp(:,3,1); sp(:,4,1)];
R=varstab*psi(1,cl(1).seg_tot);                                     % Measurement error-covariance matrix
Pp(:,1)=R;                                                          % Initial 'a priori' error cov

% First iteration (ensures correct variance,L returned)
x(:,1)=xp(:,1);        xp(:,2)=x(:,1);        P(:,1)=R;
x12(:,1)=xp12(:,1);    xp12(:,2)=x12(:,1);    Pp(:,2)=P(:,1)+Q;
L(:,1)=ones(3*fpts,1); w(:,1)=ones(3*fpts,1);
if (jackknife)
    jk(:,:,1)=zeros(fpts,5,1);      % Zero variance for single point
    %jkmean(:,:,1)=1;
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
    bias =  - psi(0,cl(ind).seg_tot) + log(cl(ind).seg_tot);
    newz=[sp(:,1,ind); sp(:,2,ind)] + bias;        % Bias correction in estimator
    % Construct measurement noise covariance (from updated spectrum)
    R=varstab*psi(1,cl(ind).seg_tot);%*(N/cl(ind).seg_len);
    % Evidence
    evidence(:,ind) = exp(-((newz-xp(:,ind)).^2)./(2*R));
    % Jazwinski algorithm (update Q)
    enewz = (newz-xp(:,ind));
    qj = ( enewz.^2 - (R+P(:,ind-1)) ); qj(qj<0)=0;
    Q = jazalph*Q + (1-jazalph)*qj;
    Pp(:,ind) = P(:,ind-1)+Q;                                       	% Update `a priori' state covar
    
    % Kalman filter auto-spectra
    [x(:,ind),P(:,ind),K,Pp(:,ind+1),xp(:,ind+1)]=kalm_filtrv(newz,Q,R,theta,H,Pp(:,ind),xp(:,ind));
    
    % Form cross-spectra from auto-spectral weights
    newz12=[sp(:,3,ind); sp(:,4,ind)] + bias;
    a=K(1:end/2); b=K((end/2+1):end);                                   % Auto-spectral weights from K (a=ch.1, b=ch.2)
    K12a=(a+b)./2; K12b=((1-a)+(1-b))./2;                               % Cross-spectral gains
    x12(:,ind)=[K12a; K12a].*newz12 + [K12b; K12b].*xp12(:,ind);        % Construct cross-spectral measure
    xp12(:,ind+1)=x12(:,ind);                                           % Define next 'a priori'
    
    % Construct weight matrix to determine effective no. segments
    Ka=[K; K12a]; Kb=[(1-K); K12b];                                     % Appended is cross-spec weighting
    w(:,1:ind)=[Kb(:,ones(1,ind-1)).*w(:,(1:(ind-1))) Ka];
    if ((~smooth) || (ind==M))
        % Force cross-spectral weights to sum=1 (as not kalman filtered)
        w2=w(:,1:ind); sumw12=sum(w(2*end/3+1:end,1:ind),2);
        w2(2*end/3+1:end,1:ind)=w2(2*end/3+1:end,1:ind)./sumw12(:,ones(1,size(w2,2)));
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
	PN=zeros(vlen,M,classstr);
    x12N=zeros(2*fpts,M,classstr);
    % Initial conditions of smoother
    xN(:,M)=x(:,M);
    PN(:,M)=P(:,M);
    x12N(:,M)=x12(:,M);
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
        A=P(:,ind)./Pp(:,ind+1);
        xN(:,ind)=x(:,ind)+A.*(xN(:,ind+1)-xp(:,ind+1));
        PN(:,ind)=P(:,ind)+A.*(PN(:,ind+1)-Pp(:,ind+1)).*A;
        % Construct smoothed cross-spectra from auto-spectral gains
        c=sqrt(A(1:end/2).*A(end/2+1:end));
        x12N(:,ind)=x12(:,ind)+[c; c].*(x12N(:,ind+1)-xp12(:,ind+1));
        
        % Adjust effective no. segments for backward sweep
        Aa=[A; c];                                              % Appended is cross-spec weighting
        Ab=[(1-A); sqrt((1-A(1:end/2)).*(1-A(end/2+1:end)))];   %  |
        % Re-normalise forward weights
        wf=[w(:,1:ind) zeros(3*fpts,M-ind)];
        wsum=sum(wf,2);
        wf=wf./wsum(:,ones(1,M));
        % Normalise conditioned cross-spectral weights
        wb=w2;
        wsum12=sum(wb(2*end/3+1:end,:),2);
        wb(2*end/3+1:end,:)=wb(2*end/3+1:end,:)./wsum12(:,ones(1,size(wb,2)));
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
	x=xN; P=PN; x12=x12N; clear('xN','PN','x12N');
end;
clear('w','w2');

% Formulate frequency domain output
x12=real(exp(x12(1:fpts,:)))+i*real(exp(x12((fpts+1):end,:)));    % (complex) cross-spectrum
sp=reshape([exp(x); x12],fpts,3,M);
params.Psp=reshape(P,fpts,2,M);
params.evidence=evidence;
clear('x','P','xp','Pp','x12');

% Generate parameters structure
params.L=L;
if (jackknife)
    params.jk=jk; clear('jk');
    params.jkmean=jkmean; clear('jkmean');
end;

% Calculate coherence & phase
warning off MATLAB:divideByZero
sp(:,4,:)=abs(sp(:,3,:)).^2./((sp(:,1,:)).*(sp(:,2,:)));
warning on MATLAB:divideByZero
sp(:,5,:)=angle(sp(:,3,:));
