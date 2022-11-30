function [sp,Psp,L]=kf_tfspm2w(N,sp,cl,Q,smooth);
%function [sp,Psp,L]=kf_tfspm2w(N,sp,cl,Q,smooth);
%
% Hybrid Kalman-Fourier analysis
% Utilises a single linear Kalman filter
%
% EXPERIMENTAL VERSION: CUBE-ROOT TRANSFORM
%
% Remaining issues
%       data length     For trial-varying experimental protocols the trials
%                       will most likely be of different length.
%
% Input parameters
%       N               Length of segments
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
%function [sp,Psp]=kf_tfspm2w(N,sp,cl,Q,smooth);

% Determine data parameters
M=size(sp,3);      % Group count
fpts=size(sp,1);   % Number of frequency points
vlen=2*fpts;       % Length of state vector

% Determine process noise covariance (Q)
if (~exist('Q'))
    Q=0.1/100;     % Default value
end;
if (~exist('smooth'))
    smooth=logical(0);
end;

% Initialise variable space
x=zeros(vlen,M);
x12=zeros(fpts,M);
P=zeros(vlen,M);
xp=zeros(vlen,M+1);
xp12=zeros(fpts,M+1);
Pp=zeros(vlen,M+1);
L=zeros(3*fpts,M);
w=zeros(3*fpts,M);

% Transform spectra to log domain
sp=sp.^(1/3);
%sp(:,1,:)=log10(sp(:,1,:)); sp(:,2,:)=log10(sp(:,2,:));
%sp(:,3,:)=log10(real(sp(:,3,:))) + i*log10(imag(sp(:,3,:)));
varstab=log10(exp(1))^2;

% Transition matrices
theta=1;                                                            % State-transition matrix
H=1;                                                                % Process-to-measurement matrix
I=ones(size(xp,1),1);                                               % Identity matrix (vector in this version)

% Initial conditions
initL=1;                                                            % No. periodograms for initial estimate
initsegs=sum([cl(1:initL).seg_tot]);                                % Segments making up initial estimate
xp(:,1)=reshape(mean(sp(:,1:2,1:initL),3),vlen,1);                  % Initial 'a priori' state vector
xp12(:,1)=reshape(mean(sp(:,3,1:initL),3),fpts,1);                  % Initial 'a priori' cross-spectra (not indexed)
R=varstab*psi(1,initsegs);                                          % Measurement error-covariance matrix
Pp(:,1)=R;                                                          % Initial 'a priori' error cov

% Kalman filter - forward sweep
for ind=1:M
    % Display progress
    disp(['Running kalman filter for group ' int2str(ind) ' of ' int2str(M)]);
    % Construct state vector
    newz=reshape(sp(:,1:2,ind),vlen,1);
    % Construct measurement noise covariance (from updated spectrum)
    R=varstab*psi(1,cl(ind).seg_tot);
    % Kalman filter auto-spectra
    [x(:,ind),P(:,ind),K,Pp(:,ind+1),xp(:,ind+1)]=kalm_filtrv(newz,Q,R,theta,H,Pp(:,ind),xp(:,ind));
    
    %x(1:end/2,ind)=log10((10.^newz(1:end/2).^K(1:end/2)).*(10.^xp(1:end/2,ind).^(1-K(1:end/2))));
    %x(end/2+1:end,ind)=log10((10.^newz(end/2+1:end).^K(end/2+1:end)).*(10.^xp(end/2+1:end,ind).^(1-K(end/2+1:end))));
    
    % Form cross-spectra from auto-spectral weights (see log book 02-05-2006)
    a=K(1:end/2); b=K(end/2+1:end);                                 % Auto-spectral weights from K (a=ch.1, b=ch.2)
    newz12=reshape(sp(:,3,ind),fpts,1);                             % Cross-periodogram (new measurement)
    K12a=sqrt(a.*b); K12b=sqrt((1-a).*(1-b));                       % Current weight 'sqrt(ab)', prev. weight 'sqrt((1-a)(1-b))'
    %K12a=(a+b)/2; K12b=((1-a)+(1-b))/2;
    x12(:,ind)=K12a.*newz12 + K12b.*xp12(:,ind);                    % Construct cross-spectral measure
    
    %x12(:,ind)=(newz12.^K12a).*(xp12(:,ind).^K12b);
    %x12(:,ind)=(real(newz12).^K12a).*(real(xp12(:,ind)).^K12b) + ...
    %           i*(imag(newz12).^K12a).*(imag(xp12(:,ind)).^K12b);
    %sp(:,4,ind)=10.^(2*real(x12(:,ind))-x(1:fpts,ind)-x(fpts+1:2*fpts,ind));
    
    xp12(:,ind+1)=x12(:,ind);                                       % Define next 'a priori'
    
%     % Construct weight matrix to determine effective no. segments
%     Ka=[K; sqrt(K(1:end/2).*K(end/2+1:end))];                       % Appended is cross-spec weighting
%     Kb=[(1-K); sqrt((1-K(1:end/2)).*(1-K(end/2+1:end)))];           %  |
%     w(:,1:ind)=[Kb(1:end*(~isempty(w)),ones(1,ind-1)).*w(:,1:ind-1) Ka];
%     if ((~smooth) | (ind==M))
%         % Force cross-spectral weights to sum=1 (as not kalman filtered)
%         w2=w(:,1:ind); sumw12=sum(w(2*end/3+1:end,1:ind),2);
%         w2(2*end/3+1:end,1:ind)=w2(2*end/3+1:end,1:ind)./sumw12(:,ones(1,size(w2,2)));
%         % Determine effective no. segments (as gains sums to 1 at each step)
%         L(:,ind)=1./sum(w2.^2,2);
%     end;
end;

% Kalman smoother - backward sweep
if (smooth)
    % Allocate memory
	xN=zeros(vlen,M);
	PN=zeros(vlen,M);
    x12N=zeros(fpts,M);
    % Initial conditions of smoother
    xN(:,M)=x(:,M);
    PN(:,M)=P(:,M);
    x12N(:,M)=x12(:,M);
    % Init backward weight matrix
    w2=w;
    % Recurse backward sweep
	for ind=M-1:-1:1
        % Display progress
        disp(['Running kalman filter for group ' int2str(ind) ' of ' int2str(M) ' (B)']);
        % Determine backward gain and combine forward/backward state estimates
        A=P(:,ind)./Pp(:,ind+1);
        xN(:,ind)=x(:,ind)+A.*(xN(:,ind+1)-xp(:,ind+1));
        PN(:,ind)=P(:,ind)+A.*(PN(:,ind+1)-Pp(:,ind+1)).*A;
        % Construct smoothed cross-spectra from auto-spectral gains
        c=sqrt(A(1:end/2).*A(end/2+1:end));
        x12N(:,ind)=x12(:,ind)+c.*(x12N(:,ind+1)-xp12(:,ind+1));
        
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
        L(:,ind)=1./sum((wb.*Aa(:,ones(1,size(wb,2)))+wf.*Ab(:,ones(1,size(wf,2)))).^2,2);
        % Update backward sweep coefficients (un-normalised bwd coeffs)
        w2=w2.*Aa(:,ones(1,size(wb,2)))+wf.*Ab(:,ones(1,size(wf,2)));
    end;
	x=xN; P=PN;
    x12=x12N;
end;
clear('w','w2');

% Formulate frequency domain output
%sp=reshape([10.^x(1:2*fpts,:); 10.^real(x(2*fpts+1:end,:))+i*(10.^imag(x(2*fpts+1:3*fpts,:)))],fpts,3,M);
sp(:,1:3,:)=reshape([x; x12].^3,fpts,3,M);
Psp=reshape(P,fpts,2,M);

% Calculate coherence
warning off MATLAB:divideByZero
sp(:,4,:)=(abs(sp(:,3,:)).^2)./((sp(:,1,:)).*(sp(:,2,:)));
warning on MATLAB:divideByZero
