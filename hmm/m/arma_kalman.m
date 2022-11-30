function [ar,ma,s2,evidence,ct,rhoavg] = arma_kalman(y,p,q,s2static,alpha)
%function [ar,ma,s2,evidence,ct,rhoavg] = arma_kalman(y,p,q,s2static,alpha)
%
%
%
% Negative `alpha' rate is interpreted as the (negative) static noise
% covariance matrix ( -W ; can be scalar or matrix form ).
%
%function [ar,ma,s2,evidence,ct,rhoavg] = arma_kalman(y,p,q,s2static,alpha)
%
% EXPERIMENTAL

% Check input parameters
if (~exist('s2static','var'))
    s2static = [];
end;
if (~exist('alpha','var'))
    alpha = 0.1;
end;
if (isvector(y))
    if (size(y,1)>1)
        y=y.';
    end;
end;

% Determine data parameters
N=length(y);

% Determine measurement noise (static moving window estimates)
if (isempty(s2static))
    disp('Estimating measurement noise variance ...');
    T = 3*p;
    npts = (1:round(T/2):(N-T));
    ntt = (1:T);
    s2static = zeros(1,length(npts));
    for n = (1:length(npts))
        m = arx(y(npts(n)+ntt),p);
        s2static(n) = m.NoiseVariance;
    end;
    s2static = sort(s2static);
    s2static = mean(s2static(max([1 round(length(s2static)/4)]):round(3*length(s2static)/4)));
end;

pq=p+q;

% Reserve memory
s2=zeros(1,N);
yhat=zeros(1,N);
ct=zeros(1,N);
s2y=zeros(1,N);
K=zeros(pq,N);
state=zeros(pq,N);
S=zeros(pq,pq,N);
rho=zeros(pq,pq,N);
rhoavg=zeros(1,N);
evidence=zeros(1,N);

% Initial estimates  (NEED TO ADJUST FOR CORRECT AR + MA INITIALISATION)
t0 = max(p+1,q);
yp = [y(1:p) sqrt(s2static)*randn(1,q)];                    %%% RANDOM -MA INITIALISATION!!!! %%%
state(:,t0) = inv(yp*yp')*yp*y(t0);                       % Ordinary-least-squares regression (1:p) onto (p+1)
Ft = -(y(p:(-1):1));
S(:,:,p+1) = s2static*Ft*(Ft.') * eye(pq);
ct = zeros(1,N);

% Recursively estimate AR parameters
percent=0;
for t=((t0+1):N)
    
    % Display progress
    if (round(t/N*100)>percent)
        percent=round(t/N*100);
        disp(['Kalman fiter ' num2str(percent) '%'])
    end;
    
    % Measurement noise (static)
    s2(t)=s2static;                                                 % Measurement/observation noise estimate
    
    % (p-dimension) time-series vector
    Ft = [ -(y((t-1):(-1):(t-p))) (y((t-1):(-1):(t-q))-yhat((t-1):(-1):(t-q))) ];	% Regression onto time-series
    yhat(t)=Ft*state(:,t-1);                                        % `a prior' state estimate
    ert=(y(t)-yhat(t));                                             % State estimate error
    
    % Evidence (before Jazwinski)
    Wt = ct(t-1) * eye(pq);
    s2ev = s2(t) + Ft*S(:,:,t-1)*(Ft.') + Ft*Wt*(Ft.');             % Estimated prediction variance
    evidence(t) = exp(-((y(t)-yhat(t)).^2)./(2*s2ev));              % PDF at y given N(yhat,s2y)
    
    % Update state noise covariance through Jazwinski's algorithm
    if ( alpha < 0 )
        % Static state noise covariance
        ct(t) = -alpha;
        Wt = ct(t) * eye(pq);
    else
        % Jazwinski algorithm
        s2state0 = s2(t) + Ft*S(:,:,t-1)*(Ft.');                    % Estimated prediction variance (assumes q=0)
        ct(t) = (ert^2-s2state0)/(Ft*(Ft.')); if (ct(t)<0), ct(t)=0; end;   % Learning rate
        ct(t) = alpha*ct(t-1)+(1-alpha)*ct(t);                      % Smoothed learning rate
        Wt = ct(t) * eye(pq);                                       % Isotropic (c.I) matrix
    end;
    
    % Kalman filter equations
    s2y(t) = s2(t) + Ft*S(:,:,t-1)*(Ft.') + Ft*Wt*(Ft.');           % Estimated prediction variance
    K(:,t) = ((S(:,:,t-1) + Wt)*(Ft.'))/s2y(t);                     % Kalman gain
    state(:,t) = state(:,t-1) + K(:,t)*ert;                         % State mean (AR coefficients)
    S(:,:,t) = S(:,:,t-1) + Wt - K(:,t)*Ft*(S(:,:,t-1) + Wt);       % State covariance
    
    % Effective learning rate
    SnI = S(:,:,t-1) + ct(t)*eye(pq);
    rho(:,:,t) = SnI/s2y(t);
    rhoavg(t) = trace(SnI)/s2y(t)/pq;
    
end;

% Calculate evidence
%warning off % div 0 on first (p+1) points
%evidence = exp(-((y-yhat).^2)./(2*s2y)) ./ sqrt(2*pi*s2y);          % PDF at y given N(yhat,s2y)
%warning on

ar = state(1:p,:);
ma = state((p+1):end,:);
