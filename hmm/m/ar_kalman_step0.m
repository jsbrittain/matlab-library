function [q,s2,evidence,ct,rhoavg] = ar_kalman_step0(y,p,s2static,alpha)
%
% Stepwise Kalman AR-filter - Initial step
%
% Negative `alpha' rate is interpreted as the (negative) static noise
% covariance matrix ( -W ; can be scalar or matrix form ).
%

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
    % Matlab routines
    disp('Estimating measurement noise variance ...');
    T = 100*p;%2*p;
    npts = (1:round(T/2):(N-T));
    ntt = (1:T);
    s2static = zeros(1,length(npts));
    for n = (1:length(npts))
        switch (1)
            case 0,
                m = arx(y(npts(n)+ntt),p);
                s2static(n) = m.NoiseVariance;
            case 1,
                addpath('~/matlab/library/thirdparty/arfit/arfit');
                [w,A,s2static(n),SBC,FPE,th]=arfit(y(npts(n)+ntt)',p,p);
        end;
    end;
    s2static = sort(s2static);
    s2static = mean(s2static(max([1 round(length(s2static)/4)]):round(3*length(s2static)/4)));
end;

% Reserve memory
s2=zeros(1,N);
yhat=zeros(1,N);
ct=zeros(1,N);
s2y=zeros(1,N);
K=zeros(p,N);
q=zeros(p,N);
S=zeros(p,p,N);
rho=zeros(p,p,N);
rhoavg=zeros(1,N);
evidence=zeros(1,N);

% Initial estimates
q(:,p+1) = inv(y(1:p)*y(1:p)')*y(1:p)*y(p+1);                       % Ordinary-least-squares regression (1:p) onto (p+1)
Ft = -(y(p:(-1):1));
S(:,:,p+1) = s2static*Ft*(Ft.') * eye(p);
ct = zeros(1,N);

% Recursively estimate AR parameters
percent=0;
for t=((p+2):N)
    
    % Display progress
    if (round(t/N*100)>percent)
        percent=round(t/N*100);
        disp(['Kalman fiter ' num2str(percent) '%'])
    end;
    
    % Measurement noise (static)
    s2(t)=s2static;                                                 % Measurement/observation noise estimate
    
    % Kalman step
    [yhat(t),evidence(t),ct(t),s2(t),s2y(t),K(:,t),q(:,t),S(:,:,t),rho(:,:,t),rhoavg(t)] = ar_kalman_step(y((t-p):t),p,alpha,ct(t-1),s2(t),q(:,t-1),S(:,:,t-1));
    
end;

% Calculate evidence
%warning off % div 0 on first (p+1) points
%evidence = exp(-((y-yhat).^2)./(2*s2y)) ./ sqrt(2*pi*s2y);          % PDF at y given N(yhat,s2y)
%warning on
