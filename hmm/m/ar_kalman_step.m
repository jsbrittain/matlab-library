function [yhat,evidence,ct,s2,s2y,K,q,S,rho,rhoavg] = ar_kalman_step(y,p,alpha,ct,s2,q,S)
%
% Stepwise Kalman AR-filter
%
% parameter 'y' should be data((t-p):t)
%
% Function call:
%   [yhat(t),evidence(t),ct(t),s2(t),s2y(t),K(:,t),q(:,t),S(:,:,t),rho(:,:,t),rhoavg(t)] = ar_kalman_step(y((t-p):t)',p,alpharate,ct(t-1),s2(t),q(:,t-1),S(:,:,t-1));
%
% Negative `alpha' rate is interpreted as the (negative) static noise
% covariance matrix ( -W ; can be scalar or matrix form ).
%

% (p-dimension) time-series vector
Ft = -(y((end-1):(-1):(end-p)));                                % Regression onto time-series
yhat = Ft*q;                                                    % `a prior' state estimate
ert = (y(end)-yhat);                                            % State estimate error

% Evidence (before Jazwinski)
Wt = ct * eye(p);
s2ev = s2 + Ft*S*(Ft.') + Ft*Wt*(Ft.');                         % Estimated prediction variance
evidence = exp(-((y(end)-yhat).^2)./(2*s2ev));                  % PDF at y given N(yhat,s2y)

% Update state noise covariance through Jazwinski's algorithm
if (~isempty(alpha))
    if ( alpha < 0 )
        % Static state noise covariance
        ct = -alpha;
        Wt = ct * eye(p);
    else
        % Jazwinski algorithm
        s2q0 = s2 + Ft*S*(Ft.');                                    % Estimated prediction variance (assumes q=0)
        ct0 = ct;
        ct = (ert^2-s2q0)/(Ft*(Ft.')); if (ct<0), ct=0; end;        % Learning rate
        ct = alpha*ct0+(1-alpha)*ct;                                % Smoothed learning rate
        Wt = ct * eye(p);                                           % Isotropic (c.I) matrix
    end;
end;

% Kalman filter equations
s2y = s2 + Ft*S*(Ft.') + Ft*Wt*(Ft.');                          % Estimated prediction variance
K = ((S + Wt)*(Ft.'))/s2y;                                      % Kalman gain
q = q + K*ert;                                                  % State mean (AR coefficients)
S = S + Wt - K*Ft*(S + Wt);                                     % State covariance

% Effective learning rate
SnI = S + ct*eye(p);
rho = SnI/s2y;
rhoavg = trace(SnI)/s2y/p;
