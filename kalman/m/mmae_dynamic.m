function [condx,condP]=mmae_dynamic(dat, Q, R, theta, H, filters);
%function [condx, condP]=mmae_dynamic(dat, Q, R, theta, H, filters);
%
% Multi-Modal Adaptive Estimator (MMAE)
% Implementation: Dynamic Model ; Multivariate
%
% Filters differ in their process noise parameter (chosen to provide
% an even spread around the given matrices)
%
% Outputs
%   condx       Conditioned estimate of state x.
%   condP       Conditioned error covariance matrix.
%
% Parameters
%   dat         Data matrix (cols=channels; rows=time).
%   Q           Process noise covariance (Constant).
%   R           Measurement noise covariance (Constant).
%   theta       State transition matrix (Constant).
%   H           State-to-measurement transition matrix (Constant).
%   filters     No. of filter to implement in the MMAE algorithm.
%
% function [condx, condP]=mmae_dynamic(dat, Q, R, theta, H, filters);

% Check input parameters
if (nargin<6)
    error(' Not enough parameters');
end;
if (nargin>6)
    error(' Too many parameters');
end;

% Determine data parameters
len=size(dat,1);
channels=size(dat,2);

% Allocate memory
xp=zeros(channels,filters);
Pp=zeros(channels,channels,filters);
Q2=zeros(channels,channels,filters);
R2=zeros(channels,channels,filters);
condx=zeros(channels,len);
condP=zeros(channels,channels,len);

% Pre-calculate variances
I=eye(channels);
theta=I; H=I;
for j=1:filters
	xp(:,j)=dat(1,:)';	        % Set xp to first data value
	Pp(:,:,j)=I;	            % Set Pp to identity                                %%%%%%%%%%%
	
    % Create an even spead of process noise statistics around the provided values
    Q2(:,:,j)=Q.*(2^(j-(filters+1)/2));
    R2(:,:,j)=R;
end;

% Recurse MMAE-dynamic algorithm
for t=1:len
	% Get next measurement
	newz=dat(t,:)';

    % Perform recursive Dynamic-MMAE estimation
    [condx(:,t), condP(:,:,t), Pp, xp, p]=mmae_dynamicr(newz, Q2, R2, theta, H, Pp, xp, channels, filters);
    
    % Report on progress
    if mod(t,1000)==0
        disp(['t=' int2str(t) ' of ' int2str(len) ' (' int2str(round(t*100/len)) '%)']);
    end;
end;
