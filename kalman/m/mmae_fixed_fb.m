function [condx, condP]=mmae_fixed_fb(dat, Q, R, theta, H, filters);
%function [condx, condP]=mmae_fixed_fb(dat, Q, R, theta, H, filters);
%
% Multi-Modal Adaptive Estimator (MMAE)
% Implementation: Fixed Model ; Multivariate ; Forward-Backward smoothing
%
% By default if 'filters' is specified but Q,R,theta,H are all two dimensional
%  (i.e.: only one matrix for each specified) then a bank of filters will be
%  generated with an evenly distributed process noise, distributed as multiples
%  of the original Q matrix.
%
% Only variable Q,R matrices supported at this time.
%
% Outputs
%   condx       Conditioned estimate of state x.
%   condP       Conditioned error covariance matrix.
%
% Parameters
%   dat         Data matrix (Rows=channels; Cols=time).
%   Q           Process noise covariance.
%   R           Measurement noise covariance.
%   theta       State transition matrix.
%   H           State-to-measurement transition matrix.
%   filters     No. of filters in filter bank.
%
% function [condx, condP]=mmae_fixed_fb(dat, Q, R, theta, filters);

% Determine input data dimensions
dat=dat';
len=length(dat);
channels=size(dat,1);
if (~exist('filters'))
    filters=max(size(Q,3),size(R,3));
end;

% Generate range of noise variance statistics as required
Q2=zeros(channels,channels,filters);
R2=zeros(channels,channels,filters);
if ((size(Q,3)~=filters) & (size(R,3)~=filters))
    for j=1:filters
        Q2(:,:,j)=Q.*(2^(j-(filters+1)/2));
    end;
else
    Q2=Q;
end;
if (size(R,3)~=filters)
    for j=1:filters
        R2(:,:,j)=R;
    end;
end;

% Allocate memory
x=zeros(channels,filters);
P=zeros(channels,channels,filters);
f=zeros(channels,channels,filters);
t_x=zeros(channels,len,filters);
t_P=zeros(channels,channels,len,filters);
condx=zeros(channels,len);
condP=zeros(channels,channels,len);

% Model conditioned estimate (Welch MMAE (c)) - includes 'a priori' estimates for the next time loop
for j=1:filters
    tic
    [t_x(:,:,j), t_P(:,:,:,j)] = kalm_fb(dat, Q2(:,:,j), R2(:,:,j), theta, H);
    disp(['FB complete for filter ' int2str(j) ' (' num2str(toc) ' secs)']);
end;

% Determine filter probability
for j=1:filters
    p(:,:,j)=(1/filters)*eye(channels);
end;

tic;
for t=1:len
	
    % Prepare for recusion
    x=reshape(t_x(:,t,:), channels, filters);
    P=reshape(t_P(:,:,t,:), channels, channels, filters);
    xp=x;
    Pp=P;
    newz=dat(:,t);
    
	% Calculate conditional probability (Welch MMAE (a))
	for chan=1:channels
		for j=1:filters
			C=H(chan,chan)*Pp(chan,chan,j)*H(chan,chan)' + R2(chan,chan,j);
			f(chan,chan,j)=mmae_dens(newz(chan), H(chan,chan)*xp(chan,j), C);
		end;
	end;
	
	% Compute model-conditioned state estimate (Welch MMAE (d))
	condx(:,t)=0;
	for j=1:filters
        condx(:,t) = condx(:,t) + p(:,:,j)*x(:,j);
	end;
	
	% Compute model-conditioned error covariance matrix(Welch MMAE (e))
	condP(:,:,t)=0;
	for j=1:filters
        e=condx(:,t) - x(:,j);
        condP(:,:,t) = condP(:,:,t) + p(:,:,j)*(P(:,:,j) + e*e');
	end;
    
    % Display progress
    if mod(t,10000)==0
        disp(['Progress: ' int2str(t) ' of ' int2str(len) ' (' int2str(100*t/len) '%) (' num2str(toc) ' secs)']);
        tic
    end;
end;
