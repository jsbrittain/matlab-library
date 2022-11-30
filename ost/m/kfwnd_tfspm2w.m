function [sp,params]=kfwnd_tfspm2w(sp,cl,K,smooth,jackknife);
%function [sp,params]=kfwnd_tfspm2w(sp,cl,K,smooth,jackknife);
%
% Kalman-Fourier Analysis routine
%
% Sliding window approach to spectral filtering
%
% Input parameters
%       sp              Spectral matrix of Fourier coefficients
%       cl              Confidence limits structure for Fourier analysis
%       K               Number of trials to average over
%       smooth          (opt) Smooth spectra (1=Yes 0=No(default))
%       jackknife       (opt) Jackknife variance (1=Yes 0=No(default))
%
% Output parameters
%       sp              Spectral matrix, cols: (Dim 3 represent group entities)
%                           1       Auto-spectra ch.1
%                           2       Auto-spectra ch.2
%                           3       Cross-spectra ch.1,2
%                           4       Coherence
%       params          Parameters structure
%
%function [sp,params]=kfwnd_tfspm2w(sp,cl,K,smooth,jackknife);

% Determine data parameters
M=size(sp,3);      % Group count
fpts=size(sp,1);   % Number of frequency points
vlen=3*fpts;       % Length of state vector

% Check input parameters
if (~exist('smooth'))
    smooth=logical(0);
end;
if (~exist('jackknife'))
    jackknife=logical(0);
end;

% Initialise variable space
x=zeros(vlen,M-K);
P=zeros(vlen,M-K);
L=K*ones(vlen,M-K);
w=zeros(vlen,M-K);
if (jackknife)
    jk=zeros(fpts,5,M-K);
end;

% Perform sliding window spectral estimation
for ind=1:(M-K)
    x(:,ind)=mean(reshape(sp(:,1:3,ind:(ind+K-1)),vlen,K),2);
end;

% Formulate frequency domain output
sp=reshape(x,fpts,3,M-K);
Psp=reshape(P,fpts,3,M-K);

% Calculate coherence
warning off MATLAB:divideByZero
sp(:,4,:)=(abs(sp(:,3,:)).^2)./((sp(:,1,:)).*(sp(:,2,:)));
warning on MATLAB:divideByZero

% Form parameters stucture
params.Psp=Psp;
params.L=L;
if (jackknife)
    params.jk=jk;
end;
