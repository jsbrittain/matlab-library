function [sp,params]=ff_tfspm2w(sp,cl,delta,smooth,jackknife);
%function [sp,params]=ff_tfspm2w(sp,cl,delta,smooth,jackknife);
%
% Exponential-decay (Forgetting-factor) analysis
%
% Input parameters
%       sp              Spectral matrix of Fourier coefficients.
%       cl              Confidence limits structure for Fourier analysis.
%       delta           (opt) Exponential decay constant
%       smooth          (opt) Smooth spectra (1=Yes 0=No(default)).
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
%function [sp,params]=ff_tfspm2w(sp,cl,delta,smooth,jackknife);

% Check input parameters
if (~exist('jackknife'))
    jackknife=logical(0);
end;

% Determine data parameters
fpts=size(sp,1);
M=size(sp,3);
vlen=3*fpts;

% Allocate memory
x=zeros(vlen,M);
xp=zeros(vlen,M+1);

% Initialise
initL=1;
xp(:,1)=reshape(mean(sp(:,:,1:initL),3),vlen,1);

% Perform forgetting-factor recursion
for ind=1:M
    newz=reshape(sp(:,:,ind),vlen,1);
    x(:,ind)=(1-delta)*newz+delta*xp(:,ind);
    xp(:,ind+1)=x(:,ind);
end;

% Reshape for output and generate coherence/phase
sp=reshape(x,fpts,3,M);
sp(:,4,:)=(abs(sp(:,3,:)).^2)./((sp(:,1,:)).*(sp(:,2,:)));
sp(:,5,:)=angle(sp(:,3,:));

% Determine effective no. segments
w=delta.^[0:M-1];
L=zeros(1,M);
for ind=1:M
    L(ind)=1/sum((w/sum(w)).^2);
end;
L=L(ones(vlen,1),:);

% Return empty covariance matrix
Psp=[];

% Form parameters stucture
params.Psp=Psp;
params.L=L;
