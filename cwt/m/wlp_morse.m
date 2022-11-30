function [param]=wlp_morse(arg1,arg2,arg3,arg4,b,g,k,dk);
%function [param]=wlp_morse(beta,gamma,k,[dk]);
%function [param]=wlp_morse(omega,scale,dt,dj,beta,gamma,k,[dk]);
%
% Return wavlet parameters
%
% Parameter returns may be reduced or complete dependent upon inputs
%
% Output parameters (Forms the structure param):
%   fourier_factor  Fourier factor (scale->Fourier period)
%   coi             e-folding time (tau_s)
%   Cg              Admissibility constant
%

% Determine input parameters
reduced=logical(0);
if (nargin<7)
    reduced=logical(1);  % Return a reduced set of parameters
    b=arg1;
    g=arg2;
    k=arg3;
    if (nargin==4)
        dk=arg4;
        K=length(dk);
    end;
else
    omega=arg1;
    scale=arg2;
    dt=arg3;
    dj=arg4;
end;

% Definition parameters
param.mother='morse';
param.beta=b;
param.gamma=g;
param.k=k;
param.paramstr={'\beta','\gamma','A'};

% Determine fourier-factor numerically (as Olhede & Walden)
omega=2*pi*(0:0.01:100);
WF=zeros(1,length(omega));
if (exist('dk')==1)
    K=length(dk);
    Krange=0:K-1;
    param.wlopt={b,g,k,dk};
    param.K=K;
else
    Krange=k;
    dk=1;
    param.wlopt={b,g,k};
end;
wk=dk./sum(dk);
for ind=1:length(Krange)
   WF=WF+wk(ind)*(wlf_morse(omega,1,1,b,g,Krange(ind)).^2);
end;
WF=WF.*omega/2/pi;
omega0=mean(omega(WF==max(WF)));
f0=omega0/2/pi;
param.fourier_factor=1/f0;
if ((omega0==omega(1)) | (omega0==omega(end)))
    warning(' f0 beyond default frequency range.');
end;

% Generalise Morse parameters
r=(2*b+1)/g;
C1=2^(-1/g)*gamma(r+1/g)/gamma(r);
C2=b*g^(-1)*2^(1/g)*gamma(r-1/g)/gamma(r);
param.r=r;
param.C1=C1;
param.C2=C2;

% Reduced calculated wavelet parameters
param.coi=Inf; %param.coi=1/sqrt(2);
param.display_coi=logical(0);

if (~reduced)
    %param.Cg=wlCg(omega,scale,dt,dj,param.mother,param.wlopt);     %%% DOESNT WORK IN THIS CASE %%%
    param.Cg=1;         %%% UNKNOWN %%%
end;
