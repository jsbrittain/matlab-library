function [param]=wlp_dog(arg1,scale,dt,dj,m);
%function [param]=wlp_dog(m);
%function [param]=wlp_dog(omega,scale,dt,dj,m);
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
reduced=0;
if (nargin<5)
    reduced=1;  % Return a reduced set of parameters
    m=arg1;
else
    omega=arg1;
end;

% Definition parameters
param.mother='dog';
param.m=m;
param.wlopt={m};

% Reduced calculated wavelet parameters
param.fourier_factor=(2*pi)/sqrt(m+.5);
param.coi=1/sqrt(2);

% Full parameters
if (~reduced)
    
    % Calculate Cg
	switch (m)         % Pre-calculated common w0 values
        case 2
            param.Cg=3.541;
        otherwise
            param.Cg=wlCg(omega,scale,dt,dj,param.mother,{m});
    end;
    
end;
