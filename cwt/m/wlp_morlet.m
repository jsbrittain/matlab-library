function [param]=wlp_morlet(arg1,scale,dt,dj,w0);
%function [param]=wlp_morlet(w0);
%function [param]=wlp_morlet(omega,scale,dt,dj,w0);
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
    w0=arg1;
else
    omega=arg1;
end;

% Definition parameters
param.mother='morlet';
param.w0=w0;
param.wlopt={w0};

% Reduced calculated wavelet parameters
param.fourier_factor=(4*pi)/(w0+sqrt(2+w0.^2));
param.coi=1/sqrt(2);

% Full parameters
if (~reduced)
    
    % Calculate Cg
	switch (w0)         % Pre-calculated common w0 values
        case 6
            param.Cg=0.776;
        otherwise
            param.Cg=wlCg(omega,scale,dt,dj,param.mother,{w0});
    end;
    
end;
