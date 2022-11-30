function Cg=wlCg(omega,scale,dt,dj,mother,wlopt);
%function Cg=wlCg(omega,scale,dt,dj,mother,wlopt);
%
% Determine Cg for arbitrary wavelet
%
% NB: This calculation involves vectors for omega and scale,
%     thus the resultant Cg is scale dependent and may be
%     inconsistent to some degree between calculations using
%     different scale/frequency resolutions
%
% Input parameters
%   omega       Frequency range
%   scale       Scale vector
%   dt          Time resolution
%   dj          Scale resolution
%   mother      Mother wavelet
%   wlOpt       Wavelet options
%
% Cg calculation performed as in Torrence et al. (1998) section 2i (Eq 13)
% (American Meteorological Society 79(1):61-78)
%

% Determine mother wavelet
wl_mother=['wl_' mother];
wlf_mother=['wlf_' mother];

% Calculate wavelet transform of a delta function at time n=0 over all scales
Wdelta=zeros(length(scale),1);
for sind=1:length(scale)
    % Compute normalised wavelet transform at scale s
    wlf=feval(wlf_mother,omega,scale(sind),dt,wlopt{:});
    Wdelta(sind)=mean(conj(wlf));
    
    % Increment scale index
    sind=sind+1;
end;

% Calculate wavelet in time domain at location 0
psi0=feval(wl_mother,0,wlopt{:});

% Compute Cg (from Torrence eq. 13)
Cg=(dj*sqrt(dt)/psi0)*sum(real(Wdelta)./sqrt(scale'));
