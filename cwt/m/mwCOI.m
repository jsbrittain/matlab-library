function wlparam=mwCOI(wlparam)
%function wlparam=mwCOI(wlparam)
%
% Determine COI using wraparound arguments
%
% Input parameter
%   wlparam     Wavelet options
%
% Output parameters
%   wlparam     Wavelet options (includes COI determination)
%
%function wlparam=mwCOI(wlparam)

% Extract variables from wlparam structure
ga=wlparam.gamma;
be=wlparam.beta;
K=wlparam.K;
r=wlparam.r;
C1=wlparam.C1;
C2=wlparam.C2;
N=length(wlparam.coi);

% Determine C based on Kth concentration ratio
Econc=0.99;
y=invbetainc(sqrt(Econc),(K-1)+1,r-1);          % Returns (C-1)/(C+1)
C=(y+1)/(1-y);
wlparam.Econc=Econc;
wlparam.C=C;

% Min and max frequency bounds (Olhede 2002, 2003I)
fminmax=C1./(2*pi*(C+[1 -1]*sqrt(C^2-1)).^(1/ga));
wlparam.fmin=fminmax(1);
wlparam.fmax=fminmax(2);

% Min and max time bounds (Olhede 2003I)
fstar=(1/pi)*2^(1/ga-1)*C1*(C*(2-ga)+sqrt(4*(ga-1)+C^2*(ga-2)^2))^(-1/ga);
tminmax=[-1 1]*C2*(C1/(2*pi*fstar))^(1-ga)*sqrt(2*C*(C1/(2*pi*fstar))^ga-(C1/(2*pi*fstar))^(2*ga)-1);
wlparam.tmin=tminmax(1);
wlparam.tmax=tminmax(2);

% Determine fupper numerically
Econc=0.995;        % Add half of remainder to account for unidirectional cumulant
ff=0:0.001:1;
psi=zeros(1,length(ff));
for k=0:K-1
    psi=psi+wlf_morse(2*pi*ff,1,wlparam.dt,be,ga,k,K);
end;
psi=psi/K;
cumwav=cumsum( abs(psi).^2/sum(abs(psi).^2) );
fupper=ff(find(abs(cumwav-Econc)==min(abs(cumwav-Econc))));
wlparam.fupper=fupper;

% Min/max plotting scales
wlparam.plotlimit.minscale=2*wlparam.dt*fupper;
wlparam.plotlimit.maxscale=N*wlparam.dt/(wlparam.tmax-wlparam.tmin);
if (isfield(wlparam,'freqs'))
    wlparam.plotlimit.maxfreq=1./(wlparam.fourier_factor*wlparam.plotlimit.minscale);
    wlparam.plotlimit.minfreq=1./(wlparam.fourier_factor*wlparam.plotlimit.maxscale);
end;

% Wraparound
tau=wlparam.scale*(wlparam.tmax-wlparam.tmin)/2/wlparam.dt;
wlparam.tau=tau;
wlparam.coi=zeros(1,ceil(N/2));
for ind=2:ceil(N/2)
    tmp=find(tau>ind);
    if isempty(tmp)
        wlparam.coi(ind)=Inf;
    else
        wlparam.coi(ind)=wlparam.scale(tmp(1));
    end;
end;
wlparam.display_coi=1;
wlparam.coi=[wlparam.coi fliplr(wlparam.coi(ceil(N/2)-floor(N/2)+1:end))];
wlparam.coi([1 end])=0;                                                     % Scale COI
if (isfield(wlparam,'freqs'))
    wlparam.coi=[Inf 1./(wlparam.fourier_factor*wlparam.coi(2:end-1)) Inf]; % Freq COI
end;
