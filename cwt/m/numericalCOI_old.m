function wlparam=numericalCOI(wlparam,N);
%
% Calculation of a numerical Cone-of-Influence
%
% Designed to be run after either a wl(f)Transform or a mw(f)Transform.
%

% Display error
error(' Numerical COI evaluation is inefficient. Consider ht_coi, mw_coi instead.');

% Check input arguments
if (nargin~=2)
    error(' Incorrect number of input arguments');
end;

% Check if dealing with a multiwavelet
multiwavelet=logical(0);
if (isfield(wlparam,'multiwavelet'))   % Check if transform is wavelet or multi-wavelet
	if (wlparam.multiwavelet)
        multiwavelet=1;
        wlopt=wlparam.wlopt;
        K=wlopt{end};
	end;
end;

% Determine scale based COI region
e2=(exp(-2))^(1/4);
wlparam.coi=[];
for ind=1:length(wlparam.scale)
    disp(['Determining numerical COI: scale ' int2str(ind) ' of ' int2str(length(wlparam.scale))]);
    N2=2^ceil(log2(N)+1);       % Fine resolution frequency grid
    halflength=floor(N2/2);
    omega=[(2*pi*(0:N2-halflength-1)/(N2*wlparam.dt)) (-(2*pi*(0:halflength-1))/(N2*wlparam.dt))]';
    if (multiwavelet)
        % Take mean of freq. domain multiwavelets up to order K
        X=zeros(N2,1);
        for k=0:K-1
            %if (isfield(wlparam,'K'))
            %    wlopt={wlparam.wlopt{1:end-1}};
            %else
                wlopt=wlparam.wlopt;
            %end;
            wlopt{end}=k;
            X=X+feval(['wlf_' wlparam.mother],omega,wlparam.scale(ind),wlparam.dt,wlopt{:});
        end;
        X=X./K;     % Take mean
    else
        % Return freq. domain wavelet
        X=feval(['wlf_' wlparam.mother],omega,wlparam.scale(ind),wlparam.dt,wlparam.wlopt{:});
    end;
    R=abs(ifft(abs(X).^2));                                                 % Autocorrelation sequence
    R=R/R(1);       % Normalise sequence (=1 at lag 0)    
    tau(ind)=mean(find(abs(R(1:end/2)-e2)==min(abs(R(1:end/2)-e2))));       % Find e-folding time
end;
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
