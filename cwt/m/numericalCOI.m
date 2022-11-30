function wlparam=numericalCOI(wlparam,N);
%function wlparam=numericalCOI(wlparam,N);
%
% Determine COI numerically based on description given in Torrence and
% Compo, 1998.
%
% This method generates a wavelet (or multiwavelet) transform with the same
% parameters as the original analysis, using all zero data with
% discontinuities at n=0,N-1.  The e-folding time is then taken at the
% point where the wavelet power drops by a factor e^(-2).
%
% Input parameters
%   wlparam     Wavelet parameters
%   N           Number of data points
%
% This method is slow but computes the numerical COI using the original
% definition and so should be correct.
%
%function wlparam=numericalCOI(wlparam,N);

% Last edited 30-01-2006 (JSB)

% Display error
%error(' Numerical COI evaluation is inefficient. Consider ht_coi, wl_coi or mw_coi instead.');

% Display progress
disp('Generating numerical COI');

% Generate data of length N with discontinuity at n=0,N-1
if (nargin<2)
    N=length(wlparam.coi);
end;
x=zeros(N,1);
x([1 end])=1;
time=[0:N-1]*wlparam.dt;

% Check for multiwavelet identifier
if (~isfield(wlparam,'multiwavelet'))
    wlparam.multiwavelet = false;
end;

% Perform wavelet transform with original arguments
if (wlparam.multiwavelet)
    if (isfield(wlparam,'slepian'))
        [S,wlparam0]=swfTransform(x,wlparam.dt,wlparam.df,wlparam.frange,wlparam.opt,wlparam.wlopt);
    elseif (isfield(wlparam,'freqs'))
        [S,wlparam0]=mwfTransform(x,wlparam.dt,wlparam.df,wlparam.frange,wlparam.opt,wlparam.mother,wlparam.wlopt);
	else
        [S,wlparam0]=mwTransform(x,wlparam.dt,wlparam.dj,wlparam.srange,wlparam.opt,wlparam.mother,wlparam.wlopt);
	end;
    % Remove wraparound plot points if they exist
    if (isfield(wlparam0,'plotlimit'))
        wlparam=rmfield(wlparam0,'plotlimit');
    end;
else
	if (isfield(wlparam,'freqs'))
        [W,wlparam0]=wlfTransform(x,wlparam.dt,wlparam.df,wlparam.frange,wlparam.opt,wlparam.mother,wlparam.wlopt);
	else
        [W,wlparam0]=wlTransform(x,wlparam.dt,wlparam.dj,wlparam.srange,wlparam.opt,wlparam.mother,wlparam.wlopt);
	end;
    S=abs(W).^2;
end;

% Generate power spectrum and normalise at each scale
maxS=max(S')';
S=S./maxS(:,ones(1,size(S,2)));

% Calculate e-folding time for the discontinuity
wlparam.coi=zeros(1,length(wlparam.scale));
for ind=1:N
    transitionpts=find(diff((S(:,ind)>=exp(-2))));
    if (isempty(transitionpts))
        wlparam.coi(ind)=0;
    else
        wlparam.coi(ind)=wlparam.scale(transitionpts(1));
    end;
end;
wlparam.display_coi=1;

% Convert to frequency if required
if (isfield(wlparam,'freqs'))
    warning off MATLAB:divideByZero
    wlparam.coi=1./(wlparam.fourier_factor*wlparam.coi);
    warning on MATLAB:divideByZero
    wlparam.coi(isinf(wlparam.coi)) = wlparam.freqs(1);
end;

% Display progress
disp('done.');
