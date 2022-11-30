function params=mt_bootstrap(p11,p22,p12,frange,params,B);
%function params=mt_bootstrap(p11,p22,p12,frange,params,[B]);
%
% Multitaper bootstrap routine for use with mt_sp2.m
% Bootstraps spectra
%
% Input parameters
%       p11         Matrix of single trial spectral estimators (ch.1)
%       p22          | ch.2
%       p12          | ch.1,2
%       frange      2-element vector of frequency offsets for
%                    cross-spectral truncation [fmin fmax]
%       params      Parameters structure
%       B           Resample count
%
%function params=mt_bootstrap(p11,p22,p12,frange,params,[B]);

% Check input parameters
if ((nargin<5) | (nargin>6))
    error(' Incorrect number of input parameters');
end;
% Resample count
if (exist('B')), params.B=B;
else             params.B=1000; end;
% Frequency truncation range
fmin=frange(1);
fmax=frange(2);

% Set parameters switch
params.bootstrap=logical(1);

% Determine data parameters
M=size(p11,2);
fcount=size(p11,1);
qfcount=size(p12,1);

% Reserve variable space
params.bs11=zeros(fcount,B);
params.bs22=zeros(fcount,B);
params.bs12=zeros(fcount,B);
params.bsq=zeros(qfcount,B);

% Display status
disp(['Bootstrapping... (resample count ' int2str(B) ')']);

% Determine if Neurospec spectral norm used (for cumulant)
if (~isfield(params,'mtparams')), mtparams=params;
else                              mtparams=params.mtparams; end;
q_norm=1; dur=params.padcount*params.rate/1000;    % Duration samples
if (mtparams.spec_norm/(params.duration*params.rate/1000)==2*pi)
    q_norm=2*pi;
end;

% Bootstrap single trial estimators
for ind=1:params.B
    % Resample with replacement
    resample=round((M-1)*rand(1,M))+1;
    bsp11=p11(:,resample);
    bsp22=p22(:,resample);
    bsp12=p12(:,resample);
    % Calculate cumulant density
	bsq(:,ind)=q_norm*real(fftshift(ifft(mean(bsp12,2))));
    % Truncate cross-spectra
    bsp12=bsp12(fmin:fmax,:);
    % Calculate mean spectra
    bs11(:,ind)=mean(bsp11,2); bs22(:,ind)=mean(bsp22,2); bs12(:,ind)=mean(bsp12,2);
end;

% Calculate coherence
bscoh=atanh(abs(bs12)./sqrt(bs11.*bs22));
% Calculate phase
bsph=angle(bs12);

% Calculate variance
params.bs11v = transpose(var(transpose(log10(bs11))));
params.bs22v = transpose(var(transpose(log10(bs22))));
params.bs12v = transpose(var(transpose(log10(bs12))));
params.bsqv = transpose(var(transpose(bsq)));
params.bscohv = transpose(var(transpose(bscoh)));
params.bsphv = transpose(var(transpose(bsph)));
