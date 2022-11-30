function [Sx,out2,Sxy,swparam]=swfTransform(varargin);
%function [Sx,[Sy,Sxy],swparam]=swfTransform(x,[y],dt,df,frange,opt,wlopt);
%
% Slepian-Wavelet Transform
%
% Input parameters
%       x           Data vector channel 1
%       y           (opt) Data vector channel 2
%       dt          Sampling interval
%       df          Frequency analysis resolution 1/(scales/octave)
%       frange      Frequency range of analysis [fmin,fmax]
%       opt         Data processing options
%                       s - single taper
%       wlopt       {p(time-bandwidth product), pc(time-bandcentre product)}
%
%function [Sx,[Sy,Sxy],swparam]=swfTransform(x,[y],dt,df,frange,opt,wlopt);

% Default parameters
bivariate=logical(0);

% Check number of input time series
if ~(length(varargin{1})>1)                     % Length of 'x'
    error(' Time series too short (Must contain multiple elements).');
end;
if (length(varargin{2})==length(varargin{1}))   % Two time series (len 'y' == len 'x')
    bivariate=logical(1);
end;

% Determine input parameters
x=varargin{1};
if (bivariate)
    y=varargin{2};
    ind=3;
else
    ind=2;
end;
dt=varargin{ind}; ind=ind+1;
df=varargin{ind}; ind=ind+1;
frange=varargin{ind}; ind=ind+1;
opt=varargin{ind}; ind=ind+1;
wlopt=varargin{ind};

% Default parameters
singletaper = false;

% Process data
opt_str = opt;
options=deblank(opt_str); opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 's'                    % Single taper only
        	singletaper = true;
        otherwise                   % Options for wavelet analysis
            error(['Unknown option -- ' opt]);
    end;
end;

% Extract parameters
p=wlopt{1}; pc=wlopt{2};
fmin=frange(1);
fmax=frange(2);
rate=1/dt;

% Determine analysis options
maxk=ceil(4*p);
J0=log2(fmin);
J1=log2(fmax);
freqs=2.^[J0:df:J1];  % fc

% Determine maximum available analysis frequencies
M=round(pc*rate./freqs);
freqs=freqs(M+1-maxk>0); M=M(M+1-maxk>0);     % Produce nwin=0

% Allocate space
Sx=zeros(length(freqs),length(x));
if (bivariate)
    Sy=zeros(size(Sx));
    Sxy=zeros(size(Sx));
end;

% Perform slepian wavelet analysis per frequency
for ind=1:length(freqs)
    disp(['Slepian wavelet freq ' int2str(ind) ' of ' int2str(length(freqs))]);
    
    % Generate (complex) slepian wavelets
    [psi0,lambda]=slepwave(M(ind),maxk,p,pc);
    psi=(psi0(:,1:2:end)+1i*psi0(:,2:2:end))/sqrt(2);
    if (singletaper)
        psi = psi(:,1);
    end;
    
    % Apply weighting (uniform at present)
    tapercount=size(psi,2);
    w=ones(tapercount,1)/tapercount;
    
    % Multitaper analysis
    for k=1:tapercount
        W1=transpose(conv(x,psi(:,k)));
        trange=floor((length(W1)-length(x))/2+1):length(W1)-ceil((length(W1)-length(x))/2);
        W1=W1(trange);
        Sx(ind,:)=Sx(ind,:)+w(k)*(abs(W1).^2);
        if (bivariate)
            W2=transpose(conv(y,psi(:,k)));
            W2=W2(trange);
            Sy(ind,:)=Sy(ind,:)+w(k)*(abs(W2).^2);
            Sxy(ind,:)=Sxy(ind,:)+w(k)*(W1.*conj(W2));
        end;
    end;
end;

% Construct parameters structure (starting with input parameters)
swparam.dt=dt;
swparam.df=df;
swparam.frange=frange;
swparam.opt=opt;
swparam.wlopt=wlopt;
swparam.mother='slepian';
swparam.paramstr={'p','p_c'};
% Additional parameters
swparam.freqs=freqs;
% Compatibility parameters (largely for plot routines)
swparam.linearscale=logical(0);
swparam.display_coi=logical(0);
swparam.tapercount=tapercount;
swparam.L=tapercount;
swparam.multiwavelet=logical(1);
swparam.slepian=logical(1);
swparam.coi=[];
swparam.fourier_factor=1;
swparam.scale=1./swparam.freqs;

% Define output argments
if (bivariate)
    out2=Sy;
else
    out2=swparam;
    clear('swparam');
end;
