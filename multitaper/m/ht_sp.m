function [sp11,sp22,sp12,params]=ht_sp(dat1,dat2,rate,opt_str)
%function [sp11,[sp22,sp12],params]=ht_sp(dat1,[dat2],rate,opt_str)
%
% Multitaper time-frequency analysis utilising Hermite functions
%
% Implementation
%   Eigenvalue weighting
%
% Input parameters
%       dat1        Time series ch.1
%       dat2        Time series ch.2
%       rate        Sampling rate
%       opt_str     Options string
%                       r<n>        rectify channels (n={0:ch.1,1:ch.2,2:ch1,2})
%                       f<n>        maximum analysis frequency (default: Nyquist)
%                       g<n>        minimum analysis frequency (used for long narrowband exminations)
%                       A<n>        localisation area (default: 12)
%                       t<msecs>    time skip for plotting (default: 50 msecs)
%                       s<n>        time scaling for elliptical discs (default: 1)
%                       n           power normalise each increment
%                       k           fixed taper count (overrides A<n>)
%                       e           return eigenspectra ( central time-point only; index by time-vector in sp<n><m> )
%                       a           re-(a)ssignment spectrogram (onto existing grid)
%                                       [ NOT YET IMPLEMENTED ]
%
%function [sp11,[sp22,sp12],params]=ht_sp(dat1,[dat2],rate,opt_str)

% Elliptical geometry
%   For an ellipse with maximum diameter 2a (max radius a) and minimum
%      diameter 2b (min radius b), its area is (pi*a*b).
%   For circle with time scaling s, elliptical deformation results in
%      max radius a=R/s and min radius b=R*s.  For (s<1) the ellipse
%      therefore expands in time and compresses in frequency.

% Check input parameters
bivariate=true;
if (nargin==3),     opt_str=rate; rate=dat2; clear('dat2');
                    bivariate=false;
elseif (nargin~=4), error(' Incorrect number of input parameters');
end;

% Determine data parameters
N=length(dat1);

% Assign default parameters
minf = 0;                       % minimum analysis freq (default: DC)
maxf=rate/2;                    % maximum analysis freq (default: Nyquist)
A=12;                           % localisation area
skip=50;                        % time skip (msecs)
scale=1;                        % time scaling (forms elliptical shape)
normalise=false;                % power normalise
fixedtapercount = [];           % Fixed taper count
return_eigenspectra = false;    % Return all eigenspectra
reassignment = false;           % TF-reassignment

% Parse options string
options=deblank(opt_str); opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'r'                    % Rectify channels
        	n=str2num(optarg);
			if ((n<0) || (n>2))
			    error(['Error in option argument -- r' optarg]);
			end;  
			if (n~=1)   % Rectify ch 1.
			    dat1=abs(dat1);
			end;  
			if (n>=1)   % Rectify ch 2.
			    dat2=abs(dat2);
			end;
        case 'f'        % Maximum frequency
            maxf=str2num(optarg);
            if (isempty(maxf))
                error(['Error in option argument -- f' optarg]);
            end;
        case 'g'        % Minimum frequency
            minf=str2num(optarg);
            if (isempty(minf))
                error(['Error in option argument -- g' optarg]);
            end;
        case 'A'        % Specify area
            A=str2num(optarg);
            if (isempty(A))
                error(['Error in option argument -- A' optarg]);
            end;
        case 't'        % Time skip
            skip=str2num(optarg);
            if (isempty(skip))
                error(['Error in option argument -- t' optarg]);
            end;
        case 's'        % Time scaling
            scale=str2num(optarg);
            if (isempty(scale))
                error(['Error in option argument -- s' optarg]);
            end;
        case 'k'        % Use specified number of taper (overrides area A)
            fixedtapercount=str2num(optarg);
            if (isempty(scale))
                error(['Error in option argument -- k' optarg]);
            end;
        case 'n'        % Power normalise
            normalise=true;
        case 'e'        % eigenspectrograms
            return_eigenspectra = true;
        case 'a'        % reassignment spectrogram
            reassignment = true;
            errror(' Reassignment spectrogram not yet implemented!');
        otherwise                   % Options for wavelet analysis
            error(['Unknown option -- ' opt]);
    end;
end;

% Determine analysis parameters
R=sqrt(2*A);                % Disc radius from area (A=(R^2)/2)
skip=skip*rate/1000;        % Time skip msecs -> samples

% Determine no. eigenspectra (K)
K=200; zeta=0.95;
j=[0:K-1]; for k=0:K-1, facj(k+1)=factorial(k); end;
lambda=1-exp(-R^2/2)*cumsum(2.^(-j).*R.^(2*j)./facj);
E=lambda.^2;                % Energy concentration ratio
if (isempty(fixedtapercount))
    K=max(find(E>=zeta));       % Maximum energy conc >= zeta
else
    K=fixedtapercount;
end;
E=E(1:K); lambda=lambda(1:K);
weights=E/sum(E);           % Normalised weights

% % Check options
% if ((normalise) && (K~=1))
%     warning(' Normalisation removed - this should only be used with single tapers.');
%     normalise = false;
% end;

% Form hermite functions
t=[-(N-1):N]';              % Cover complete time axis
t=t/rate/scale;             % Time deformation to ellipical regions
H=zeros(length(t),K);
h=zeros(length(t),K);
facj=[];
for k = (0:K-1)
    % Form hermite polynomials (recursive structure)
    if (k==0),     H(:,1)=ones(length(t),1);
    elseif (k==1), H(:,2)=2*t;
    else           H(:,k+1)=2*t.*H(:,k)-2*(k-1)*H(:,k-1);
    end;
    % Generate hermite functions
    h(:,k+1)=H(:,k+1).*exp(-t.^2/2)/(pi^(1/4)*sqrt(2^k*gamma(k+1)));
end;

% Determine maximum frequency
freqs=rate*[0:N/2]/N;
fmin=max(1,find(abs(freqs-minf)==min(abs(freqs-minf)))); % pts
fmax=min(floor(N/2),find(abs(freqs-maxf)==min(abs(freqs-maxf)),1,'first')); % pts
freqs=freqs(fmin:fmax);

% Allocate memory
t=0:skip:N-1;
if (~return_eigenspectra)
    Tcount = length(t);
else
    Tcount = K;
end;
sp11=zeros(fmax-fmin+1,Tcount);
if (bivariate)
    sp22=zeros(fmax-fmin+1,Tcount);
    sp12=zeros(fmax-fmin+1,Tcount);
end;

% Zero mean data
dat1=dat1-mean(dat1);
if (bivariate)
    dat2=dat2-mean(dat2);
end;

fprintf('Tapers: %d\n<',K);
percent=0;

tind = (1:length(t));
if ( return_eigenspectra )
    tind = round(length(t)/2);
end;

for ind=tind
    if (round(100*ind/length(t))>percent)
        percent = round(100*ind/length(t));
        if ( (mod(percent,10)==0) & (percent<100) )
            fprintf('%d',percent);
        else
            fprintf('.');
        end
    end;
    %disp(['Time location ' int2str(t(ind)*1000/rate) ' of ' int2str(N*1000/rate) 'msecs']);
    for k=(0:K-1)
        % Windowed Fourier transform
        trange=round((1:N)+N-t(ind));                   % Determine Hermite window offset
        F1=fft(h(trange,k+1).*dat1);                    % Windowed FFT
        if (normalise)                                  % Power normalise each windowed segment
            F1=F1/sqrt(sum(abs(F1).^2))/length(F1);     %  | includes correction for fft/ifft conversion
        end;                                            %  |
        F1=F1(fmin:fmax);                               % Reduce to freq range of interest
        if (bivariate)                                  % Repeat for ch.2
            F2=fft(h(trange,k+1).*dat2);                %  | FFT
            if (normalise)                              %  | Power normalise
                F2=F2/sqrt(sum(abs(F2).^2))/length(F2); %  |  |
            end;                                        %  |  |
            F2=F2(fmin:fmax);                           %  | Reduce freq range
        end;
        % Store eigenspectra in time-axixs?
        if ( return_eigenspectra )
            % Replace time index (limited to the central sample) with taper number
            sind = k+1;
        else
            % Usual time-driven analysis
            sind = ind;
        end;
        % Form spectra
        sp11(:,sind)=sp11(:,sind)+weights(k+1)*abs(F1).^2;
        if (bivariate)
            sp22(:,sind)=sp22(:,sind)+weights(k+1)*abs(F2).^2;
            sp12(:,sind)=sp12(:,sind)+weights(k+1)*F1.*conj(F2);
        end;
    end;
end;
fprintf('>\n');

% Assign arguments to parameter structure
params.freqs=freqs;             % Frequency vector
params.L=1/sum(weights.^2);     % Equiv. no. eigenspectra
params.A=A;                     % Area of concentration
params.E=E;                     % Energy contribution from eigenspectra
params.h=h;                     % Hermite functions
% Compatibility parameters for wavelet plot routines
params.mother='hermite';        % Mother (wavelet) function
params.wlopt={scale,A};         % Options {scaling, area}
params.paramstr={'scale','A'};  % Parameters string (corresponds to options)
params.linearscale=true;        % Display on a linear (not log) freq scale
params.display_coi=false;       % Do not display a COI
params.coi=[];                  % Provide empty COI vector
params.rate = rate;

% Truncate output if not bivariate analysis
if (~bivariate)
    sp22=params; clear('params');
end;
