function [sp,params]=kf_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);
%function [sp,params]=kf_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);
%
% Optimal Spectral Tracking (OST) analysis (Type 2)
%
% Input parameters
%       dat1            Data vector channel 1
%       dat2            Data vector channel 2
%       trig            Trigger times
%       offset          Trial offset (msecs) (scalar or vector)
%       duration        Trial duration (msecs) (scalar or vector)
%       rate            Sampling rate
%       seg_pwr         Segment power for Fourier transform (no padding for seg_pwr=0)
%       opt_str         Options string
%                           r<0|1|2>    Rectify channels
%                           W<W>        Multitaper time-bandwidth (default 2.5, 0-trial avg.)
%                           Q<q>        Process noise (default 0.1)
%                           f<fmax>     Maximum frequency (default Nyquist)
%                           s           Use Kalman smoother
%                           n           Power normalise segments
%                           k<n>        KF method (see implementations below; default 2)
%                           j           Jackknife variance
%                           m           Mains artifact removal
%                           a           Alpha-rate adaptation of Jazwinski algorithm (log-transform only for now)
%                           l           Legacy mode (Neurospec compatible)
%
% Output parameters
%       sp              Matrix of spectral components, cols: (3rd dim are trials)
%                           1           Auto-spectra ch.1
%                           2           Auto-spectra ch.2
%                           3           Cross-spectra ch.1,2
%                           4           Coherence
%       params          Parameters structure
%
% KF implementations
%       0 - null        Return periodogram estimates (parameters structure will not fully correlate)
%       1 - kf          Time-domain filtering (Inaccurate distributional assumptions)
%       2 - kflog       Log-transform filtering
%       3 - ff          Exponential decay (forgetting-factor) (Q<q> becomes Q<delta(<1)>)
%       4 - kfwnd       Sliding window (Q<q> becomes Q<K>, number of trials)
%       5 - kfroot      Root-transform (optimal choice depending on analysis parameters)
%
%function [sp,params]=kf_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);

% Check input parameters
if (isempty(seg_pwr))
    seg_pwr = 0;
end;
if ((2^seg_pwr<duration*rate/1000) & (seg_pwr~=0))
    error(' seg_pwr too small for data');
end;

% Default parameters
multitaper=false;           % Perform multitaper analysis
smooth=false;               % Backward kalman smoothing
jackknife=false;            % Jackknife variance
varlen=false;               % Variable length segments
normalise=false;            % Energy normalisation per segment
mains=false;                % Mains supression
method=2;                   % KF method
Q=0.1;                      % Process noise
plotf=rate/2;               % Maximum analysis frequency
mainsline=50;               % Mains line frequency
jazalph=[];                 % Jazwinski alpha-rate
legacy = false;             % Legacy (NeuroSpec) mode

% Parse opt_str for Kalman-Fourier related options
options=deblank(opt_str);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'r'                    % Rectify channels
            n=str2num(optarg);
			if ((n<0) | (n>2))
			    error(['Error in option argument -- r' optarg]);
			end;
			if (n~=1)   % Rectify ch 1.
			    dat1=abs(dat1-mean(dat1));
			end;  
			if (n>=1)   % Rectify ch 2.
			    dat2=abs(dat2-mean(dat2));
			end;
        case 'W'                    % Apply multitaper method
            if (isempty(optarg))
                W=2.5;
            else
                W=str2num(optarg);
            end;
            if (W>0)
                multitaper=true;
            end;
        case 'Q'                    % Q%age (time-domain kf only)
            Q=str2num(optarg);
        case 'f'                    % maximum plot frequency
            plotf=str2num(optarg);
        case 'm'                    % Supress mains
            mains=true;
        case 's'                    % smoother
            smooth=true;
        case 'n'                    % power normalise segments
            normalise=true;
        case 'k'                    % KF method
            method=str2num(optarg);
            if (isempty(method))
                error([' Error in parameter ' opt]);
            end;
        case 'j'                    % jackknife
            jackknife=true;
        case 'a'                    % Jazwinksi alpha-rate
            jazalph=str2num(optarg);
        case 'l'                    % Legacy (NeuroSpec) mode
            legacy = true;
    	otherwise                   % Pass on to analysis
            opt_str=[opt_str ' ' opt];
    end;
end;

% Apply options to data
[dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);

% Determine truncated frequency range
if (seg_pwr==0), seg_len=max(duration*rate/1000);       % min possible (memory saving - but slower)
else             seg_len=2^seg_pwr; end;                % power of 2 (incorporates zero padding)
freqs=rate*[1:floor(seg_len/2)]/seg_len;
f_max=max(find(freqs<=plotf));          % Determine freq. plot range
freqs=freqs(1:f_max);

% Ensure triggered segments fall within data length
inrange=((trig+offset*rate/1000)>0 & (trig+(offset+duration)*rate/1000)<=size(dat1,1));
trig=trig(inrange);
if (length(offset)>1), offset=offset(inrange); end;
if (length(duration)>1), duration=duration(inrange); end;

% Determine data parameters
N=max(duration*rate/1000);
M=length(trig);
switch (length(offset))
    case 1, offset=offset*ones(M,1);
    case M, % Do nothing
    otherwise error(' Offset must be a scalar or vector of the same size as trig.');
end;
switch (length(duration))
    case 1, duration=duration*ones(M,1);
    case M, varlen=1;
    otherwise error(' Duration must be a scalar or vector of the same size as trig.');
end;

% Determine tapers
if ((~varlen) & (multitaper))
    % Get tapers for MTM (Using 3rd party software by E.Breitenberger)
    if (W==0.5)
        [E,V]=dpss(N,1);
        E=E(:,1); V=V(1);
    else
        [E,V]=dpss(N,W);
    end;
else
	V=1;
end;

% Allocate memory
sp=zeros(f_max,3,M);

% Perform multitaper analysis on each trial
disp('Generating single segment estimators');
for seg=1:M
    % Determine time range
    T=round(duration(seg)*rate/1000);
    trange=round(trig(seg))+offset(seg)*rate/1000+[0:T-1];
    if (varlen & multitaper)                % If variable length segment duration
        %disp(['MTM analysis: ' int2str(seg) ' of ' int2str(M)]);
        [E,V]=dpss(T,W);
    end;
    % Set Parcevals normalisation for rectangular windowing (consistent with MTM levels)
    if (~multitaper)
        E=(1/sqrt(T))*ones(T,1);
    end;
    % Set spectral norm (Only works if multitaper OFF)
    spec_norm = 1;
    if (legacy)
        spec_norm=2*pi;%*length(trange);
    end;
    % Extract segment from data
    dat11=dat1(trange); dat22=dat2(trange);
    % Zero mean segments
    dat11=dat11-mean(dat11); dat22=dat22-mean(dat22);
    % Power normalise segments
    if (normalise)
        dat11=dat11/sqrt(sum(dat11.^2)/T);
        dat22=dat22/sqrt(sum(dat22.^2)/T);
    end;
    % Generate spectra over tapers
    for ind=1:size(E,2)
        F1=fft(E(:,ind).*dat11,seg_len); F1=F1(2:f_max+1);
        F2=fft(E(:,ind).*dat22,seg_len); F2=F2(2:f_max+1);
        sp(:,1,seg)=sp(:,1,seg)+(abs(F1).^2)/size(E,2)/spec_norm;
        sp(:,2,seg)=sp(:,2,seg)+(abs(F2).^2)/size(E,2)/spec_norm;
        sp(:,3,seg)=sp(:,3,seg)+(F1.*conj(F2))/size(E,2)/spec_norm;
        
        %sp(:,5,seg)=F1/sqrt(size(E,2)); sp(:,6,seg)=F2/sqrt(size(E,2));
    end;
    % Construct confidence limit structure
    cl(seg).seg_tot=1/sum((V/sum(V)).^2);   % no. segments per estimate
    cl(seg).seg_len=T;                      % samples per segment
    cl(seg).pad_len=seg_len;
end;

% Hybrid Kalman-Fourier analysis
switch (method)
    case 0  % null (return sp unprocessed)
        sp(:,4,:)=abs(sp(:,3,:)).^2./(sp(:,1,:).*sp(:,2,:));
        sp(:,5,:)=angle(sp(:,3,:));
        params.Psp=psi(1,cl(1).seg_tot)*ones(f_max,2,M);
        params.L=ones(f_max,3,M);
    case 1  % time-domain
        [sp,params]=kf_tfspm2w(sp,cl,Q/100,smooth,jackknife);
    case 2  % log-domain
        [sp,params]=kflog_tfspm2w_3(single(sp),cl,Q,smooth,jackknife,jazalph);
    case 3  % exponential decay
        [sp,params]=ff_tfspm2w(sp,cl,Q,smooth,jackknife);
    case 4  % sliding window
        [sp,params]=kfwnd_tfspm2w(sp,cl,Q,smooth,jackknife);
        M=size(sp,3);
    case 5  % root-transform
        [sp,params]=kfroot_tfspm2w(single(sp),cl,Q,smooth,jackknife,jazalph);
    case 6  % log-domain with covariance process noise structure
        [sp,params]=kflog_tfspm2w_covar(single(sp),cl,Q,smooth,jackknife,jazalph);
    case 7  % log-domain spectra, atanh |coherency|
        sp(:,3,:) = (sp(:,3,:)./sqrt(sp(:,1,:).*sp(:,2,:)));
        % Bias subtracted
        %sp(:,3,:) = sqrt( abs(sp(:,3,:)./sqrt(sp(:,1,:).*sp(:,2,:))).^2 - 1/size(E,2));
        
        [sp,params]=kflog_tfspm2w_4(single(sp),cl,Q,smooth,jackknife,jazalph);
        sp(:,3,:) = abs( sp(:,3,:) ).^2;
    otherwise
        error(' Unknown KF method specified.');
end;

% Mains supression
if (mains)
    lineind=dsearchn(freqs',mainsline);
    % Spectra
    sp(lineind,1:3,:)=0.5*(sp(lineind-2,1:3,:)+sp(lineind+2,1:3,:));
    sp(lineind-1,1:3,:)=0.5*(sp(lineind-2,1:3,:)+sp(lineind-3,1:3,:));
    sp(lineind+1,1:3,:)=0.5*(sp(lineind+2,1:3,:)+sp(lineind+3,1:3,:));
    % Smooth elements in upper hermetian section of cross spectral estimate.
    % This data used in ifft() to generate cumulant. Data is complex conjugate.
    %sp(segsize-lineind+2,3,:)=conj(sp(lineind,3,:));
    %sp(segsize-lineind+3,3,:)=conj(sp(lineind-1,3,:));
    %sp(segsize-lineind+1,3,:)=conj(sp(lineind+1,3,:));
    % Recalculate phase and coherence
    sp(:,4,:)=(abs(sp(:,3,:)).^2)./(sp(:,1,:).*sp(:,2,:));
    sp(:,5,:)=angle(sp(:,3,:));
end;

% Consolidate effective no. segments
params.L=reshape(params.L,f_max,3,M);           % Effective No. segments per (freq,channel,trial)
tapercount=1;
if (multitaper)
    tapercount=size(E,2);
    params.L=single(double(params.L)*tapercount);
end;
params.Lprime=squeeze(min(params.L,[],2));      % Estimate L' as min of spectral L (most conservative)

% Form parameters structure (starting with input arguments)
params.trig=trig;
params.offset=offset; %ms
params.duration=duration; %ms
params.rate=rate;
params.seg_pwr=seg_pwr;
params.opt_str=opt_str;
params.Q=Q;
params.tapercount=tapercount;
params.smooth=smooth;
params.jackknife=jackknife;
params.normalise=normalise;
params.inrange=inrange;
% Additional parameters
params.freqs=freqs;
params.M=M;
params.method=method;
% Optional parameters
if (exist('W'))
    params.NW=W;
end;
