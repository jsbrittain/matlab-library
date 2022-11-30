function [sp,Psp,freqs]=kf_spm2w_mmd(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);
%function [sp,Psp,freqs]=kf_spm2w_mmd(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);
%
% Hybrid Kalman-Fourier analysis, type 2
%
% Dynamic Multi-Modal Adaptive Estimator
%   (NOT YET WORKING PROPERLY)
%
% Input parameters
%       dat1            Data vector channel 1
%       dat2            Data vector channel 2
%       trig            Trigger times
%       offset          Trial offset (msecs)
%       duration        Trial duration (msecs)
%       rate            Sampling rate
%       seg_pwr         Segment power for Fourier transform
%       opt_str         Options string
%                           m<W>        Multitaper analysis (using W; default 2.5)
%                           Q<%age>     Process noise adaptivity percentage (default 0.1%)
%                           f<fmax>     Maximum frequency (default Nyquist)
%                           s           Use Kalman smoother
%
% Output parameters
%       sp              Matrix of spectral components, cols: (3rd dim are trials)
%                           1       Auto-spectra ch.1
%                           2       Auto-spectra ch.2
%                           3       Cross-spectra ch.1,2
%                           4       Coherence
%       Psp             Estimated error covariance matrix for sp
%       freqs           Vector of frequency values
%
%function [sp,Psp,freqs]=kf_spm2w_mmd(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);

% Default parameters
multitaper=logical(0);
smooth=logical(0);
Qpercentage=0.1;
plotf=rate/2;

% Parse opt_str for Kalman-Fourier related options
options=deblank(opt_str);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'm'                    % Apply multitaper method
            multitaper=logical(1);
            if (isempty(optarg))
                W=2.5;
            else
                W=str2num(optarg);
            end;
        case 'Q'
            Qpercentage=str2num(optarg);
        case 'f'
            plotf=str2num(optarg);
        case 's'
            smooth=logical(1);
    	otherwise                   % Pass on to analysis
            opt_str=[opt_str ' ' opt];
    end;
end;

% Segmented Fourier analysis
if (multitaper)
    % Apply options to data
    [dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);
    N=duration*rate/1000;
	M=length(trig);
	[E,V]=dpss(N,W);
	fmax=floor(N/2);
	sp=zeros(fmax,3,M);
	for seg=1:M
        disp(['MTM analysis: ' int2str(seg) ' of ' int2str(M)]);
        trange=round(trig(seg))+offset*rate/1000+[0:duration*rate/1000-1];
        for ind=1:size(E,2)
            F1=fft(E(:,ind).*dat1(trange),2^seg_pwr); F1=F1(1:fmax);
            F2=fft(E(:,ind).*dat2(trange),2^seg_pwr); F2=F2(1:fmax);
            sp(:,1,seg)=(1/ind)*((ind-1)*sp(:,1,seg)+abs(F1).^2);
            sp(:,2,seg)=(1/ind)*((ind-1)*sp(:,2,seg)+abs(F2).^2);
            sp(:,3,seg)=(1/ind)*((ind-1)*sp(:,3,seg)+F1.*conj(F2));
        end;
        cl(seg).seg_tot=1;
	end;
else
    S=1;    % Trials per group
    [f,t,cl,sp,M]=kf_fourier2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);
end;

% Hybrid Kalman-Fourier analysis
[sp,Psp]=kf_tfspm2w_mmd(sp,cl,Qpercentage/100,smooth);

% Reduce plot range to that specified by fmax
seg_len=2^seg_pwr;
freqs=rate*[1:floor(seg_len/2)]/seg_len;
f_max=min(find(freqs>=plotf));          % Determine freq. plot range
freqs=freqs(1:f_max);
sp=sp(1:f_max,:,:);                     % Reduce spectral matrices
Psp=Psp(1:f_max,:,:);
