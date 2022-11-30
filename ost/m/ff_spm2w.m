function [sp,freqs,L]=ff_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);
%function [sp,freqs,L]=ff_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);
%
% Forgetting factor analysis
%
% Reduced vector implmentation
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
%                           d<delta>    Delta forgetting factor (default 0.95)
%                           f<fmax>     Maximum frequency (default Nyquist)
%                           s           Use Kalman smoother
%
% Output parameters
%       sp              Matrix of spectral components, cols: (3rd dim are trials)
%                           1           Auto-spectra ch.1
%                           2           Auto-spectra ch.2
%                           3           Cross-spectra ch.1,2
%                           4           Coherence
%       Psp             Estimated error covariance matrix for sp
%       freqs           Vector of frequency values
%
%function [sp,freqs,L]=ff_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);

% Default parameters
multitaper=logical(0);
smooth=logical(0);
delta=0.95;
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
        case 'd'
            delta=str2num(optarg);
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
    % Get tapers for MTM (Using 3rd party software by E.Breitenberger)
    N=duration*rate/1000;
	M=length(trig);
	[E,V]=dpss(N,W);
	fmax=floor(N/2);
	sp=zeros(fmax,3,M);
    % Perform multitaper method on each trial
	for seg=1:M
        disp(['MTM analysis: ' int2str(seg) ' of ' int2str(M)]);
        trange=round(trig(seg))+offset*rate/1000+[0:duration*rate/1000-1];
        for ind=1:size(E,2)           % Gaussian envelope if truncated to 1 taper             %%% <--- %%%
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

% Forgetting-factor analysis
fpts=size(sp,1);
vlen=3*fpts;
x=zeros(vlen,M);
xp=zeros(vlen,M+1);
initL=1;
xp(:,1)=reshape(mean(sp(:,:,1:initL),3),vlen,1);
for ind=1:M
    newz=reshape(sp(:,:,ind),vlen,1);
    x(:,ind)=(1-delta)*newz+delta*xp(:,ind);
    xp(:,ind+1)=x(:,ind);
end;
sp=reshape(x,fpts,3,M);
sp(:,4,:)=(abs(sp(:,3,:)).^2)./((sp(:,1,:)).*(sp(:,2,:)));  % Calculate coherence

% Reduce plot range to that specified by fmax
seg_len=2^seg_pwr;
freqs=rate*[1:floor(seg_len/2)]/seg_len;
f_max=min(find(freqs>=plotf));          % Determine freq. plot range
freqs=freqs(1:f_max);
sp=sp(1:f_max,:,:);                     % Reduce spectral matrices

% Determine effective no. segments
L=cumsum(delta.^[0:M-1]);
if (multitaper)
    L=L*size(E,2);
end;
