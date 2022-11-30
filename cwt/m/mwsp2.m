function [Wsp,wlparam]=mwsp2(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);
%function [Wsp,wlparam]=mwsp2(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);
%
% Multi-wavelet trial-averaging wavelet spectral analysis: Type 2
%
% Calculate auto-spectra, cross-spectra, coherence and phase for
% two channels.
%
% Recursive implementation
%
% Input parameters
%   Data parameters
%       dat1        Time series 1
%       dat2        Time series 2
%       trig        Trigger
%       offset      Offset (ms)
%       duration    Duration (ms)
%       rate        Sampling rate
%   Analysis parameters
%       df          Frequency resolution (1/(scales/octave))
%       frange      Frequency range (may be [] for defaults)
%       opt_str     Options string
%                       r<0,1,2>    rectify (none, ch1, ch1&2)
%                       p           pad with zeros
%   Wavelet parameters
%       mother      Mother wavelet
%       wlopt       Wavelet options
%
% Ouput parameters
%       Wsp - Spectral coefficients (3rd dimension)
%                   1 - Auto-spectra ch.1
%                   2 - Auto-spectra ch.2
%                   3 - Cross-spectra
%                   4 - Coherence
%                   5 - Phase
%       wlparam     Wavelet parameters
%
%function [Wsp,wlparam]=mwsp2(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);

% Last edited 30-01-2006 (JSB)

% Parse options string
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
			    dat1=abs(dat1);
			end;  
			if (n>=1)   % Rectify ch 2.
			    dat2=abs(dat2);
			end;
        otherwise                   % Options for wavelet analysis
            opt_str=[opt_str opt optarg ' '];
    end;
end;

% Determine data parameters
dt=1/rate;
time=((offset*rate/1000):((offset+duration)*rate/1000-1));

% Segregate data into active regions
N=length(dat1);
offset=offset*rate/1000;        % Convert ms -> samples
duration=duration*rate/1000;    %
trigstart=trig+offset;
trigstart=trigstart((trigstart>0) & (trigstart<=(N-duration-offset+1)));     % Remove incomplete segments
trigcount=length(trigstart);
for ind=1:trigcount
    segdat1(:,ind)=dat1(trigstart(ind):(trigstart(ind)+duration-1));
    segdat2(:,ind)=dat2(trigstart(ind):(trigstart(ind)+duration-1));
end;
segcount=trigcount;

% Wavelet transform the individual segments
[Wsp1,wlparam]=mwfTransform(segdat1(:,1),dt,df,frange,opt_str,mother,wlopt);   % Determine freqs
wlparam=numericalCOI(wlparam,duration);
Wsp=zeros(length(wlparam.freqs),length(time),5);
for seg=1:segcount
    disp([' Wavelet transform segment ' int2str(seg) ' of ' int2str(segcount)]);
    [Wsp1,Wsp2,Wsp12,param]=mwfTransform(segdat1(:,seg),segdat2(:,seg),dt,df,frange,opt_str,mother,wlopt);
    
    % Recursively form auto/cross spectra
    Wsp(:,:,1)=(1/seg)*((seg-1)*Wsp(:,:,1)+Wsp1);
    Wsp(:,:,2)=(1/seg)*((seg-1)*Wsp(:,:,2)+Wsp2);
    Wsp(:,:,3)=(1/seg)*((seg-1)*Wsp(:,:,3)+Wsp12);
end;

% Wavelet coherence
Wsp(:,:,4)=(abs(Wsp(:,:,3)).^2)./(Wsp(:,:,1).*Wsp(:,:,2));

% Wavelet phase
Wsp(:,:,5)=angle(Wsp(:,:,3));

% Truncate data to COI cutoff
[wlparam,Wsp]=mwTruncate(wlparam,Wsp);
