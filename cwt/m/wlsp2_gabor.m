function [Wsp,wlparam]=wlsp2_gabor(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);
%function [Wsp,wlparam]=wlsp2_gabor(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);
%
% Spectrum Gabor atom analysis Type 2
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
%       dt          Time resolution
%       df          Frequency resolution
%       frange      Frequency range (may be [] for defaults)
%       opt_str     Options string
%                       r<0,1,2>    rectify (none, ch1, ch1&2)
%                       p           pad with zeros
%   Wavelet parameters
%       mother      Mother wavelet
%       wlopt       Wavelet options
%
% Ouput parameters
%   Wsp - Spectral coefficients
%       Col 1:  Auto-spectra ch.1
%       Col 2:  Auto-spectra ch.2
%       Col 3:  Cross-spectra
%       Col 4:  Coherence
%       Col 5:  Phase
%
%function [Wsp,wlparam]=wlsp2_gabor(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);

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
            opt_str=[opt_str opt ' '];
    end;
end;

% Determine data parameters
dt=1/rate;
time=((offset*rate/1000):((offset+duration)*rate/1000));

% Segregate data into active regions
N=length(dat1);
offset=offset*rate/1000;        % Convert ms -> samples
duration=duration*rate/1000;    %
trigstart=trig+offset;
trigstart=trigstart((trigstart>0) & (trigstart<(N-duration-offset+1)));     % Remove incomplete segments
trigcount=length(trigstart);
for ind=1:trigcount
    segdat1(:,ind)=dat1(trigstart(ind):(trigstart(ind)+duration));
    segdat2(:,ind)=dat2(trigstart(ind):(trigstart(ind)+duration));
end;
segcount=trigcount;

% Wavelet transform the individual segments
[Wsp1,wlparam]=wlfGaborTransform(segdat1(:,1),dt,df,frange,opt_str,mother,wlopt);   % Determine freqs
Wsp=zeros(length(wlparam.freqs),length(time),5);
for seg=1:segcount
    disp([' Gabor transform segment ' int2str(seg) ' of ' int2str(segcount)]);
    [Wsp1]=wlfGaborTransform(segdat1(:,seg),dt,df,frange,opt_str,mother,wlopt);
    [Wsp2]=wlfGaborTransform(segdat2(:,seg),dt,df,frange,opt_str,mother,wlopt);
    
    % Recursively form auto/cross spectra
    Wsp(:,:,1)=(1/seg)*((seg-1)*Wsp(:,:,1)+(abs(Wsp1).^2));
    Wsp(:,:,2)=(1/seg)*((seg-1)*Wsp(:,:,2)+(abs(Wsp2).^2));
    Wsp(:,:,3)=(1/seg)*((seg-1)*Wsp(:,:,3)+(Wsp1.*conj(Wsp2)));
end;

% Wavelet coherence
Wsp(:,:,4)=(abs(Wsp(:,:,3)).^2)./(Wsp(:,:,1).*Wsp(:,:,2));

% Wavelet phase
Wsp(:,:,5)=angle(Wsp(:,:,3));
