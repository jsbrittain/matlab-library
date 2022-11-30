function [W,wlparam]=wlsp2(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);
%function [W,wlparam]=wlsp2(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);
%
% Spectrum wavelet analysis Type 2
%
% Calculate auto-spectra, cross-spectra, coherenc and phase for
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
%       df          Frequency resolution
%       frange      Frequency range (may be [] for defaults)
%       opt_str     Options string (r0,r1,r2=rectify; p=pad with zeros)
%   Wavelet parameters
%       mother      Mother wavelet
%       wlopt       Wavelet options
%
% Ouput parameters
%   W - Spectral coefficients
%       Col 1:  Auto-spectra ch.1
%       Col 2:  Auto-spectra ch.2
%       Col 3:  Cross-spectra
%       Col 4:  Coherence
%       Col 5:  Phase
%
%function [W,freqs,coi]=wlsp2(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);

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
time=((offset*rate/1000):((offset+duration)*rate/1000));

% Segregate data into active regions
offset=offset*rate/1000;
duration=duration*rate/1000;
trigstart=trig+offset;
trigcount=length(trigstart);
for ind=1:trigcount
    segdat1(:,ind)=dat1(trigstart(ind):(trigstart(ind)+duration));
    segdat2(:,ind)=dat2(trigstart(ind):(trigstart(ind)+duration));
end;
segcount=trigcount;

% Wavelet transform the individual segments
[W1,wlparam]=wlfTransform(segdat1(:,1),dt,df,frange,opt_str,mother,wlopt);   % Determine freqs
W=zeros(length(wlparam.freqs),length(time),5);
for seg=1:segcount
    disp([' Wavelet transform segment ' int2str(seg) ' of ' int2str(segcount)]);
    [W1]=wlfTransform(segdat1(:,seg),dt,df,frange,opt_str,mother,wlopt);
    [W2]=wlfTransform(segdat2(:,seg),dt,df,frange,opt_str,mother,wlopt);
    
    % Recursively form auto/cross spectra
    W(:,:,1)=(1/(seg+1))*(seg*W(:,:,1)+(abs(W1).^2));
    W(:,:,2)=(1/(seg+1))*(seg*W(:,:,2)+(abs(W2).^2));
    W(:,:,3)=(1/(seg+1))*(seg*W(:,:,3)+(W1.*conj(W2)));
end;

% Wavelet coherence
W(:,:,4)=(abs(W(:,:,3)).^2)./(W(:,:,1).*W(:,:,2));

% Wavelet phase
W(:,:,5)=angle(W(:,:,3));
