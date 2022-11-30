function [time,q,c95]=kf_nae(dat,rate,trig,duration,onset,offset,opt_str)
%function [time,q,c95]=kf_nae(dat,rate,trig,duration,onset,offset,opt_str)
%
% Generate Normalised Activity Envelope (NAE) for variable length duration
% segments
%
% This routine generates the average Fourier transform from a triggered
% region (including zero padding in the frequency domain) and inverse
% transforms to display an average (normalised length) activity envelope.
% Segment energy is also normalised to prevent longer envelopes swamping
% the analysis.
%
% Power normalisation applied so that each trial has sample variance 1,
% thus the mean of L trials has sample variance 1/L and 95% conf. under
% null hypothesis is "0 +/- 1.96*sqrt(1/L)".  Produces same shape as energy
% normalisation.
%
% Input parameters
%       dat             Time series
%       rate            Sampling rate
%       trig            Trigger times (samples)
%       duration        Duration vector (msecs)
%       onset           Pre-duration offset (percentile of duration)
%       offset          Post-duration offset (percentile of duration)
%       opt_str         Options string
%                           r - Rectify time series
%                           u - Do not power normalise
%                           c -Summary figure plot of construction process
%
% Output parameters
%       time            Time axis (normalised to triggered region)
%       q               Average activity envelope of energy normalised segments
%       c95             Upper and lower 95% confidence interval about 0 (independence)
%
%function [time,q,c95]=kf_nae(dat,rate,trig,duration,onset,offset,opt_str)

% Check input parameters
if (nargin~=7)
    error(' Incorrect number of input parameters');
end

% Default parameters
powernorm=true;         % Power normalise
rectify=false;
diagram=false;

% Parse options string
options=deblank(opt_str); opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'u'                    % Do not power normalise
            powernorm=false;
	    case 'r'                    % Rectify channels
        	rectify=true;
        case 'c'                    % Construction diagram
        	diagram=true;
        otherwise                   % Options for wavelet analysis
            error(['Unknown option -- ' opt]);
    end
end

% Conserve original data for diagram
if (diagram)
    udat=dat;
end
    
% Rectify data
if (rectify)
    dat=abs(dat-nanmean(dat));
end

% Check that duration is a vector
if (length(duration)==1)
    duration=duration*ones(length(trig),1);
end

% Crop triggers to only include those segments fully defined in the time range
dur=duration*rate/1000;
trigpos=(((trig+round(onset*dur))>0) & ((trig+round(offset*dur))<=length(dat)));
dur=dur(trigpos); trig=trig(trigpos);

% Determine parameters
maxdur=max(duration)*rate/1000;     % samples
trigcount=length(trig);
N=round((offset-onset)*maxdur+1);

% Reserve variable space
F11=zeros(N,1);

% Construction diagram
if (diagram)
    figure;
    freqs=(rate/2)*linspace(-1,1,N);
end

% Time vector (also required in construction diagram)
time=((onset*maxdur)+[0:N-1])/maxdur;

% Form spectra
%q=zeros(N,1);
actualtrigcount = 0;
for ind=1:trigcount
    % Determine data segment
    trange=trig(ind)+round(onset*dur(ind))+[0:round((offset-onset)*dur(ind))];
    if ( any(isnan(dat(trange))) )
        % Reject epochs containing NaN samples
        continue;
    end
    actualtrigcount = actualtrigcount + 1;
    % Zero pad Fourier transform (zero mean and energy normalise segment first)
    dat1=dat(trange)-mean(dat(trange));     % Zero mean
    if (powernorm)
        dat1=dat1/sqrt(sum(dat1.^2)/length(dat1)); % Power normalise (req. for conf. limits)
    end
    F1=fft(dat1)/sqrt(length(dat1));        % Fourier transform
    zcount=N-length(F1);                    % Zero pad count
    F1=[F1(1:ceil(length(F1)/2)); ...       % Zero pad in freq domain
        zeros(zcount,1); ...                %  |
        F1(ceil(length(F1)/2+1):end)];      %  |
    % Apply power correction factor for zero padding
    F1=F1*sqrt(N/length(dat1));
    % Form average Fourier transform
    F11=F11+F1;%/trigcount;     % Divisor later to account for NaN segments
    % Apply power correction factor and average in time-domain
    %q=q+(sqrt(N/length(dat1)))*real(ifft(F1)*sqrt(length(F1)))/trigcount;
    
    % Construction diagram
    if (diagram)
        if ((ind==1) || (ind==5) || (ind==15))
            tt=[0:length(trange)-1]*1000/rate;
            if (ind==1)
                % EMG burst
                subplot(3,5,1); plot(tt,udat(trange)-mean(udat(trange)),'k');
                % Rectified EMG burst
                subplot(3,5,2); plot(tt,dat(trange),'k');
                % FFT (magnitude FFT with highlighted zero padding)
                subplot(3,5,3); plot(freqs,abs(F1),'k');
                % Inverse FFT
                subplot(3,5,4); plot(time,real(ifft(F1))*sqrt(length(F1)),'k');
            elseif (ind==5)
                subplot(3,5,6); plot(tt,udat(trange)-mean(udat(trange)),'k');
                subplot(3,5,7); plot(tt,dat(trange),'k');
                subplot(3,5,8); plot(freqs,abs(F1),'k');
                subplot(3,5,9); plot(time,real(ifft(F1))*sqrt(length(F1)),'k');
            elseif (ind==15)
                subplot(3,5,11); plot(tt,udat(trange)-mean(udat(trange)),'k');
                subplot(3,5,12); plot(tt,dat(trange),'k');
                subplot(3,5,13); plot(freqs,abs(F1),'k');
                subplot(3,5,14); plot(time,real(ifft(F1))*sqrt(length(F1)),'k');
            end
        end
    end
end
F11 = F11/actualtrigcount;

% Determine envelope
q=real(ifft(F11)*sqrt(length(F11)));

% Construct 95% confidence interval (Valid only for power normalisation)
c95=1.96/sqrt(trigcount);

% Construction process
if (diagram)
    subplot(3,5,10); kf_psp_nae(time,q,c95);
end
