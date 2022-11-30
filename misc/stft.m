function S=stft(dat,rate,opt,freqs);
%function S=stft(dat,rate,[opt,[freqs]]);
%
% STFT
%
% Input parameters
%   dat         Data vector
%   rate        Sampling rate (Hz)
%   opt         Options
%                   i           Treat as an infinite series
%                   w<n,opts>   Window function
%                                   Boxcar:     0,<duration-ms>,<step-ms>
%                                   Gaussian:   1,<stddev-ms>,<step-ms>
%   freqs       Set of frequencies for evaluation (Inf series)
%
% Output parameters
%   S           Time-dependent spectra
%
%function S=stft(dat,rate,[opt,[freqs]]);

% Check input parameters
if (nargin<3)
    error(' Incorrect number of input parameters');
end;
if (~exist('opt'))
    opt='';
end;
N=length(dat);

% Default parameters
data_window=0;      % Boxcar
infinite_series=0;  % Treat as an infinite series

% Parse options string
options=deblank(opt);
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'i'            % Treat as an infinite series
        	infinite_series=1;
        case 'w'            % Choose window envelope
            [tmpstr,remainder]=strtok(optarg,',');
            data_window=str2num(tmpstr);
            switch (data_window)
                case 0      % -- Boxcar
                    [tmpstr,remainder]=strtok(remainder,',');
                    wnd_duration=str2num(tmpstr)*rate/1000;     % Duration
                    [tmpstr,remainder]=strtok(remainder,',');
                    wnd_step=str2num(tmpstr)*rate/1000;         % Step
                    wnd_count=(N-wnd_duration)/wnd_step+1;
                case 1      % -- Gaussian
                    [tmpstr,remainder]=strtok(remainder,',');
                    wnd_sd=str2num(tmpstr)*rate/1000;           % Std.dev.
                    [tmpstr,remainder]=strtok(remainder,',');
                    wnd_step=str2num(tmpstr)*rate/1000;         % Step
                    wnd_count=N/wnd_step;
                    infinite_series=1;
            end;
    	otherwise
            error (['Illegal option -- ',opt]);
    end;
end;

% Convert msecs to samples
time=[0:wnd_count-1]*wnd_step;
dt=1/rate;

% Evaluate STFT per window
for ind=1:wnd_count    
    % Data window
    switch (data_window)
        case 0      % Boxcar
            window=logical(zeros(1,N));
            window((ind-1)*wnd_step+1:(ind-1)*wnd_step+wnd_duration)=1;
        case 1      % Gaussian
            window=exp(-((([0:N-1]-wnd_step*(ind-1))/rate).^2)/(2*wnd_sd^2));
    end;
    
    % Evaluate Fourie transform
    if (~exist('freqs'))
        if (infinite_series)        
            % Evaluate FFT over data length with a window on top
            stmp=transpose(fft(dat.*window));
            Spect(:,ind)=abs(stmp(1:end/2+1)).^2;
            freqs=rate*[0:N/2]/N;
        else
            % Evaluate FFT over windowed region only
            stmp=transpose(fft(dat(window)));
            Spect(:,ind)=abs(stmp(1:wnd_duration/2+1)).^2;
            freqs=rate*[0:wnd_duration/2]/wnd_duration;
        end;
    else
        % Evaluate at a set of predefined frequencies
        for indf=1:length(freqs)
            Spect(indf,ind)=abs(dt*sum(dat.*window.*exp(-i*2*pi*freqs(indf)*[0:(length(dat)-1)]/rate))).^2;
		end;
    end;
end;

% Check output requirements
if (nargout~=1)
    contourf(time,freqs,Spect);
else
    S=Spect;
end;
