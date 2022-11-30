function stftSummary(time,offset,dat,f,frange,ch,logplot);
%function stftSummary(time,offset,dat,f,[ch,[logplot]]);
%
% Compatibility function to plot STFT in a scalogram style with
% time series and global spectra
%
% Input parameters
%   time        Time vector
%   offset      Offset
%   dat         Time series
%   f           Frequency-domain parameters from NeuroSpec
%   frange      Frequency range for plotting (may be empty)
%   ch          (opt) Channel for plotting (1=time,2=a11,3=a22,4=ch,5=ph; default 1)
%   logplot     (opt) Display power on a log10 scale
%                     (0=none, 1=spectrogram, 2=spectrogram plus global power)
%

% Determine input parameters
if (nargin<7)
    logplot=1;
end;
if (nargin<6)
    ch=2;
end;

% Reduce f matrix within frange limits
if (~isempty(frange))
    fmin=frange(1);
    fmax=frange(2);
    min_pt=dsearchn(f(:,1,1),fmin);
    max_pt=dsearchn(f(:,1,1),fmax);
    f=f(min_pt:max_pt,:,:);
end;

% Convert spectra to linear scale
power=10.^flipud(reshape(f(:,ch,:),size(f,1),size(f,3)));
if (logplot==2)
    power=log10(power);
end;
glPower=flipud(mean(power,2));
if (logplot==1)
    power=log10(power);
end;

% Generate a wlparam for STFT so we can use the scalosummary function
wlparam.display_coi=0;
wlparam.linearscale=1;
wlparam.freqs=f(:,1,1);
freqs=f(:,1,1);
freqlim=[min(wlparam.freqs) max(wlparam.freqs)];

% Plot time series
figure;
ax(1)=subplot(3,1,1); ax(2)=subplot(3,1,2); ax(3)=subplot(3,1,3);
subplot(ax(1));
set(ax(1),'position',[.1 .85 .7 .1]);
plot(time,dat,'k-');
xlim([min(time) max(time)]);

% Plot (scalo-)spectrogram
subplot(ax(2));
set(ax(2),'position',[.1 .1 .7 .7]);
scalogram(offset,power,wlparam);

% Plot global power
subplot(ax(3));
set(ax(3),'position',[.85 .1 .1 .7]);
plot(glPower,freqs,'k-');
set(ax(3), 'ylim', freqlim,'yticklabel',[],'ygrid','on');
fig2a4l;
