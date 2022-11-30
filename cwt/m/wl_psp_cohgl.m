function glPower=wl_psp_cohgl(time,Wsp,wlparam,trange);
%function glPower=wl_psp_cohgl(time,Wsp,wlparam,[trange]);
%
% Plot scalogram of coherence & global wavelet power
%
% Input parameters
%   time        Time vector (ms)
%   Wsp         Wavelet spectral coefficients
%   wlparam     Wavelet and transform parameters
%   trange      (Opt) Time range [tmin tmax] (in ms)
%
% Output parameters
%   ax          2-element vector containing axis handles
%
%function wl_psp_cohgl(time,Wsp,wlparam,[trange]);

% Check input parameters
t_limits=0;
if (nargin>3)
    tmin=trange(1);
    tmax=trange(2);
    t_limits=1;
end;

% Reduce time within trange
if (t_limits)
    % Determine start and stop samples
    tstart=dsearchn(time',tmin);
    tstop=dsearchn(time',tmax);
    tt=(tstart:tstop);
    
    % Reduce time of input parameters
    time=time(tt);
    Wsp=Wsp(:,tt,:);
    wlparam.coi=wlparam.coi(tt);
    wlparam.gb_coi=wlparam.gb_coi-tstart;
end;

% Determine scalogram plot and global power
coh=Wsp(:,:,4);
glPower=wl_cohgl(Wsp,wlparam);

% Setup axes
figure;
ax(1)=subplot(2,1,1);
ax(2)=subplot(2,1,2);

% Scalogram
subplot(ax(1));
set(ax(1),'position',[.1 .1 .7 .8]);
scalogram(time,coh,wlparam);
title(['coh: ' wlparam.label]);

% Global Wavelet power with and without COI
freqs=wlparam.freqs;
freqlim=[min(wlparam.freqs) max(wlparam.freqs)];
if (~wlparam.linearscale)
    freqs=log2(freqs);
    freqlim=log2(freqlim);
end;
subplot(ax(2));
set(ax(2),'position',[.85 .1 .1 .8]);
plot(glPower,freqs,'k');
set(ax(2), 'ylim', freqlim,'yticklabel',[],'ygrid','on');
title('Global coherence');
