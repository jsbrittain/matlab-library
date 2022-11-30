function mwpsp_summary(time,dat1,dat2,Sx,Sy,Sxy,wlparam,sp11,sp22,sp12,params,title1,title2,xspec,displaycoherence);
%function mwpsp_summary(time,dat1,dat2,Sx,Sy,Sxy,wlparam,sp11,sp22,sp12,params,title1,title2,xspec,displaycoherence);
%
% Multi-wavelet data, spectra and coherence plotting routine.
% Includes allocation of 95% confidence limits on coherence.
% Also compatible with multitaper hermite functions.
%
% Input parameters
%       time              Time vector
%       dat1               Data for channel 1
%       dat2               Data for channel 2
%       Sx                 Multi-wavelet spectra ch.1
%       Sy                 Multi-wavelet spectra ch.2
%       Sxy                Multi-wavelet cross-spectra ch.1,2
%       wlparam            Multi-wavelet transform parameters
%       sp11               Multitaper spectrum ch.1 (from mt_ routines)
%       sp22               Multitaper spectrum ch.2
%       sp12               Multitaper cross-spectrum
%       params             Multitaper parameters
%       title1             Title for channel 1
%       title1             Title for channel 2
%       xspec              X-limit for multitaper spectra and coherence
%       displaycoherence   Display coherence (!!!not yet implemented!!!)
%
%function mwpsp_summary(time,dat1,dat2,Sx,Sy,Sxy,wlparam,sp11,sp22,sp12,params,title1,title2,xspec,displaycoherence);

% Determine statistical properties
L=wlparam.L;
speccl=[];%0.851/sqrt(L);           % PBMB (6.2)
specopt='b';
R95=(1-0.05^(1/(L-1)));             % PBMB (6.6)

% Generate numerical COI if required
if ((isempty(find(wlparam.coi(2:end-1)~=0))) & (~wlparam.linearscale))
    N=size(Sx,2);
    wlparam=numericalCOI(wlparam,N);
end;

% Truncate parameters
if (~wlparam.linearscale)
    [wlparam,Sx,Sy,Sxy]=mwTruncate(wlparam,Sx,Sy,Sxy);
end;
coh=abs(Sxy).^2./(Sx.*Sy);

% Determine parameters string
paramstr=['(' upper(wlparam.mother(1)) lower(wlparam.mother(2:end))];
for ind=1:length(wlparam.wlopt)
    if (ind==1)
        paramstr=[paramstr ' '];
    else
        paramstr=[paramstr ', '];
    end;
    paramstr=[paramstr wlparam.paramstr{ind} '=' num2str(wlparam.wlopt{ind})];
end;
paramstr=[paramstr ')'];

% Plot wavelet spectra
figure;
%colormap(flipud(gray));
ah(1)=subplot(5,1,1);       % Data ch.1
plot(time,dat1,'k');
xlim(time([1 end]));
ylim([min(dat1) max(dat1)]);
title(title1);
ah(2)=subplot(5,1,3);       % Data ch.2
plot(time,dat2,'k');
xlim(time([1 end]));
ylim([min(dat2) max(dat2)]);
title(title2);
newah=subplot(5,1,2);       % Autospectra ch.1
wlpsp_a(time,Sx,wlparam,speccl,specopt);
grid on
subplot(5,1,4);             % Autospectra ch.2
wlpsp_a(time,Sy,wlparam,speccl,specopt);
grid on
subplot(5,1,5);             % Coherence
wlpsp_coh(time,coh,wlparam,R95,'b');%[wlparam.freqs' R95*ones(size(coh,1),1)],'b');
grid on
fig2a4l;

% Reposition raw data plots (alignment accounts for colourbars)
oldpos=get(ah,'pos');
newpos=get(newah,'pos');
for ind=1:length(ah)
    set(ah(ind),'pos',[oldpos{ind}(1) oldpos{ind}(2) newpos(3) oldpos{ind}(4)]);
end;
