function ah = mwpsp(time,dat1,dat2,Sx,Sy,Sxy,wlparam,opts)
%function mwpsp(time,dat1,dat2,Sx,Sy,Sxy,wlparam,[opts])
%
% Multi-wavelet data, spectra and coherence plotting routine.
% Includes allocation of 95% confidence limits on coherence.
% Also compatible with multitaper hermite functions.
%
% Input parameters
%       time        Time vector
%       dat1        Data for channel 1
%       dat2        Data for channel 2
%       Sx          Multi-wavelet spectra ch.1
%       Sy          Multi-wavelet spectra ch.2
%       Sxy         Multi-wavelet cross-spectra ch.1,2
%       wlparam     Multi-wavelet transform parameters
%
%function mwpsp(time,dat1,dat2,Sx,Sy,Sxy,wlparam,[opts])

% Check input parameters
if (~exist('opts'))
    opts = '';
end;

% Determine statistical properties
L=wlparam.L;
speccl=[];%0.851/sqrt(L);           % PBMB (6.2)
specopt=['c ' opts];
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

% Create spectral time-vector
sptime=(1:size(Sx,2))/size(Sx,2)*(time(end)-time(1))+time(1);

% Plot wavelet spectra
figure;
%colormap(flipud(gray));
ah(1)=subplot(5,1,1);       % Data ch.1
    plot(time,dat1,'k');
    xlim(time([1 end]));
    ylim([min(dat1) max(dat1)]);
    title(['Channel 1 ' paramstr]);
ah(2)=subplot(5,1,2);       % Autospectra ch.1
    wlpsp_a(sptime,Sx,wlparam,speccl,specopt);
    %grid on
ah(3)=subplot(5,1,3);       % Data ch.2
    plot(time,dat2,'k');
    xlim(time([1 end]));
    ylim([min(dat2) max(dat2)]);
    title('Channel 2');
ah(4)=subplot(5,1,4);             % Autospectra ch.2
    wlpsp_a(sptime,Sy,wlparam,speccl,specopt);
    %grid on
ah(5)=subplot(5,1,5);             % Coherence
    wlpsp_coh(sptime,coh,wlparam,R95,'c');%[wlparam.freqs' R95*ones(size(coh,1),1)],'b');
    %grid on
    %fig2a4l;

% Reposition raw data plots (alignment accounts for colourbars)
% oldpos=get(ah,'pos');
% newpos=get(ah(2),'pos');
% for ind = [ 1 3 ]
%     set(ah(ind),'pos',[oldpos{ind}(1) oldpos{ind}(2) newpos(3) oldpos{ind}(4)]);
% end;

% newpos = get(ah(5),'pos');
% for k = [1 3]
%     pos = get(ah(k),'pos');
%     set(ah(k),'pos',[pos(1) pos(2) newpos(3) pos(4)]);
% end;

if ( exist('linkaxes','file') )
    linkaxes(ah,'x');
end;
