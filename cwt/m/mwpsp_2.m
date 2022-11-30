function mwpsp_2(time,dat1,dat2,orig_dat1,orig_dat2,Sx,Sy,Sxy,wlparam)
%function mwpsp_2(time,dat1,dat2,orig_dat1,orig_dat2,Sx,Sy,Sxy,wlparam)
%
% Multi-wavelet data, spectra and coherence plotting routine.
% Includes allocation of 95% confidence limits on coherence.
% Also compatible with multitaper hermite functions.
%
% Displays both processed and unprocessed data in channel plots.
%
% Input parameters
%       time        Time vector
%       dat1        Data for channel 1 (processed)
%       dat2        Data for channel 2 (processed)
%       orig_dat1   Original data for channel 1
%       orig_dat2   Original data for channel 2
%       Sx          Multi-wavelet spectra ch.1
%       Sy          Multi-wavelet spectra ch.2
%       Sxy         Multi-wavelet cross-spectra ch.1,2
%       wlparam     Multi-wavelet transform parameters
%
%function mwpsp_2(time,dat1,dat2,orig_dat1,orig_dat2,Sx,Sy,Sxy,wlparam)

% Determine statistical properties
L=wlparam.L;
speccl=[];%0.851/sqrt(L);           % PBMB (6.2)
specopt='';
R95=(1-0.05^(1/(L-1)));             % PBMB (6.6)

% Generate numerical COI if required
if ((isempty(find(wlparam.coi(2:end-1)~=0, 1))) && (~wlparam.linearscale) && (wlparam.display_coi))
    N=size(Sx,2);
    wlparam=numericalCOI(wlparam,N);
end;

% Truncate parameters
if ((~wlparam.linearscale) && (wlparam.display_coi))
    [wlparam,Sx,Sy,Sxy]=mwTruncate(wlparam,Sx,Sy,Sxy);
end;
coh=abs(Sxy).^2./(Sx.*Sy);

% Determine parameters string
paramstr=['(' upper(wlparam.mother(1)) lower(wlparam.mother(2:end)) ' '];
for ind=1:length(wlparam.wlopt)
    if (ind~=1), paramstr=[paramstr ', ']; end;
    paramstr=[paramstr wlparam.paramstr{ind} '=' num2str(wlparam.wlopt{ind})];
end; paramstr=[paramstr ')'];

% Plot wavelet spectra
figure; %fig2a4l;
%colormap(flipud(gray));
ah(1)=subplot(5,1,1);       % Data ch.1
    plot(time,dat1); hold('on');
    if (~isempty(orig_dat1)), plot(time,orig_dat1); end;
    xlim(time([1 end]));
    ylim([min([dat1; orig_dat1]) max([dat1; orig_dat1])]);
    set(gca,'ytick',[]);
    title(['Channel 1 ' paramstr]);
ah(2)=subplot(5,1,3);       % Data ch.2
    plot(time,dat2); hold('on');
    if (~isempty(orig_dat2)), plot(time,orig_dat2); end;
    xlim(time([1 end]));
    ylim([min([dat2; orig_dat2]) max([dat2; orig_dat2])]);
    set(gca,'ytick',[]);
    title('Channel 2');
subplot(5,1,2);       % Autospectra ch.1
    wlpsp_a(time,Sx,wlparam,speccl,specopt);
    %grid('on');
subplot(5,1,4);             % Autospectra ch.2
    wlpsp_a(time,Sy,wlparam,speccl,specopt);
    %grid('on');
newah=subplot(5,1,5);             % Coherence
    wlpsp_coh(time,coh,wlparam,[wlparam.freqs' R95*ones(size(coh,1))],'');
    %grid('on');

% Reposition raw data plots (alignment accounts for colourbars)
oldpos=get(ah,'pos');
newpos=get(newah,'pos');
for ind=1:length(ah)
    set(ah(ind),'pos',[oldpos{ind}(1) oldpos{ind}(2) newpos(3) oldpos{ind}(4)]);
end;

if ( exist('linkaxes','file') )
    linkaxes(get(gcf,'children'),'x');
end;
