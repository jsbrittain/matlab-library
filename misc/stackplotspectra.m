function h = stackplotspectra(time,dat,fmax,chlabels)
%function stackplotspectra([time],dat,fmax,[chlabels])
%
% Stackplot on a single axis
%
%function stackplotspectra([time],dat,fmax,[chlabels])

% Check input arguments
if (~exist('chlabels'))
    chlabels = {};
end;
if (~exist('fmax'))
    fmax = [];
end;
if (~exist('dat'))
    dat = time;
    time = [];
end;
if (isempty(time))
    time = (1:size(dat,1));
end;

% Frequency analysis
duration = 1000; rate = 1/(time(2)-time(1)); opt_str = 'W0.5';
for ind = (1:size(dat,2))
    [sp11(:,ind),~,~,params] = mt_sp0(dat(:,ind),dat(:,ind),duration,rate,opt_str);
end;

% Construct figure
figure;
subplot(1,8,(1:7));
    stackplotsc( time, dat, chlabels );
subplot(1,8,8);
    stackplotsc( params.freqs, sp11, chlabels );
    set(gca,'yticklabel',[]);
    if (~isempty(fmax)), xlim([0 fmax]); end;
    