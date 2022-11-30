function h = stackplotsc(time,dat,chlabels,gap)
%function stackplotsc([time],dat,[chlabels],gap)
%
% Stackplot on a single axis
%
%function stackplotsc([time],dat,[chlabels],gap)

% Check input parameters
if (~exist('dat'))
    dat=time;
    time = [];
end;
if (isempty(time))
    time=[0:size(dat,1)-1];
end;
if (~exist('chlabels'))
    chlabels = [];
end;
if (~exist('gap','var'))
    gap = [];
end;

% Default channel names (numbers)
if (isempty(chlabels))
    chlabels = cell(size(dat,2));
    for k = (1:size(dat,2))
        chlabels{k} = num2str(k);
    end;
end;
if (isempty(gap))
    gap = 0;
end;

% Add offset to each data channel
if (range(dat(:,1))==0)
    dat(:,1) = zeros(size(dat(:,1)));
else
    dat(:,1) = dat(:,1) - nanmin(dat(:,1));
    dat(:,1) = dat(:,1)./nanmax(dat(:,1)) + gap;
end;
for ind=(2:size(dat,2))
    if (range(dat(:,ind))==0)
        dat(:,ind) = zeros(size(dat(:,ind)));
    else
        dat(:,ind) = dat(:,ind) - min(dat(:,ind));
        dat(:,ind) = dat(:,ind)./max(dat(:,ind)) + ind - 1 + ind*gap;
    end;
    pause(0);   % just in-case we need to break!
end;

% Stackplot
%stackplot(time,dat,chlabels);
h = plot(time,dat,'k');
axis('tight'); set(gca,'ytick',[]);
title('STACK PLOT');

% Add channel labels
set(gca,'YTick',(1+gap)*(1:size(dat,2))-0.5);
set(gca,'YTickLabel',chlabels);
