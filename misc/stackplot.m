function [h,yticks] = stackplot(time,dat,chlabels)
%function stackplot([time],dat,[chlabels])
%
% Stackplot on a single axis
%
%function stackplot([time],dat,[chlabels])

% Check input parameters
if (~exist('dat'))
    dat=time;
    time=[0:size(dat,1)-1];
end;
if (~exist('chlabels'))
    chlabels = [];
end;

% Default channel names (numbers)
if (isempty(chlabels))
    chlabels = cell(size(dat,2));
    for k = (1:size(dat,2))
        chlabels{k} = num2str(k);
    end;
end;

dat(isnan(dat)) = 0;

% Add offset to each data channel
dat(:,1)=dat(:,1)-min(dat(:,1));
for ind=(2:size(dat,2))
    dat(:,ind) = dat(:,ind)-min(dat(:,ind)) + max(dat(:,ind-1));
    pause(0);   % just in-case we need to break!
end;

% Produce stackplot
h = plot(time,dat,'k');
axis('tight'); set(gca,'ytick',[]);
title('STACK PLOT');

% Add channel labels
if (~isempty(chlabels))
    yticks = min(dat,[],1) + (max(dat,[],1)-min(dat,[],1))/2;
    set(gca,'YTick',yticks);
    set(gca,'YTickLabel',chlabels);
end;
