function h=erpplot_lfpset(erptime,epochs,individual,confscaling,smoothing,colstr,facealpha,ylims)
%function h=erpplot_lfpset(erptime,epochs,individual,confscaling,smoothing,colstr,facealpha,ylims)
%
% Plot mono/bi-polar grid
%
% Input parameters
%
%   ylims   [ymono ybi] max
%
%function h=erpplot_lfpset(erptime,epochs,individual,confscaling,smoothing,colstr,facealpha,ylims)

% Check input parameters
if (~exist('individual'))
    individual=false;
end;
if (~exist('confscaling'))
    confscaling=1.96;
end;
if (~exist('smoothing'))
    smoothing=1;
end;
if (~exist('colstr'))
    colstr='b';
end;
if (~exist('facealpha'))
    facealpha=0.7;
end;

chlabels = { 'LFP-01', ...
             'LFP-12', ...
             'LFP-23', ...
             'LFP-0',  ...
             'LFP-1',  ...
             'LFP-2',  ...
             'LFP-3'   };

% Plot monopolars
ymono=0;
subplot(4,2,1); erpplot(erptime,epochs(:,7,:),individual,confscaling,smoothing,colstr,facealpha);
                title(chlabels{7}); ymono=max([ymono abs(ylim)]);
subplot(4,2,3); erpplot(erptime,epochs(:,6,:),individual,confscaling,smoothing,colstr,facealpha);
                title(chlabels{6}); ymono=max([ymono abs(ylim)]);
subplot(4,2,5); erpplot(erptime,epochs(:,5,:),individual,confscaling,smoothing,colstr,facealpha);
                title(chlabels{5}); ymono=max([ymono abs(ylim)]);
subplot(4,2,7); erpplot(erptime,epochs(:,4,:),individual,confscaling,smoothing,colstr,facealpha);
                title(chlabels{4}); ymono=max([ymono abs(ylim)]);
if (exist('ylims'))
    ymono=ylims(1);
end;
for ind=[1:2:7]
    subplot(4,2,ind);
    ylim(ymono*[-1 1]);
    %set(gca,'xminorgrid','on');
    %set(gca,'yminorgrid','on');
    grid('on');
    box('on');
end;

% Plot bipolars
ybi=0;
subplot(3,2,2); h=erpplot(erptime,epochs(:,3,:),individual,confscaling,smoothing,colstr,facealpha);
                title(chlabels{3}); ybi=max([ybi abs(ylim)]);
subplot(3,2,4); erpplot(erptime,epochs(:,2,:),individual,confscaling,smoothing,colstr,facealpha);
                title(chlabels{2}); ybi=max([ybi abs(ylim)]);
subplot(3,2,6); erpplot(erptime,epochs(:,1,:),individual,confscaling,smoothing,colstr,facealpha);
                title(chlabels{1}); ybi=max([ybi abs(ylim)]);
if (exist('ylims'))
    ybi=ylims(2);
end;
for ind=[2:2:6]
    subplot(3,2,ind);
    ylim(ybi*[-1 1]);
    %set(gca,'xminorgrid','on');
    %set(gca,'yminorgrid','on');
    grid('on');
    box('on');
end;

if (nargout==0)
    clear('h');
end;
