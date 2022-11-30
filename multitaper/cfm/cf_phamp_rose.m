function cf_phamp_rose( phaseangle, rmag )
% Rose plot

smoothing = 5;

x0 = phaseangle;
y0 = (rmag-min(rmag))/(max(rmag)-min(rmag));

x = interp1((1:length(x0)),x0,(1:0.1:length(x0)));
y = interp1((1:length(x0)),y0,(1:0.1:length(x0)));

h0=polar( x([end 1:end]), y([end 1:end]), 'b' ); set(h0,'linewidth',1,'color',[0.4 0.8 1.0]); hold('on');
    x=get(h0,'xdata'); y=get(h0,'ydata');
fill([x zeros(size(x))],[y zeros(size(y))],'b','linestyle','none','facecolor',[0.4 0.8 1.0],'facealpha',0.5); hold('on');
h=polar( x0([end 1:end]), smooth( (y0([end 1:end])-min(y0))/(max(y0)-min(y0)), smoothing)', 'r');
    set(h, 'linewidth', 2 );
h=polar( x0([end 1:end]), (mean(y0)-min(y0))/(max(y0)-min(y0))*ones(1,length(x0)+1), 'k');
    set(h, 'linewidth', 2 );
