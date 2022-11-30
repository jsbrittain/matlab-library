function cf_phamp_rose_conf( phaseangle, rmag, rsd, count )
% Rose plot

% Extract metric of interest
x0 = phaseangle;
y0 = (rmag-min(rmag))/(max(rmag)-min(rmag));
z0 = (rsd-min(rmag))/(max(rmag)-min(rmag));

% Interpolate
x = interp1((1:length(x0)),x0,(1:0.1:length(x0)),'');
y = interp1((1:length(x0)),y0,(1:0.1:length(x0)),'');
z = interp1((1:length(x0)),z0,(1:0.1:length(x0)),'');

% Plot
h0=polar( x([end 1:end]), y([end 1:end]), 'b' );
    set(h0,'linewidth',1,'color',[0.4 0.8 1.0]); hold('on');
    x=get(h0,'xdata'); y=get(h0,'ydata');
fill([x zeros(size(x))],[y zeros(size(y))],'b','linestyle','none','facecolor',[0.4 0.8 1.0],'facealpha',0.5); hold('on');
h=polar( x0([end 1:end]), smooth( (y0([end 1:end])-min(y0))/(max(y0)-min(y0)), 5)', 'r');
    set(h, 'linewidth', 2 );
h=polar( x0([end 1:end]), (mean(y0)-min(y0))/(max(y0)-min(y0))*ones(1,length(x0)+1), 'k');
    set(h, 'linewidth', 2 );

% Plot confidence intervals
polar( x0, y0-1.96*z0/sqrt(count), 'r' );
polar( x0, y0+1.96*z0/sqrt(count), 'r' );
