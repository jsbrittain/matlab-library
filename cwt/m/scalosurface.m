function scalogram(time,period,scale,power,coi);

% Check inputs
if (nargin<4)
    error('Not enough input parameters.');
end;
if (nargin>4)
    plot_coi=1;
else
    plot_coi=0;
end;

% Contour plot the wavelet transform with a Fourier frequency y-axis
surface(time,log2(period),power);
ax(1)=gca;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
set(ax(1),'YLim',log2([min(period) max(period)]),'YTick',log2(Yticks(:)),'YTickLabel',Yticks,...
    'XLim',[min(time) max(time)],'ygrid','on', 'xgrid','on');
ylabel('Frequency (Hz)');

% Add COI
if (plot_coi)
    % COI=0 for first and last points
    hold on;
    plot(time(2:(length(time)-1)),log2(1./(coi(coi>0))),'k');
end;

view([-149 54]);
