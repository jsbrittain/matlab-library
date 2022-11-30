function joystick_plot_rangle_err(rangle,group)

if (~exist('group','var'))
    group = [];
end;
if (isempty(group))
    group=16;
end;
scaling = 1;

hold('on');
% plot(rangle,'k'); xlim([0 length(rangle)]); plot(xlim,[0 0],'k');
t=(1:group:length(rangle));
for k=(1:length(t))
    tt=t(k)+(0:(group-1)); tt=tt(tt<length(rangle));
    h=plot(mean(tt),nanmean(rangle(tt)),'ro');
    %h=plot(mean(tt)+[-0.5 0.5],nanmean(rangle(tt))*[1 1],'r');
    h(end+1) =  plot(mean(tt)*[1 1],nanmean(rangle(tt))+scaling*[-1 1]*nanstd(rangle(tt)),'r');
    set(h,'linewidth',2);
end;
