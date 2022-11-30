function joystick_plot_raw(joyphase)

targetcount = 8;

hold('on');
for n=(1:length(joyphase))
    plot(joyphase(n).data(:,2),joyphase(n).data(:,3),'k');
    axis('equal');
end;
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'b','linewidth',2);
plot([0 0],[0 0],'r.','markersize',20);
plot(cos([0:(2*pi/targetcount):(2*pi-0.01)]),sin([0:(2*pi/targetcount):(2*pi-0.01)]),'r.','markersize',20);
