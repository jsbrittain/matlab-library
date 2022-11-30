function add_middle_yaxis

hold on;
plot([0 0],ylim,'k','linewidth',1);
tickpos = get(gca,'ytick');
ticklength = get(gca,'ticklength');
for n = (1:length(tickpos))
    plot(0.0075*range(xlim)*[-1 1],tickpos(n)*[1 1],'k','linewidth',1);
end;
set(gca,'ytick',[]);
