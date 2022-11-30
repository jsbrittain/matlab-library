function add_middle_xaxis

hold on;
plot(xlim,[0 0],'k','linewidth',1);
tickpos = get(gca,'xtick');
ticklength = get(gca,'ticklength');
for n = (1:length(tickpos))
    plot(tickpos(n)*[1 1],0.02*range(ylim)*[-1 1],'k','linewidth',1);
end;
set(gca,'xtick',[]);
