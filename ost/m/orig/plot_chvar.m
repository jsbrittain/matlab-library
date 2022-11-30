function plot_chvar(f,x,varx);
% function plot_chvar(f,x,varx);
% Plot coherence with variance bounds

limits=1.96*sqrt(varx);
plot(f,x,'k-',f,x-limits,'k:',f,x+limits,'k:');
