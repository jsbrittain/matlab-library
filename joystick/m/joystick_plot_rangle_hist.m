function [p,bins,moments]=joystick_plot_rangle_hist(rangle,histgrid)

% Check input parameters
if (~exist('histgrid'))
    histgrid=[];
end;

% Display histogram
[p,bins,moments] = joystick_rangle_hist(rangle,histgrid);

bar(bins,p);
hold('on');
plot([0 0],ylim,'r');
xlim(bins([1 end]));
