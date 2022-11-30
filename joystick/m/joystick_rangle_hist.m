function [p,bins,moments]=joystick_rangle_hist(rangle,histgrid)

% Check input parameters
if (~exist('histgrid'))
    histgrid=[];
end;

% Default values
if (isempty(histgrid))
    histgrid=(-180:5:180);
end;

% Construct
[p,bins]=hist(rangle,histgrid);
p=p/sum(p);     % Convert to PDF

% Determine central moments
N = length(bins); dt=bins(2)-bins(1);
moments(1) = sum( bins.*p );      % Mean
moments(2) = sum( ((bins-moments(1)).^2).*p );       % Variance
moments(3) = sum( (((bins-moments(1))/(moments(2)^(1/2))).^3).*p ); % Skew
