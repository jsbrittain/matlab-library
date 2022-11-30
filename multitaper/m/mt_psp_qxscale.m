function h=mt_psp_qxscale(params,sp11,sp22,xscale);
%function h=mt_psp_qxscale(params,sp11,sp22,xscale);
%
% Plot routine for use with multitaper analysis
%
% Plot cumulant density on its full time-base then overlay with an
% enlargment over [xscale(1) xscale(2)].  Useful in presenting global
% activity and primary peak detail in a single plot.
%
% Input parameters
%       params      Parameters structure (contains cumulant details)
%       sp11        Spectrum ch.1 (used for confidence limits)
%       sp22        Spectrum ch.2 (used for confidence limits)
%       xscale      (opt) Enlargment x-axis (msecs), i.e. [-50 50]
%                         Default: [], no enlargment axis
%
% Output parameter
%       h           (opt) Axes handles (2)
%
%function h=mt_psp_qxscale(params,sp11,sp22,xscale);

% Check input parameters
if ((nargin<3) | (nargin>4))
    error(' Incorrect number of input parameters');
end;
if (~exist('xscale'))
    xscale=[];
end;

% Plot cumulant on full-time base
mt_psp_q(params,sp11,sp22);
if (isempty(xscale)), return; end;
h(1)=gca; ylims=ylim; box('off'); title('');
set(h(1),'yticklabel',{});

% Overlay second axis over cumulant plot
h(2)=axes('position',get(h(1),'position'), ...
          'xaxislocation','top', ...
          'yaxislocation','right', ...
          'color','none');
line(params.qlags,params.q,'color','k','linewidth',2);
xlim(xscale); ylim(ylims);
box('off');
set(h(2),'yticklabel',{});

% Maintain alignment when rescaling figures (i.e. for printing)
align(h,'distribute','distribute');
