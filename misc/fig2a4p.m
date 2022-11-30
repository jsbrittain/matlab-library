function fig2a4p(fh);
% function fig2a4p(fh);
%
% Print setup is configured to full page A4 Portrait.
%

if nargin<1
    fh=gcf;
end;

% Fill print page for A4 Portrait (default: cm)
top = 0.63;
left = 0.63;
width = 19.72;
height = 28.41;

for ind=1:length(fh)
	set(fh(ind), 'PaperOrientation', 'portrait');
	punits=get(fh(ind), 'PaperUnits');
	set(fh(ind), 'PaperUnits', 'centimeters');
	rect = [left, top, width, height];
	set(fh(ind), 'PaperPosition', rect);
    set(fh(ind), 'PaperUnits', punits);
end;
