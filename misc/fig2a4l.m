function fig2a4l(fh);
% function fig2a4l(fh);
%
% Print setup is configured to full page A4 Landscape
%

if nargin<1
    fh=gcf;
end;

% Fill print page for A4 Portrait (default: cm)
top = 0.63;
left = 0.63;
width = 28.41;
height = 19.72;

% rect = [left, bottom, width, height]
set(fh, 'PaperOrientation', 'landscape');
punits=get(fh, 'PaperUnits');
set(fh, 'PaperUnits', 'centimeters');
rect = [left, top, width, height];
set(fh, 'PaperPosition', rect);
set(fh, 'PaperUnits', punits);
