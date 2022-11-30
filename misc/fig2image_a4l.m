function fig2image_a4l(fh, filename, imageformat);
% function fig2image_a4l(fh, filename, imageformat);
%
% Print a figure specified by its handle (fh) to a jpeg file specified by
% filename.  Print setup is configured to full page A4 Landscape before
% print commences.
%

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

% Print to jpeg
print(fh, filename, ['-d' imageformat]);
