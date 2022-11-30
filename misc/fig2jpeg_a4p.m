function fig2jpeg_a4p(fh, filename);
% function fig2jpeg_a4p(fh, filename);
%
% Print a figure specified by its handle (fh) to a jpeg file specified by
% filename.  Print setup is configured to full page A4 Portrait before
% print commences.
%

% Fill print page for A4 Portrait (default: cm)
top = 0.63;
left = 0.63;
width = 19.72;
height = 28.41;

% rect = [left, bottom, width, height]
punits=get(fh, 'PaperUnits');
set(fh, 'PaperUnits', 'centimeters');
rect = [left, top, width, height];
set(fh, 'PaperPosition', rect);
set(fh, 'PaperUnits', punits);

% Print to jpeg
print(fh, filename, '-djpeg');
