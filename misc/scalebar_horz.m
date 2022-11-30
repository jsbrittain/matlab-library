function scalebar_horz(x,y,width,label);
%function scalebar_horz(x,y,width,label);
%
% Add a horizontal scale bar to the current axis
%
% Input parameters
%       x           Centre x co-ordinate
%       y           Centre y co-ordinate
%       width       Width of bar (in x-axis units)
%       label       Text label (displayed centralised below scalebar)
%
%function scalebar_horz(x,y,width,label);

% Prepare to overlay scalebar on plot
holdon=ishold;
hold('on');
% Draw scalebar (horizontally)
plot(x+[-1 1]*width/2,y*[1 1],'color','k','linewidth',3);
% Add label (string includes an empty line to shift label below scalebar)
text(x,y,{'',label},'HorizontalAlignment','center');
% Return axis to previous hold condition
if (~holdon)
    hold('off');
end;
