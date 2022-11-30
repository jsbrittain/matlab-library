function newh=inaxes(x,y,width,height,h)
%
% Construct axis within confines of the current axis
%
% Input parameters
%       x           x-position
%       y           y-position
%       width       width
%       height      height
%
% All positional elements are stated in percentage terms with reference to
% the current axis
%

% Convert percentage terms to decimal
x=x/100;
y=y/100;
width=width/100;
height=height/100;

% Check for axis handle
if (~exist('h'))
    h=gca;
end;

% Determine new axis location
pos=get(h,'pos');
newx=pos(1)+x*pos(3);
newy=pos(2)+y*pos(4);
newwidth=width*pos(3);
newheight=height*pos(4);
newh=axes('Position',[newx newy newwidth newheight],...
          'Parent',get(h,'Parent'));
axes(newh);

% Lock subaxis plot to parent handle
set(h,'activepositionproperty','position');
set(newh,'activepositionproperty','position');
