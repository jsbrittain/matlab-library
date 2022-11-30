function pos=electrode_trace(pos1,pos2,dist)
%function pos=electrode_trace(pos1,pos2,dist)
%
% Determine point <dist> units along electrode trajectory
%
% Originally designed for identifying stimulation contacts along an
% electrode trajectory, usually identified by post-op MR or after fusion
% with CT.
%
% Input parameters
%       pos1            Top of electrode [ x1 y1 z1 ]
%       pos2            Any other point on electrode trajectory [ x1 y1 z1 ]
%       dist            Distance locale required from tip
%
% Output
%       pos             (x,y,z) tuple of point locale
%
%function pos=electrode_trace(pos1,pos2,dist)

% Uses 2D trigonometry

% Extract input paramters
x1=pos1(1); y1=pos1(2); z1=pos1(3);
x2=pos2(1); y2=pos2(2); z2=pos2(3);

% Determine distances
x12=x2-x1;
y12=y2-y1;
z12=z2-z1;

% (x,y) locales
h2=sqrt(x12^2+y12^2);       % Euclidean distance
alpha=atan(x12/y12);        % Angle of triangle
h2=h2-dist;                 % Hypotenuse after dist subtraction
x12b=sin(alpha)*h2;         % Width of new (shortened) triangle
y12b=sqrt(h2^2-x12b^2);     % Height of new triangle

% (y,z)-locale
h2=sqrt(y12^2+z12^2);       % Euclidean distance
alpha=atan(y12/z12);        % Angle of triangle
h2=h2-dist;                 % Hypotenuse after dist subtraction
z12b=cos(alpha)*h2;         % Height of new triangle

% Output result
pos(1)=x2-x12b;
pos(2)=y2-y12b;
pos(3)=z2-z12b;
