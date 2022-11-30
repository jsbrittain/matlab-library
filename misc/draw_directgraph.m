function h = draw_directgraph( A, xy, labels, nodeProp, arrowProp )
%
%
% Input parameters
%       A           Connectivity matrix
%       xy          Node rendering locations ( N x 2{x,y} )
%       labels      Node labels
%       nodeProp    Node rendering properties
%                       Single struct or array (length N) of struct of
%                       'FaceColor', 'EdgeColor', 'EdgeWidth'
%       arrowProp   Arrow rendering properties (cell array of param-value pairs)
%                       e.g. 'FaceColor', 'EdgeColor'
%

% Input parser
p = inputParser;
addRequired( p, 'A', @isnumeric );
addRequired( p, 'xy', @isnumeric );
addOptional( p, 'labels', [] );
parse(p,A,xy,labels);

if (~exist('nodeProp'))
    nodeProp = {};
end;
if (~exist('arrowProp'))
    arrowProp = {};
end;

if (isempty('labels'))
    labels = {'1','2','3'};
end;

% Drawing properties
r = 0.35;
rscale = 3;

% Set NaN connection weights to zero
A(isnan(A))=0;

% Render
%cla;
hold('on');
axis('image','off');
h = [];
for k = (1:size(xy,1))
    h(end+1) = fillcircle( xy(k,1), xy(k,2), r, labels{k}, nodeProp{:} );
end;

% Default arrow shape
shape = [0.2,0.2,0.15,0.05];

for to = (1:size(A,1))
    for from = (1:size(A,2))
        
        if (A(to,from)<=0)
            continue;
        end;
        
        z = (xy(to,1)+1i*xy(to,2)) - (xy(from,1)+1i*xy(from,2));
        offset_perp = [ (r/rscale)*cos(3*pi/2-angle(z)) -(r/rscale)*sin(3*pi/2-angle(z)) ]; % Perpendicular
        offset_parallel = 1*[ r*cos(angle(z)) r*sin(angle(z)) ];                % Parallel
        
        h(end+1) = arrows( xy(from,1) + offset_perp(1) + offset_parallel(1), ...
                xy(from,2) + offset_perp(2) + offset_parallel(2), ...
                max([0 abs(z)-2*r]), (360/2/pi)*(pi/2-angle(z)),  ...
                [ shape([1 2 3 4])*A(to,from) ], arrowProp{:} );
            
    end;
end;
