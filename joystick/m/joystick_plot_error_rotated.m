function joystick_plot_error_rotated(joyphase,leftarc,rightarc)

if (~exist('leftarc'))
    leftarc = 0;
end;
if (~exist('rightarc'))
    rightarc = 0;
end;

hold('on');

% Fill phase offset area
if ((leftarc>0) || (rightarc>0))
    xfill = [0]; yfill = [0];
    if (rightarc>0)
        xfill = [ xfill cos( (-rightarc*2*pi/360:0.01:0) + pi/2 ) ];
        yfill = [ yfill sin( (-rightarc*2*pi/360:0.01:0) + pi/2 ) ];
    end;
    if (leftarc>0)
        xfill = [ xfill cos( (0:0.01:leftarc*2*pi/360) + pi/2 ) ];
        yfill = [ yfill sin( (0:0.01:leftarc*2*pi/360) + pi/2 ) ];
    end;
    xfill = [xfill 0]; yfill = [yfill 0];
    fill(xfill,yfill,'g','linestyle','none','facealpha',0.5);
end;

% Plot single traces
for n=(1:length(joyphase))
    
    xpos = [ 0 cos( (joyphase(n).target-joyphase(n).responseangle)*2*pi/360 + pi/2 ) ];
    ypos = [ 0 sin( (joyphase(n).target-joyphase(n).responseangle)*2*pi/360 + pi/2 ) ];
    
    plot(xpos,ypos,'k');
    axis('equal');
end;
% Draw circle
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'b','linewidth',2);
plot([0 0],[0 1],'r.','markersize',20);
