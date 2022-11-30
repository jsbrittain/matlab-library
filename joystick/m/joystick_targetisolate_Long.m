function joystick=joystick_targetisolate_Long(joystick)
%function joyphase=joystick_targetisolate_Long(joystick)
%
% Isolate trials with of `target' type
%
%function joyphase=joystick_targetisolate_Long(joystick)

% Check input parameters
if (~exist('opt'))
    opt=0;
end;

ix = strcmp({joystick.type},'TARGET');
ix0 = find(ix,1,'first');
for n = ix0:(length(joystick)-1)
    joystick(n).data = [ joystick(n).data; joystick(n+1).data ];
end;
joystick = joystick(ix);

