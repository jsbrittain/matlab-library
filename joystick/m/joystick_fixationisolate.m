function joystick=joystick_fixationisolate(joystick)
%function joyphase=joystick_fixationisolate(joystick)
%
% Isolate trials with of `fixation' type
%
%function joyphase=joystick_fixationisolate(joystick)

% Check input parameters
if (~exist('opt'))
    opt=0;
end;

% Reject non-specific phase offsets
n=1;
while (n<=length(joystick))
    % Reject trials not conforming to criteria
    if (~strcmp(joystick(n).type,'FIXATION'))
        joystick=[joystick(1:(n-1)) joystick((n+1):end)];
    else
        % Increment counter
        n=n+1;
    end;
end;
