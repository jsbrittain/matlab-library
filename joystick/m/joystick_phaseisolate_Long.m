function joystick=joystick_phaseisolate_Long(joystick,phaseoffset,opt)
%function joyphase=joystick_phaseisolate_Long(joystick,phaseoffset,[opt])
%
% Isolate trials with specific phase componenets
%
% opt   0   all trials (default)
%       1   target trials
%       2   fixation trials
%
%function joyphase=joystick_phaseisolate_Long(joystick,phaseoffset,[opt])

% Check input parameters
if (~exist('opt'))
    opt=0;
end;

% Reject trial type
switch (opt)
    
    case 1,     % Target trials only
        joystick = joystick_targetisolate_Long(joystick);
        
    case 2,     % Fixation trials only
        joystick = joystick_fixationisolate(joystick);
        
end;


% Reject non-specific phase offsets
n=1;
while (n<=length(joystick))
    % Reject trials not conforming to criteria
    if (joystick(n).phaseoffset~=phaseoffset)
        joystick=joystick([1:(n-1) (n+1):end]);
    else
        % Increment counter
        n=n+1;
    end;
end;
