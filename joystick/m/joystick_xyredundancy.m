function joystick=joystick_xyredundancy(joystick)
%function joyphase=joystick_xyredundancy(joystick)
%
% Remove XY redundant points
%
%function joyphase=joystick_xyredundancy(joystick)

% Reject redundant items
for n=(1:length(joystick))
    ind=2;
    while (ind<=size(joystick(n).data,1))
        % Reject trials not conforming to criteria
        if ( ( joystick(n).data((ind-1),2) == joystick(n).data(ind,2) ) & ...
             ( joystick(n).data((ind-1),3) == joystick(n).data(ind,3) ) )
                joystick(n).data=joystick(n).data([1:(ind-1) (ind+1):end],:);
        else
            % Increment counter
            ind=ind+1;
        end;
    end;
end;
