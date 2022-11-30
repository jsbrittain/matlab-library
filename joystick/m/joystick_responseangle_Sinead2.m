function [rangle,joystick]=joystick_responseangle_Sinead2( joystick, filterpts100 )
%function [rangle,joystick]=joystick_responseangle( joystick, filterpts100 )
%
%
% Calculate error angle by finding the maximum velocity point of a trace
% and calculating the angle from that point to the origin.
%
%function [rangle,joystick]=joystick_responseangle( joystick, filterpts100 )

if (~exist('filterpts100'))
    filterpts100 = 4;
end;

% Parameters
minreactiontime = 50;       % msecs
maxRadius = 1.0;%0.2;
%filterpts100 = 15;           % Filter points (for 100 Hz, rescaled later)
crossingTh = 0.06;

% Recurse joystick data determining response angle
rangle=nan(length(joystick),1);
for n=(1:length(joystick))
        
    % By default, reject
    joystick(n).reject = true;
    
    % Check sample rate (per trial)
    if (isfield(joystick(n),'rate'))
        filterpts = max(1,ceil(filterpts100*joystick(n).rate/100));
    end;
    
    % Determine baseline range
    time = joystick(n).data(:,1) - joystick(n).data(1,1);
    startpos = mean(joystick(n).data((1:5),(2:3)),1);
    dist = sqrt(sum(joystick(n).data(:,2:3).^2,2));
    vel = abs(diff(dist,[],1));
    
    % Reject trials beginning outside of a 0.2 radius from the origin
    if (~isempty(find(dist(1:5)>maxRadius,1,'first')))
        continue;
    end;
    
    % Filter and extract smoothed velocity statistics
    b=ones(1,filterpts)/filterpts; a=1;
    try
        fvel = filtfilt(b,a,vel);
    catch
        disp(['Excluding trial ' num2str(n) ' (filter error)']);
        continue;
    end;
    startpt = find( (fvel > (max(fvel)/14)) & (time(1:end-1)>minreactiontime), 1, 'first' );
    endpt = find( (vel(startpt+3:end) < max(vel)/12) , 1, 'first' ) + startpt + 2 ; %& (time(startpt:end-1)<time(end-10)) & (time(startpt:end-1)>time(startpt+2))
   %  endpt = find( fvel(startpt:end) < (max(fvel)/4), 1, 'first' ) +
   %  startpt - 1;
    if (isempty(endpt)), endpt = length(vel); end;
    
    % Find max velocity point and determine baseline - max-vel gradient (degrees)
    [maxvel,velpt] = max(vel(startpt:endpt)); %[maxvel,velpt] = max(fvel(startpt:endpt));
    velpt0 = velpt + startpt - 1;
    velpt = find(vel(startpt:endpt) > crossingTh,1,'first') + startpt - 1;
    if (isempty(velpt))
        velpt=velpt0;
    end;
    % Reject velocity points inside of a 0.1 radius from the origin (false starts)
%     if (dist(velpt)<0.1)
%         continue;
%     end;
    
if (isempty(velpt))
    joystick(n).reject = true;
    joystick(n).responseangle = NaN;
    joystick(n).rt = NaN;
    joystick(n).rtrange = [NaN NaN];
    joystick(n).rtvelpt = [];
    joystick(n).maxvel = NaN;
    rangle(n) = NaN;
else
    % Record values
    joystick(n).reject = false;
    joystick(n).responseangle = mod( 90 - (atan2(joystick(n).data(velpt,3)-startpos(2),joystick(n).data(velpt,2)-startpos(1))*360/2/pi), 360);
    joystick(n).rt = time(velpt);
    joystick(n).rtrange = [startpt endpt];
    joystick(n).rtvelpt = velpt;
    joystick(n).maxvel = maxvel;
    rangle(n) = joystick(n).responseangle - joystick(n).target;
end
    
    if ( rangle(n) > 180 )
        % Wrap response vector
        rangle(n) = rangle(n) - 360;
    end;
    
end;
