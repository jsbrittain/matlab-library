function [rangle,joy] = joystick_responseangle_correct_Muriel_SV( rangle, joy, filterpts100, startoffset )
%function [rangle,joy] = joystick_responseangle_correct( rangle, joy, filterpts100, startoffset )
%
%
%
%function [rangle,joy] = joystick_responseangle_correct( rangle, joy, filterpts100, startoffset )

if (~exist('filterpts100'))
    filterpts100 = [];
end;
if (isempty(filterpts100))
    filterpts100 = 4;
end;
if (isempty(startoffset))
    startoffset = 0;
end;

% Correct individual trials
figure; n = 1;
if (true), set(gcf,'WindowStyle','docked'); figure(gcf); end;
while ( n <= length(joy) )
    
    % Variables for response angle calculation
    startpos = mean(joy(n).data((1:5)+startoffset,(2:3)),1);
    
    
    % Adjust filter length
    filterpts = max(1,ceil(filterpts100*joy(1).rate/100));
    b=ones(1,filterpts)/filterpts; a=1;
    
    % Variables for plotting only
    dist = sqrt(sum(joy(n).data(:,2:3).^2,2));
    vel = abs(diff(dist,[],1));
    fvel = filtfilt(b,a,vel);
    if (joy(n).reject), statusstr = ' - REJECTED'; else statusstr = ''; end;
    
    % Find max velocity point and determine baseline - max-vel gradient (degrees)
    if (~isfield(joy,'rtrange'))
        joy(n).rtrange = [];
    end;
    if (~isfield(joy,'rtvelpt'))
        joy(n).rtvelpt = [];
    end;
    if (~isfield(joy,'responseangle'))
        joy(n).responseangle = NaN;
    end;
    if (~isempty(joy(n).rtvelpt))
        rangle(n) = mod( 90 - (atan2(joy(n).data(joy(n).rtvelpt,3)-startpos(2),joy(n).data(joy(n).rtvelpt,2)-startpos(1))*360/2/pi), 360);
        joy(n).responseangle = rangle(n);
        rangle(n) = joy(n).responseangle - joy(n).target;
        if ( rangle(n) > 180 )
            % Wrap response vector
            rangle(n) = rangle(n) - 360;
        end;
    end;
    
    % Display trial
    clf;
    subplot(3,3,[1 2 4 5]); box('on'); hold('on');
        plot(0,0,'ko');
        plot(joy(n).data(:,2),joy(n).data(:,3),'b');
        plot(startpos(1),startpos(2),'gx');
        plot(joy(n).data(joy(n).rtvelpt,2),joy(n).data(joy(n).rtvelpt,3),'rx');
        plot([0 cos((90-joy(n).responseangle)/360*2*pi)],[0 sin((90-joy(n).responseangle)/360*2*pi)],'r');
        plot(cos((90-joy(n).target)/360*2*pi),sin((90-joy(n).target)/360*2*pi),'ko');
        axis('equal'); xlim([-1 1]*1.3); ylim([-1 1]*1.3);
        if (joy(n).reject), set(gca,'color',[1.0 0.8 0.8]); else set(gca,'color',[0.8 1.0 0.8]); end;
        title([num2str(n) ' (TARGET ' num2str(joy(n).target) '^o, ACHIEVED ' num2str(joy(n).responseangle) '^o)' statusstr]);
    subplot(3,3,3); box('on');
        if (~isempty(vel))
            plot(vel);
        end; hold('on');
        if (~isempty(fvel))
            plot(fvel);
        end;
        if (~isempty(joy(n).rtrange))
            plot(joy(n).rtrange(1)*[1 1],ylim,'g');
            plot(joy(n).rtrange(2)*[1 1],ylim,'g');
        end;
        if (~isempty(joy(n).rtvelpt))
            plot(joy(n).rtvelpt,fvel(joy(n).rtvelpt),'ro');
        end;
        xlim([0 125]);
    
    % Get user choice and respond
    choice = input( ' continue (ENTER), back (b), peak velocity (p), search range (s), reject (r), filter length (f), save & exit (x) ? ', 's' );
    if (isempty(choice))         % ENTER
        n = n + 1;
    elseif strcmpi(choice,'p')   % Adjust peak velocity
        disp('Select peak velocity ...');
        [x,y] = ginput(1);
        joy(n).rtvelpt = round(x);
    elseif strcmpi(choice,'s')   % Correct search range ( start / end points )
        disp('Select start point ...');
        [startpt,y] = ginput(1); startpt = round( startpt );
        disp('Select end point ...');
        [endpt,y] = ginput(1); endpt = round( endpt );
        joy(n).rtrange = [startpt endpt];
        % Auto-select max velocity between adjusted points
        [maxvel,velpt] = max(vel(startpt:endpt));
        velpt = velpt + startpt - 1;
        joy(n).rtvelpt = velpt;
        joy(n).maxvel = maxvel;
    elseif strcmpi(choice,'r')   % Toggle rejection
        joy(n).reject = (~joy(n).reject);
    elseif strcmpi(choice,'b')   % Back
        n = n - 1;
        if (n<1)
            disp('  YOU ARE AT THE FIRST TRIAL. Type ''x'' to exit.');
            n = 1;
        end;
    elseif strcmpi(choice,'f')   % Change filter length
        disp(' WARNING: This changes the filter length on ALL TRIALS.');
        filterpts100 = str2num( input(['  Enter new filter length (currently ' num2str(filterpts100) ' pt): '],'s') );
    elseif strcmpi(choice,'x')   % Save and exit
        break;
    else
        disp([' Unknown option: ' choice]);
    end;
    
    
     
end;

% Fill response angle for rejected trials
rangle([joy.reject]) = nan;

% Clear figure
close(gcf);
