%%

addpath ~/matlab/library/misc/
addpath ~/matlab/library/joystick/m/

filename = { 'Vera_baseline.txt',
             'Vera_training.txt',
             'Vera_consolidation.txt',
             'Vera_aftereffect.txt' };
phaseoffset = [0 60 60 0];

removezeros = true;
resample = 1000;
for ind=(1:length(filename))
    joystick{ind} = joystick_load(filename{ind},removezeros,resample);
end;

%%

clear rangle

selection = 1;      % Targets

for ind=(1:length(filename))
    joyphase{ind} = joystick_phaseisolate(joystick{ind},phaseoffset(ind),selection);
    [rangle{ind},joyphase{ind}]=joystick_responseangle(joyphase{ind});
end;

%% Plot angular error

group=16;

figure; fig2a4p;
ylims = [-75 75];
maxx = 0;

for ind=(1:length(filename))
    subplot(length(filename),1,ind);
        joystick_plot_rangle_err(rangle{ind},group); ylim(ylims);
        maxx = max(maxx,length(rangle{ind}));
        title(strrep(filename{ind}(1:end-4),'_','\_'));
end;
for ind=(1:length(filename))
    subplot(length(filename),1,ind);
    xlim([0 maxx+1]);
end;

%% Plot all trials

figure;
ylims = [-1.5 1.5];
for ind=(1:length(joyphase))
    
    joy=joyphase{ind};
    trialstart = 1;
    
    triallist = (trialstart:length(joyphase{ind}));
    
    subplot(length(joyphase),4,(ind-1)*4+1);
        joystick_plot_raw(joy(triallist));
        ylim(ylims); ylabel(strrep(filename{ind}(1:end-4),'_','\_'));
    
    subplot(length(joyphase),4,(ind-1)*4+2);
        joystick_plot_raw_rotated(joy(triallist),60,60);
        ylim(ylims);
        
    subplot(length(joyphase),4,(ind-1)*4+3);
        joystick_plot_error_rotated(joy(triallist),60,60);
        ylim(ylims);
    
    subplot(length(joyphase),4,(ind-1)*4+4);
        [p,bins,moments]=joystick_plot_rangle_hist(rangle{ind}(triallist));
        xlabel(['\mu=' num2str(moments(1)) ', \sigma=' num2str(sqrt(moments(2))) ', \lambda_1=' num2str(moments(3))]);
        
end;
subplot(length(joyphase),4,1); title('RESPONSE');
subplot(length(joyphase),4,2); title('ROTATED RESPONSE');
subplot(length(joyphase),4,3); title('ROTATED ERROR');
subplot(length(joyphase),4,4); title('ERROR HISTOGRAM');

%% Plot individual trial

clear rangle

joy=joyphase{1};

figure;
for n=(1:length(joy))
    
    % Determine baseline range
    startpos = mean(joy(n).data((1:5),(2:3)),1);
    dist = sqrt(sum(joy(n).data(:,2:3).^2,2));
    vel = abs(diff(dist,[],1));
    
    % Filter and extract smoothed velocity statistics
    N=5; b=ones(1,N)/N; a=1; fvel = filtfilt(b,a,vel);
    startpt = find( fvel > (max(fvel)/4), 1, 'first' );
    endpt = find( fvel(startpt:end) < (max(fvel)/4), 1, 'first' ) + startpt - 1;
    if (isempty(endpt)), endpt = length(vel); end;

    % Find max velocity point and determine baseline - max-vel gradient (degrees)
    [dummy1,velpt] = max(fvel(startpt:endpt));
    velpt = velpt + startpt - 1;
    rangle(n) = mod( 90 - (atan2(joy(n).data(velpt,3)-startpos(2),joy(n).data(velpt,2)-startpos(1))*360/2/pi), 360);
    
    clf;
    subplot(3,3,[1 2 4 5]); hold('on');
        plot(0,0,'ko');
        plot(joy(n).data(:,2),joy(n).data(:,3),'b');
        plot(startpos(1),startpos(2),'gx');
        plot(joy(n).data(velpt,2),joy(n).data(velpt,3),'rx');
        plot([0 cos((90-rangle(n))/360*2*pi)],[0 sin((90-rangle(n))/360*2*pi)],'r');
        plot(cos((90-joy(n).target)/360*2*pi),sin((90-joy(n).target)/360*2*pi),'ko');
        axis('equal'); xlim([-1 1]*1.3); ylim([-1 1]*1.3);
        title([num2str(n) ' (TARGET ' num2str(joy(n).target) '^o, ACHIEVED ' num2str(rangle(n)) '^o)']);
    subplot(3,3,3);
        plot(vel); hold('on'); plot(fvel);
        plot(startpt*[1 1],ylim,'g');
        plot(endpt*[1 1],ylim,'g');
        plot(velpt,fvel(velpt),'ro');

    pause;
end;
