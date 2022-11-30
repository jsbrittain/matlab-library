addpath(genpath('/home/jsb/Dropbox/Matlab/library'));
addpath(genpath('c:/users/jsb/Dropbox/Matlab/library'));


% Simulation parameters
secs = 100;
dt = 1e-3;
rate = 1/dt;
init = 2;

% Common encoder parameters
tau = 2.5e-2;
th = 1;

% Synchronisation level
synchrony = 0.1;

N = 100;
unsync_cells = (1:N);
sync_cells = zeros(1,floor(synchrony*N));
fb_cells = [];%(N+1);
for k = (1:floor(synchrony*N))
    kk = randi(length(unsync_cells),1,1);
    sync_cells(k) = unsync_cells(kk);
    unsync_cells = unsync_cells([1:(kk-1) (kk+1):end]);
end;

muc = 1.18; sigmac = 0.1;
mui = 0;%.1;%1.02;
sigmai = 0;%.8;%0.065;

% Run
T = 10*rate;
v = zeros(N+1,T); x = zeros(N+1,1);
v(:,1) = rand(size(v(:,1)));
v(:,100) = rand(size(v(:,1)));
firing = false(N+1,T);
progress = 0; fprintf('\n Sim: ');
fired = false(N+1,1);
for t = (1:T-1)
    
    % Change-point
    if (t==T/2)
        muc = 2.4; sigmac = 0.1;
    end;
    
    % Display progress
    if (round(100*t/(T-1))>progress)
        progress=round(100*t/(T-1));
        %disp(['    ' int2str(progress) '% complete.']);
        if ( mod(progress,10)==0 )
            fprintf('%g',progress);
        else
            fprintf('.');
        end;
    end;
    % Independent / synchronised input
    x(unsync_cells) = muc+sigmac*randn(length(unsync_cells),1);
    x(sync_cells)   = muc+sigmac*randn(1);
%     % Feedback cell
%     x(fb_cells)     = sum( 0.1*fired );
%     if (firing(fb_cells,t-100))     % Firing with delay
%         x(1:N) = x(1:N) + (1/N);
%     end;
    % Additional random input
    x(1:N) = x(1:N) + mui+sigmai*randn(N,1);
    % Euler forward step (1ms)
    v(:,t+1) = (v(:,t)+(dt/tau)*x)/(1+dt/tau);
    % Firing times
    fired = v(:,t+1) > th;
    v(fired,t+1) = v(fired,t+1) - th;
    firing( fired, t+1 ) = true;
end;
output = mean( firing, 1 );
fprintf(' %%\n');

figure;
subplot(5,5,[2:5 7:10 12:15 17:20]);
    hold on;
    for k = (1:N)
        if (ismember(k,sync_cells))
            colstr = 'r.';
        else
            colstr = 'k.';
        end;
        plot( find(firing(k,:)), k, colstr );
    end;
    title(['Mean firing rate: ' num2str(mean(mean(firing,2)*rate)) ' Hz']);
    colormap(flipud(gray(2)));
subplot(5,5,22:25);
    %plot( smooth( output, 200 ) );
    plot( smooth( mean(firing(unsync_cells,:),1), 250 ), 'k' ); hold on
    plot( smooth( mean(firing(sync_cells,:),1), 250 ), 'r' );
subplot(5,5,21);
    if false && (exist('mt_sp0','file'))
        [sp11,sp22,sp12,params]=mt_sp0(output,output,1000,rate,'W1');
        plot( params.freqs, 10*log10(sp11), 'k' );
    else
        Hs = spectrum.welch('hamming',2048);
        psd(Hs,output,'Fs',rate);
        %xlim([1 50]);
    end;
