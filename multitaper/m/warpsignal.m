function z = warpsignal(v, lambda)
%function z = warpsignal(v, lambda)
%
%
%
%function z = warpsignal(v, lambda)

if ( lambda == 0  )
    z = v;
    return;
end;

% Warp signals
bap = [-lambda 1]; aap = [1 -lambda];                               % All-pass filter
ir = ( filter([-lambda 1],[1 -lambda],[1; zeros(size(v,1),1)]) );   % Compute impulse response
H = fft(ir, size(v,1));                                             % FFT for circular convolution
z0 = v; z = z0;

progress = 0;
display_progress = false;

for n = (2:size(v,1))         % Samples
    
    % Determine time range for trial
    if (display_progress)
        if (floor(n/size(v,1)*100)>progress)
            progress=floor(n/size(v,1)*100);
            disp(['Multi-taper trial ' int2str(progress) '% (' int2str(size(v,1)) ' samples)']);
            pause(0);
        end;
    end;
    
    for k = (1:size(v,2))       % Channels
        switch ( 2 )
            case 1,     % Linear convolution
                z0(:,k) = filter(bap,aap,z0(:,k));
            case 2,     % Circular convolution
                z0(:,k) = real(ifft(fft(z0(:,k)).*H));
        end;
        z(n,k) = z0(end,k);
    end;
end;
