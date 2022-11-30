function preprocdata = cf_preproc( x, y, passbandx, passbandy, dropoff, rate, amp_norm_secs, xMagThreshold, smoothing, filterchoice )
%
%
% Most preprocessing and amplitude analysis done on 'x', whereas 'y'
% provides the phase comparitor.
%
%

% Check input parameters
if (~exist('filterchoice','var'))
    filterchoice = [];
end;

% Default parameters
if (isempty(passbandy))
    passbandy = passbandx;
end;
if (isempty(smoothing))
    smoothing = 1;
end;
if (isempty(filterchoice))
    filterchoice = 0;
end;

%% Preprocessing

% Update display
%disp( 'Preprocessing...' );

if (~isempty(passbandx))
    
    switch ( filterchoice )
        
        case 0,         % Bespoke Kaiser FIR linear phase filter
            % Specify filter
            d = fdesign.bandpass;
            d = fdesign.bandpass(passbandx(1)-dropoff,passbandx(1),passbandx(2),passbandx(2)+dropoff,d.Astop1,d.Apass,d.Astop2,rate);
            h = design(d,'kaiserwin');      % Kaiser window FIR filter (linear phase)
            Num = h.Numerator;
            lNum = length(Num);
            % Filter first channel
            xf = filtfilt(Num,1,x);
            xf = xf - mean(xf);
            % Filter TACS channel with secondary filter
            if (~isempty(find(passbandx~=passbandy, 1)))
                % Specify filter
                d = fdesign.bandpass(passbandy(1)-dropoff,passbandy(1),passbandy(2),passbandy(2)+dropoff,d.Astop1,d.Apass,d.Astop2,rate);
                h = design(d,'kaiserwin');      % Kaiser window FIR filter (linear phase)
                Num = h.Numerator;
                lNum = max([lNum length(Num)]);
            end;
            % Filter second channel
            yf = filtfilt(Num,1,y);
            yf = yf - mean(yf);
            
        case 1,         % Butterworth (IIR) filter
            [bh,ah] = butter(8,passbandx(1)/(rate/2),'high');
            [bl,al] = butter(8,passbandx(2)/(rate/2),'low');
            xf = filtfilt(bh,ah,filtfilt(bl,al,x));
            xf = xf - mean(xf);
            % Filter TACS channel with secondary filter
            if (~isempty(find(passbandx~=passbandy, 1)))
                [bh,ah] = butter(8,passbandy(1)/(rate/2),'high');
                [bl,al] = butter(8,passbandy(2)/(rate/2),'low');
            end;
            yf = filtfilt(bh,ah,filtfilt(bl,al,y));
            yf = yf - mean(yf);
            lNum = 2*length(bh);
            
    end;
    
    % Chop edges
    xf = xf(lNum:end-lNum);
    yf = yf(lNum:end-lNum);
else
    % No filtering
    xf = x;
    yf = y;
end;

% Hilbert and smooth Ax
xfh = hilbert(xf);
xfh = smooth(xfh, smoothing);

% Accelerometer magnitude
magxh = abs(xfh);

% Clip out data where magnitude of `x' does not meet threshold
if (~isempty(xMagThreshold))
    % Clip x-vector variables
    x   =   x(magxh >= xMagThreshold);
    xf  =  xf(magxh >= xMagThreshold);
    xfh = xfh(magxh >= xMagThreshold);
    % Clip y-vector variables
    y   =   y(magxh >= xMagThreshold);
    yf  =  yf(magxh >= xMagThreshold);
    % Clip magnitude itself
    magxh = magxh(magxh > xMagThreshold);
end;

% Running-average amplitude normalisation
magxh_orig = magxh;
if (~isempty(amp_norm_secs))
    switch (0)      %%% NB: BOTH IDENTICAL %%%
        case 0,         % Normalise envelope
            magxh = magxh./runavg(magxh,amp_norm_secs*rate);
        case 1,         % Normalise data based on envelope
            xf = xf./runavg(magxh,amp_norm_secs*rate);
            xfh = hilbert(xf);
            xfh = smooth(xfh, smoothing);
            magxh = abs(xfh);
    end;
end;

% Accelerometer (instantaneous frequency)
ifxh = diff(angle(xfh)); ifxh = ifxh([(1:end) end]);
ifxh(ifxh<0) = ifxh(ifxh<0) + 2*pi;
ifxh = ifxh*rate/2/pi;
% Reset severe deflections to mean and smooth
%ifxh(ifxh>(2*passbandx(2))) = mean(ifxh);
%ifxh = smooth(ifxh,smoothing);

%% Collapse (local function) workspace into structure and output

vnames = who;
preprocdata = struct;

for n=(1:length(vnames))
    preprocdata.(vnames{n}) = eval(vnames{n});
end;
