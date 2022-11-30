function joystick=joystick_load_Angel(filename,removezeros,resample)
%function joystick=joystick_load(filename,[removezeros,[resample]])
%
% Load joystick trajectories for LabVIEW motor adaptation task coded by
% JSB (11/06/09).
%
%
% `resample' rate should be multiples of 1 msec (ie. 200, 500, 1000 Hz)
%
%function joystick=joystick_load(filename,[removezeros,[resample]])

% Parameters
if (~exist('removezeros','var'))
    removezeros = true;
end;
if (~exist('resample','var'))
    resample = [];
end;

% Initialise data structure
joystick=[];
index=0;

% Open file
fid=fopen(filename,'r');

% Check filename
if (fid==-1)
    error(' File not found');
end;

% Recurse file
disp(' Parsing file ...');
while (~feof(fid))
    
    % Read line
    line=fgetl(fid);
    
    % Blank line - ignore
    if (isempty(line))
        continue;
    end;
    
    % Comment - ignore
    if (strcmp(line(1),'#'))
        continue;
    end;
    
    % Event message
    if (strcmp(strtok(line),'MSG'))
        
        % Increment trial index
        index=index+1;
        
        % Parse message
        remain=line;
        for n=(1:5)
            [token{n},remain]=strtok(remain);
        end;
        % Populate structure
        joystick(index).type=token{2};
        joystick(index).phaseoffset=str2num(token{5});
        switch (joystick(index).type)
            case 'TARGET'
                joystick(index).target=str2num(token{3});
            otherwise
                joystick(index).target=0;
        end;
        joystick(index).data=[];
        
        % Continue loop
        continue;
    end;
    
    % Data
    numeric=str2num(line);
    if (length(numeric)==3)
        % Append data to structure
        joystick(index).data=[joystick(index).data; numeric];
    else
        % Unknown file contents
        warning([' Unknown file contents: ' line]);
    end;

end;
disp(' File loaded.');

% Remove trailing zeros
if (removezeros)
    for n=(1:length(joystick))
        joystick(n).data = joystick(n).data(joystick(n).data(:,1)>0,:);
    end;
end;

% Remove duplicate time entries
for n=(1:length(joystick))
    duplicates = find(diff(joystick(n).data(:,1))==0);
    joystick(n).data = joystick(n).data(setdiff((1:size(joystick(n).data,1)),duplicates+1),:);
end;

% Remove
for n=(1:length(joystick))
    k = find(diff(joystick(n).data(:,1))<0,1,'first');
    if (~isempty(k))
        joystick(n).data = joystick(n).data(1:(k-1),:);
        warning([' NON-CAUSAL TIME ENTRY (n=' num2str(n) ') !!!!!!!']);
    end;
end;

% Resample data
if (~isempty(resample))
    dt = round(1000/resample);
    for n=(1:length(joystick))
        joystick(n).rate = resample;
        xi = (1:dt:max(joystick(n).data(:,1))).';
        % Linear interpolate (prevents wild spline oscillations later)
        datainterp1 = interp1(joystick(n).data(:,1),joystick(n).data(:,2:3),xi,'linear');   % No start points
        firstNonNan = find(~isnan(joystick(1).data(:,2)),1,'first');                        % NaNs at beginning
        % Spline interpolate (smoothes linear interpolation and provides start / end points)
        datainterp1 = interp1(joystick(n).data(firstNonNan:end,1),joystick(n).data(firstNonNan:end,2:3),xi,'spline');
        joystick(n).data = xi;
        joystick(n).data(:,2:3) = datainterp1;
    end;
end;

% Close file
fclose(fid);
