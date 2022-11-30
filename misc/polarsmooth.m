function h = polarsmooth( z, smoothing, opt_str )
%
% z is complex
%
% Parameters
%       z           Complex valued vectors
%       smoothing   Filter width in samples (non-uniform) or degrees (if interpolating)
%       opt_str     Options string
%                       m    Min subtract (subtract min magnitude from population)
%                       n    Normalise max magnitude to unity
%                       d    Decimal places for plot magnitude rounding
%                       c    Circle radius
%                               Special options (NaN-none, 0-fit to mean smoothed magnitude)
%                               Radius scaled appropriately for `n' and `m' options;
%                               if you do not want scaling state radius as -ve value
%                       i<n> Smooth over a (uniform) interpolated grid (samples have non-uniform spacing otherwise)
%                               n = interpolation method (0-none, 1-nearest, 2-linear, 3-spline)
%                               Interpolated to 360o, smoothing parameter then specified in degrees
%

% Check input parameters
if (~exist('smoothing','var'))
    smoothing = [];
end;
if (~exist('opt_str','var'))
    opt_str = '';
end;

% Default parameters
if (isempty(smoothing))
    smoothing = 10;
end;
minmagsubtract = false;         % Subtract min magnitude to maximize dynamic range of plot
maxmagnorm = false;             % Normalise max magnitude to unity
magdps = 0;                     % Decimal places to round magnitude axis to
interpmethod = 0;               % Smooth over a (uniform) interpolated grid
circleradius = 0;               % Circle radius (NaN = none)

% Parse options string
options=deblank(opt_str); opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'm'                    % Min magnitude subtract
        	minmagsubtract = true;
        case 'n'                    % High magnitude normalise (normalise max mag to unity)
        	maxmagnorm = true;
        case 'd'                    % Decimal places for plot magnitude rounding
        	magdps = str2num(optarg);
            if (isempty(magdps))
                error(['Error in option argument -- ' opt]);
            end;
        case 'c'                    % Circle radius
        	circleradius = str2num(optarg);
            if (isempty(circleradius))
                error(['Error in option argument -- ' opt]);
            end;
        case 'i'                    % Smooth over a (uniform) interpolated grid
        	interpmethod = str2num(optarg);
            if (isempty(interpmethod))
                error(['Error in option argument -- ' opt]);
            end;
        otherwise                   % Options for wavelet analysis
            error(['Unknown option -- ' opt]);
    end;
end;

% Subtract the min magnitude (helps use the whole dynamic range of the plot)
if ( minmagsubtract )
    minabsz = min(abs(z));
    z = (abs(z)-minabsz).*exp(1i*angle(z));
else
    minabsz = 0;
end;
% Normalise max magnitude to unity
if ( maxmagnorm )
    maxabsz = max(abs(z));
    z = (abs(z)/maxabsz).*exp(1i*angle(z));
else
    maxabsz = 1;
end;

% Sort data by angle
[~,ix] = sort(angle(z));
zsorted = z( ix );

% Interpolate data before smoothing
if ( interpmethod > 0 )
    switch ( interpmethod )
        case 1, method = 'nearest';
        case 2, method = 'linear';
        case 3, method = 'spline';
        otherwise
            error(' Unknown interpolation method specified.');
    end;
    newangles = (-pi:(2*pi/360):pi);
    newmag = interp1( angle(zsorted), abs(zsorted), newangles, method );
    zsortedi = newmag.*exp(1i*newangles);
else
    zsortedi = zsorted;
end;

% Smooth data (circularly)
smoothz = smooth( repmat(zsortedi(:),3,1), smoothing );
smoothz = smoothz((end/3+1):2*end/3);
% Sort smoothed data
[~,ix] = sort(angle(smoothz));
smoothzsorted = smoothz( ix );

% Plot figure
%figure;
if (isnan(magdps))
    delete(polar(max(abs(z)))); % Use rose to format plot
else
    delete(polar(ceil(max(abs(z))*(10^-magdps))*(10^magdps)));   % Use rose to format plot
end;
hold('on');
h(1) = plot( zsorted([1:end 1]), 'b.-' );
h(2) = plot( smoothzsorted([1:end 1]), 'r.-', 'linewidth', 2 );

% Plot circle
if (~isnan(circleradius))
    if (circleradius==0)
        % Fit circle to smooth vector
        circleradius = mean(abs(smoothz));
    elseif (circleradius>0)
        circleradius = (circleradius-minabsz)/maxabsz;
    else
        circleradius = -circleradius;
    end;
    h(3) = plot( circleradius*cos(-pi:0.01:pi), circleradius*sin(-pi:0.01:pi), 'k--' );
end;
