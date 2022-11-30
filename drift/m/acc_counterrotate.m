function [txm,tym,tzm,mx,my,mz,pitch,roll] = acc_counterrotate( dat, rate, acccalib, skip )
%[txm,tym,tzm] = acc_counterrotate( mx, my, mz, acccalib )


% Default calibration data (from calib_acc.smr; note absolute values may change if biometrics amp zeroed!!!)
if (~exist('acccalib','var'))
    acccalib.xrange = [ 1.9404 2.0404 ];    % Axis range (-dir, +dir) relative to gravity (+ve upwards)
    acccalib.yrange = [ 1.9433 2.0445 ];    %  |
    acccalib.zrange = fliplr( [ 2.0943 1.9945 ] );    %  |
    % Base orientation (laying flat; z oriented upwards )
    acccalib.base   = [ mean(acccalib.xrange) mean(acccalib.yrange) acccalib.yrange(2) ];
end;
if (~exist('skip','var'))
    skip = 1;
end;

ch_ax = 1;
ch_ay = 2;
ch_az = 3;
rate = round(rate);

% Extract DC shifts (orientation changes)
mx = smooth( median_filter( dat(:,ch_ax), rate, true, skip ), rate );
my = smooth( median_filter( dat(:,ch_ax), rate, true, skip ), rate );
mz = smooth( median_filter( dat(:,ch_ax), rate, true, skip ), rate );

% Normalisation
nmx = 2*(mx-acccalib.xrange(1))/range(acccalib.xrange)-1;
nmy = 2*(my-acccalib.yrange(1))/range(acccalib.yrange)-1;
nmz = 2*(mz-acccalib.zrange(1))/range(acccalib.zrange)-1;

% Determine pitch / roll - normalisation matters!
roll = atan2( nmy, nmz );                       % Roll
pitch = atan2( nmx, sqrt(nmy.^2+nmz.^2) );      % Pitch

% Combine acc-3 after drift removal
tx = dat(:,ch_ax) - mx;
ty = dat(:,ch_ay) - my;
tz = dat(:,ch_az) - mz;

% Convert 3 tremor axes to spherical coorindates
r = sqrt( sum([ tx ty tz ].^2,2) );
theta = acos( tz ./ r );
phi = atan2( ty, tx );

% Remove Acc-3 base orientation (due to accelerometer shifts)
theta = theta + pitch;          %%% Addition done because it looked right!
phi = phi + roll;               %%% !!!

% Convert back to tremor coordinates
txm = r.*sin(theta).*cos(phi);
tym = r.*sin(theta).*sin(phi);
tzm = r.*cos(theta);
