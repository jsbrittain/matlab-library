function [ smoothing_kernel, x0, y0, z0 ] = kernel_3dgauss( varargin )
%smoothing_kernel = kernel_3dgauss( x, y, z, sd )
%
% Construct 3D Gaussian kernel
%
% Input parameters
%       x       (Nx-vec) x-grid points
%       y       (Ny-vec) y-grid points
%       z       (Nz-vec) z-grid points
%       sd      (scalar) Std dev for isotropic Gaussian kernel
%
% Outputs a smoothing_kernel (3D grid) on the provided grid set yet
%  restricted to the 95% limits to reduce kernel size and aid subsequent
%  convolution operations.
%
% Notes:
%   Assumes an isotropic kernel and regular, equally spaced grid points
%    in x, y & z.
%   If all inputs are absent (or empty) then script runs in demo mode.
%
%smoothing_kernel = kernel_3dgauss( x, y, z, sd )

% Default parameters
verbose            = true;
display_kernel     = false;
sd_extents         = 2;
display_isosurface = false;

% Input parser
p = inputParser;
addOptional( p, 'x', [] );
addOptional( p, 'y', [] );
addOptional( p, 'z', [] );
addOptional( p, 'sd', [] );
addParameter( p, 'verbose', verbose );
addParameter( p, 'display_kernel', display_kernel );
addParameter( p, 'sd_extents', sd_extents );
addParameter( p, 'Sigma', [] );
addParameter( p, 'display_isosurface', display_isosurface );
parse(p,varargin{:});
% Read-out values
x = p.Results.x;
y = p.Results.y;
z = p.Results.z;
sd = p.Results.sd;
verbose = p.Results.verbose;
display_kernel = p.Results.display_kernel;
sd_extents = p.Results.sd_extents;
Sigma = p.Results.Sigma;
display_isosurface = p.Results.display_isosurface;

% Default parameters
if ( (isempty(x)) || (isempty(y)) || (isempty(z)) )
    x = (-10:0.1:10);
    y = x; z = x;
    display_kernel = true;
end;
if (isempty(sd))
    sd = 1.0;
end;

% Check if std dev isotropic (scalar) or vector-valued
if (~isempty(Sigma))
    % Full covariance matrix specified
    if ~all( size(Sigma) == [ 3 3 ] )
        error('Covariance matrix badly specified, should be 3 x 3.');
    end;
    sd(1) = sqrt( Sigma(1,1) );
    sd(2) = sqrt( Sigma(2,2) );
    sd(3) = sqrt( Sigma(3,3) );
    if ( verbose ), disp('Using specified full covariance matrix...'); end;
else
    if (isscalar(sd))
        sd = repmat(sd,1,3);
    end;
    assert(length(sd)==3,'Parameter `sd` not scalar (isotropic), 3-vector (anisotropic).');
    Sigma = diag( sd.^2 );
end;

% Re-centre and restrict
x0 = x - mean(x); x0 = x0( (x0 > -sd_extents*sd(1)) & (x0 < sd_extents*sd(1)) );
y0 = y - mean(y); y0 = y0( (y0 > -sd_extents*sd(2)) & (y0 < sd_extents*sd(2)) );
z0 = z - mean(z); z0 = z0( (z0 > -sd_extents*sd(3)) & (z0 < sd_extents*sd(3)) );

% Construct kernel on existing grid
[XX,YY,ZZ] = meshgrid(x0,y0,z0);        % Re-centre grid; outputs ( Ny x Nx x Nz )
XX = permute( XX, [ 2 1 3 ] );          %  |
YY = permute( YY, [ 2 1 3 ] );          %  |
ZZ = permute( ZZ, [ 2 1 3 ] );          %  |
smoothing_kernel = zeros(length(x0),length(y0),length(z0));
if ( verbose ), progress = 0; fprintf('Constructing kernel over %g points (%g x %g x %g) from (%g x %g x %g) grid, sd = (%.1f, %.1f, %.1f)\n <',numel(XX),length(x0),length(y0),length(z0),length(x),length(y),length(z),sd); end;
for k = (1:numel(XX))
    if ( verbose )
        newprogress = round(100*k/numel(XX));
        if ( newprogress > progress )
            progress = newprogress;
            if (( mod(progress,10) == 0 ) && ( progress ~= 100 ))
                fprintf('%g',progress);
            else
                fprintf('.');
            end;
        end;
    end;
    smoothing_kernel(k) = mvnpdf( [ XX(k) YY(k) ZZ(k) ], [0 0 0], Sigma );
end;
smoothing_kernel = smoothing_kernel/sum(smoothing_kernel(:));
if (verbose), fprintf('>\n'); end;

% Plot kernel
if ( display_kernel )
    image3d_ortho( x0, y0, z0, smoothing_kernel );
    if ( display_isosurface )
        subplot(224);
            cla; hold('on');
            cmap = lines(3);
            p = patch( isosurface( permute(XX,[2 1 3]), permute(YY,[2 1 3]), permute(ZZ,[2 1 3]), permute(smoothing_kernel,[2 1 3]), quantile( smoothing_kernel(:), 0.75 ) ) );
            p.FaceColor = cmap(1,:); p.EdgeColor = 'None'; p.FaceAlpha = 0.4;
            p = patch( isosurface( permute(XX,[2 1 3]), permute(YY,[2 1 3]), permute(ZZ,[2 1 3]), permute(smoothing_kernel,[2 1 3]), quantile( smoothing_kernel(:), 0.50 ) ) );
            p.FaceColor = cmap(2,:); p.EdgeColor = 'None'; p.FaceAlpha = 0.6;
            p = patch( isosurface( permute(XX,[2 1 3]), permute(YY,[2 1 3]), permute(ZZ,[2 1 3]), permute(smoothing_kernel,[2 1 3]), quantile( smoothing_kernel(:), 0.25 ) ) );
            p.FaceColor = cmap(3,:); p.EdgeColor = 'None'; p.FaceAlpha = 0.8;
            
            axis('image'); xlabel('x'); ylabel('y'); zlabel('z');
            daspect([1 1 1]);
            view(3); axis('tight');
            camlight; lighting('gouraud');
    end;
end;
