function polarRing( scale, varargin )

inner_ring = 1.0;
outer_ring = 3.0;
colormapstr = 'parula';
resdeg = 1;     % 1 deg resolution
normalise = true;
linewidth = 2;
common_props = {};
direction = 1;      % +1 = anticlockwise, -1 for clockwise
offset = 0;

p = inputParser;
addRequired( p, 'scale', @isnumeric );
addOptional( p, 'inner_ring', inner_ring, @isnumeric );
addOptional( p, 'outer_ring', outer_ring, @isnumeric );
addOptional( p, 'colormap', colormapstr, @isnumeric );
addOptional( p, 'res', resdeg, @isnumeric );
addOptional( p, 'normalise', normalise, @isscalar );
addOptional( p, 'linewidth', linewidth, @isnumeric );
addOptional( p, 'common_props', common_props, @iscell );
addOptional( p, 'direction', direction, @isnumeric );
addOptional( p, 'offset', offset, @isnumeric );
parse(p,scale,varargin{:});

inner_ring = p.Results.inner_ring;
outer_ring = p.Results.outer_ring;
colormapstr = p.Results.colormap;
resdeg = p.Results.res;
normalise = p.Results.normalise;
linewidth = p.Results.linewidth;
common_props = p.Results.common_props;
direction = p.Results.direction;
offset = p.Results.offset;

if ( normalise )
    scale = scale/nanmax(scale);
end;

[~,ix] = sort( scale );
[~,order] = sort(ix);

N = length(scale);
res = 2*pi*resdeg/360;

cmap = eval([colormapstr '(N)']);
colormap(cmap);
facecolor = mat2cell( cmap, ones(N,1), 3 );
common_properties = { 'linewidth', linewidth, common_props{:} };

hold('on');
for n = (1:N)
    theta = direction * ( offset*2*pi/360 + (n-1)*2*pi/N + (0:res:1)*2*pi/N );
    x = [ inner_ring*cos(theta) (inner_ring + scale(n)*(outer_ring-inner_ring))*cos(fliplr(theta)) inner_ring*cos(theta(1)) ];
    y = [ inner_ring*sin(theta) (inner_ring + scale(n)*(outer_ring-inner_ring))*sin(fliplr(theta)) inner_ring*sin(theta(1)) ];
    fill( x, y, 'r', 'facecolor', facecolor{order(n)}, common_properties{:} );
end;

axis('image','off');
xlim(outer_ring*[-1 1]);
ylim(outer_ring*[-1 1]);
