function image3d_ortho( x, y, z, I, varargin )

% Defaults
projection = 'nanmean';
flipx = false;
flipy = false;
flipz = false;

% Layout specific defaults (also prototypes)
layout = [];
switch ( 'grid' )
    case 'grid',
        layout.type = 'grid';
        layout.grid = [ 2 2 ];
        layout.gridpos = [ 1 2 3 ];
    case 'col',
        layout.type = 'col';
        layout.ncols = 6;
        layout.colpos = 1;
    case 'row',
        layout.type = 'row';
        layout.nrows = 4;
        layout.rowpos = 1;
end;

% Input parser
p = inputParser;
addRequired( p, 'x', @isnumeric );
addRequired( p, 'y', @isnumeric );
addRequired( p, 'z', @isnumeric );
addRequired( p, 'I', @isnumeric );
addParameter( p, 'projection', projection );
addParameter( p, 'flipx', flipx );
addParameter( p, 'flipy', flipy );
addParameter( p, 'flipz', flipz );
addParameter( p, 'layout', layout );
parse(p,x,y,z,I,varargin{:});
% Read-out values
projection = p.Results.projection;
flipx = p.Results.flipx;
flipy = p.Results.flipy;
flipz = p.Results.flipz;
layout = p.Results.layout;

% Check if axes specified
if (isempty(x)), x = (1:size(I,1)); end;
if (isempty(y)), y = (1:size(I,2)); end;
if (isempty(z)), z = (1:size(I,3)); end;

% Select layout method
figure;
switch ( layout.type )
    case 'grid',
        ah(1) = subplot(layout.grid(1),layout.grid(2),layout.gridpos(1));
        ah(2) = subplot(layout.grid(1),layout.grid(2),layout.gridpos(2));
        ah(3) = subplot(layout.grid(1),layout.grid(2),layout.gridpos(3));
    case 'col',
        ah(1) = subplot(3,layout.ncols,layout.colpos);
        ah(2) = subplot(3,layout.ncols,layout.colpos+layout.ncols);
        ah(3) = subplot(3,layout.ncols,layout.colpos+2*layout.ncols);
end;

% Select projection method
switch ( projection )
    case 'mean',    op = @(x,dim) squeeze(mean(x,dim));
    case 'nanmean', op = @(x,dim) squeeze(nanmean(x,dim));
    case 'max',     op = @(x,dim) squeeze(max(x,[],dim));
    case 'nanmax',  op = @(x,dim) squeeze(nanmax(x,[],dim));
    case 'min',     op = @(x,dim) squeeze(min(x,[],dim));
    case 'nanmin',  op = @(x,dim) squeeze(nanmin(x,[],dim));
end;

% Display
axes(ah(1));
    imagesc(x,z,op(I,2).');
    axis('image');
    xlabel('x'); ylabel('z');
    if ( flipx ), set(gca,'xdir','reverse'); end;
    if ( flipz ), set(gca,'ydir','normal'); end;
axes(ah(2));
    imagesc(y,z,op(I,1).');
    axis('image');
    xlabel('y'); ylabel('z');
    if ( flipy ), set(gca,'xdir','reverse'); end;
    if ( flipz ), set(gca,'ydir','normal'); end;
axes(ah(3));
    imagesc(x,y,op(I,3).');
    axis('image');
    xlabel('x'); ylabel('y');
    if ( flipx ), set(gca,'xdir','reverse'); end;
    if ( flipy ), set(gca,'ydir','normal'); end;

% Link axes
linkaxes(ah([1 2]),'y');
linkaxes(ah([1 3]),'x');
