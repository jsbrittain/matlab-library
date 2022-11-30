
[X,Y,Z] = sphere(10);
X = [ X(:) Y(:) Z(:) ];

X = randn(1e2,5);

[h,y,v1] = plotn( X, '.' );
axis('image');
