function z = sortangle( z )

[~,ix] = sort( angle(z) );
z = z(ix);
