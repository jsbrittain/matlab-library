function conv_kernel = conv_kernel_circle( n )

if (mod(n,2)==0)
    error('Circular convolution kernel must be of odd-size');
end;

% Circular convolution kernel
conv_kernel = zeros( n );
width = ceil(size(conv_kernel,1)/2);
height = ceil(size(conv_kernel,2)/2);
[ii,jj] = ind2sub( size(conv_kernel), 1:numel(conv_kernel) );
dist = sqrt((ii-width).^2 + (jj-height).^2);
conv_kernel( dist<width ) = 1;
