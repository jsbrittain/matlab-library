function p = wl_ttest_convolve( p, rejection_kernel, alpha )

if (~exist('alpha','var'))
    alpha = 0.05;
end;

rejection_kernel = rejection_kernel/sum(rejection_kernel(:));
cout = convolve2( double(p<alpha), rejection_kernel, 'wrap' );
p( cout < 1 ) = 1;
