function x = fft_comb( x, rate, fc, width )
%function fft_comb( x, rate, fc, width )
%
%function fft_comb( x, rate, fc, width )

if (~exist('width','var'))
    width = [];
end
if (isempty(width))
    width = 1;
end

fcc = (fc:fc:floor(rate/2-width-1));
freqs = repmat(fcc',1,2)+repmat(width*[-1 1],length(fcc),1);
for ind = (1:size(x,2))
    x(:,ind) = fft_notch( x(:,ind), rate, freqs );
end
