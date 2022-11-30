function out = runavg( data, windowSize )
%
% This function is more efficient for long windowSize arguments than
% Matlab's 'filter' routine.
%
%

if ( isvector(data) )
    data = data(:);
else
    out = zeros(size(data));
    for k = (1:size(data,2))
        out(:,k) = runavg(data(:,k),windowSize);
    end
    return;
end

if (windowSize == 1 )
    out = data;
    return;
end

% Halve averaging window since we apply this filter fwd then bwd
windowSize = floor(windowSize/2);

% Sum and subtract (forward)
out = cumsum(data,'omitnan')/windowSize;
out(windowSize+1:end) = out(windowSize+1:end) - out(1:end-windowSize);
% Correct weighting on edge (beginning) samples
out(1:windowSize) = (out(1:windowSize)*windowSize)./(1:windowSize)';

% Sum and subtract (backwards)
data=rot90(out,2);                    % Flip running-average
out = cumsum(data,'omitnan')/windowSize;
out(windowSize+1:end) = out(windowSize+1:end) - out(1:end-windowSize);
% Correct weighting on edge (end) samples
out(1:windowSize) = (out(1:windowSize)*windowSize)./(1:windowSize)';

% Flip back
out = rot90(out,2);
