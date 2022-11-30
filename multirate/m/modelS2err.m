function f = modelS2err( z, attained, attainedvar, target, clamptrials )
% Multi-rate model

% Run model
x = modelS2( z, target, clamptrials );

% Chi-square goodness-of-fit (set var=1 for sum-of-squares)
f = nansum( (( attained - x).^2)./attainedvar );
