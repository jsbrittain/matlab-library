function sp = wl_ttest2_permute( stftdata, paired, logtransform, iters )
%function sp = wl_ttest2_permute( stft, paired, logtransform, iters )
%
%
%
% Input parameters
%       stft        Spectrogram ( time x freq x trial x condition )
%                   Note that trial can be a patient average
%
% NOTE: THIS ROUTINE PROBABILITY CAN'T COPE WITH UNEVEN DATA (I.E. UNEQUAL
% TRIAL NUMBERS PER CONDITION).
%
% NOTE: THIS ROUTINE ONLY DEALS WITH 2 CONDITIONS AT A TIME FOR NOW.
%
% To be used with
%   sp_multcomp_correction.m
%
%function sp = wl_ttest2_permute( stft, paired, logtransform, iters )

% Check input parameters
if (~exist('paired','var'))
    paired = true;
end;
if (~exist('logtransform','var'))
    logtransform = [];
end;
if (~exist('iters','var'))
    iters = [];
end;

% Default parameters
if (isempty(paired))
    paired = true;
end;
if (isempty(logtransform))
    logtransform = false;
end;
if (isempty(iters))
    iters = 100;
end;

% Re-organise data
stftdata = reshape( stftdata, size(stftdata,1), size(stftdata,2), size(stftdata,3)*size(stftdata,4) );
trials = size(stftdata,3);
r = trials/2;

% Initial (unpermuted) estimate
sp = zeros( size(stftdata,1), size(stftdata,2), iters+1 );
disp([ ' Permutation 0/' num2str(iters) ]);
sp(:,:,1) = wl_ttest2( stftdata(:,:,1:r), stftdata(:,:,(r+1):end), paired );

% Permutation analysis
for iter = 2:(iters+1)
    
    % Display progress
    disp([ ' Permutation ' num2str(iter-1) '/' num2str(iters) ]);
    
    % Choose permutation metho
    if ( paired )
        % paired
        paired = true;
        groupA = (1:r) + r*(rand(1,r)>0.5);
        groupB = groupA + r;
        groupB(groupB>trials) = groupB(groupB>trials) - trials;
    else
        % unpaired
        paired = false;
        list = (1:trials);        list = list( randperm( trials ) );
        groupA = false(1,trials); groupA( list(1:r) ) = true;
        groupB = ~groupA;
    end;
    
    % Image stats
    sp(:,:,iter) = wl_ttest2( stftdata(:,:,groupA), stftdata(:,:,groupB), paired, logtransform );
    
end;
