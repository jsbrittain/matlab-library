function params = cf_phamp_seg_mag( preprocdata, trimlevel, bincount, permutation_count, segs )
%
%
% Most preprocessing and amplitude analysis done on 'x', whereas 'y'
% provides the phase comparitor.
%
%

% Expand proprocessed data to (local) workspace
struct2vars( preprocdata );

% Reset severe deflections to mean and smooth
ifxh(ifxh>(2*passbandx(2))) = mean(ifxh);
ifxh = smooth(ifxh,smoothing);

%% Analyse TACS data, then randomly permute phase and continue re-analysing

% Reset phase
yfh = hilbert(yf);
yfh = smooth(yfh, smoothing);

% Reserve memory
rmag = zeros(permutation_count+1, segs, bincount);
rmagsd = rmag; rif = rmag; rcount = rif;

% Corresponding Phase vector
phaseangle = (0:(bincount-1))/bincount*2*pi - pi;

% Determine magnitude level (quantised into `segs' portions)
maglevel = zeros(length(magxh),1);
[~,ix] = sort(magxh);
iix = floor( (0:(length(ix)-1))*segs/length(ix) );
maglevel(ix) = iix + 1;

% Iterate through permutations
for ind = (1:(permutation_count+1))
    
    % Update display
    if (ind==1)
        disp('Deriving phase-amplitude ...')
    else
        disp([' Permutation trial ' num2str(ind-1) ' of ' num2str(permutation_count) ])
    end;
    
    % Compute phase angle difference (and quantize)
    anglediff = angle( xfh./yfh );
    anglegrid = round( (anglediff+pi)/2/pi*bincount );
    anglegrid(anglegrid==bincount) = 0;                     % Wrap first/last half-elements together

    % Recurse segments
    for k = (1:segs)
        
        % Mean amplitude per angle
        for n = (0:(bincount-1))
            rmag(ind,k,n+1) = trimmean( magxh((anglegrid==n) & (maglevel==k)), trimlevel );
            rmagsd(ind,k,n+1) = std( magxh((anglegrid==n) & (maglevel==k)) );
            rif(ind,k,n+1) = trimmean( ifxh((anglegrid==n) & (maglevel==k)), trimlevel );
            rcount(ind,k,n+1) = sum((anglegrid==n) & (maglevel==k));
        end;

    end;
    
    % Phase shuffle predictor for next iteration
    yfh = yfh( randperm(length(yfh)) );
    
end;

% Populate parameters structure
params.phaseangle = phaseangle;
params.rmag = squeeze(rmag(1,:,:));
params.rmagsd = squeeze(rmagsd(1,:,:));
params.rif = squeeze(rif(1,:,:));
params.rcount = squeeze(rcount(1,:,:));
params.perm.rmag = rmag(2:end,:,:);
params.perm.rmagsd = rmagsd(2:end,:,:);
params.perm.rif = rif(2:end,:,:);
params.perm.rcount = rcount(2:end,:,:);
params.settings.amp_norm_secs = amp_norm_secs;
params.settings.permutation_count = permutation_count;
params.settings.bincount = bincount;

% Display update
disp('done.');
