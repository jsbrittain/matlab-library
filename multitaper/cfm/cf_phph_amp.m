function params = cf_phph_amp( preprocdata, trimlevel, bincount, permutation_count )
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
rmag = zeros(permutation_count+1, bincount, bincount);
rmagsd = rmag; rif = rmag; rcount = rif;

% Corresponding Phase vector
phaseangle = (0:(bincount-1))/bincount*2*pi - pi;

% Iterate through permutations
for k = (1:(permutation_count+1))
    
    % Update display
    if (k==1)
        disp('Deriving phase-phase amplitude ...')
    else
        disp([' Permutation trial ' num2str(k-1) ' of ' num2str(permutation_count) ])
    end;
    
    % Compute phase angle difference (and quantize)
    anglegridx = round( (angle( xfh )+pi)/2/pi*bincount );
    anglegridx(anglegridx==bincount) = 0;                   % Wrap first/last half-elements together
    anglegridy = round( (angle( yfh )+pi)/2/pi*bincount );
    anglegridy(anglegridy==bincount) = 0;                   % Wrap first/last half-elements together
    
    % Mean amplitude per angle
    for n1 = (0:(bincount-1))
        disp([ ' bin ' num2str(n1) ' of ' num2str(bincount) ]);
        for n2 = (0:(bincount-1))
            rmag(k,n1+1,n2+1)   = trimmean( magxh((anglegridx==n1) & (anglegridy==n2)), trimlevel );
            rmagsd(k,n1+1,n2+1) = std( magxh((anglegridx==n1) & (anglegridy==n2)) );
            rif(k,n1+1,n2+1)    = trimmean( ifxh((anglegridx==n1) & (anglegridy==n2)), trimlevel );
            rcount(k,n1+1,n2+1) = sum((anglegridx==n1) & (anglegridy==n2));
        end;
    end;

    % Phase shuffle predictor for next iteration
    yfh = yfh( randperm(length(yfh)) );
    
end;

% Populate parameters structure
params.phaseangle = phaseangle;
params.rmag = squeeze( rmag(1,:,:) );
params.rmagsd = squeeze( rmagsd(1,:,:) );
params.rif = squeeze( rif(1,:,:) );
params.rcount = squeeze( rcount(1,:,:) );
params.perm.rmag = rmag(2:end,:,:);
params.perm.rmagsd = rmagsd(2:end,:,:);
params.perm.rif = rif(2:end,:,:);
params.perm.rcount = rcount(2:end,:,:);
params.settings.amp_norm_secs = amp_norm_secs;
params.settings.permutation_count = permutation_count;
params.settings.bincount = bincount;

% Display update
disp('done.');
