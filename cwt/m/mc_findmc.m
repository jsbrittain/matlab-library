function cl=mc_findmc(wlparam);
%
% Find appropriate confidence limits for any specified analysis from the
% monte-carlo store.
%

% Load confidence limits structure
load('mc_results');
mc_index=0;
mc_trials=0;
for ind=1:length(mc_results)
    if (mc_results{ind}.wlparam.linearscale==0), res='dj'; else res='df'; end;
    if ((mc_results{ind}.wlparam.linearscale==wlparam.linearscale) & ...
        (mc_results{ind}.wlparam.dt==wlparam.dt) & ...
        (length(mc_results{ind}.wlparam.coi)==length(wlparam.coi)) & ...
        (mc_results{ind}.wlparam.(res)==wlparam.(res)) & ...
        (mc_results{ind}.trials>mc_trials))
            mc_index=ind;
    end;
end;
if (mc_index==0)
    disp(' Warning: No montecarlo simulations have been performed with the specified analysis parameters.');
    disp('          Coherence figures will not include 95% confidence limits.');
    cl=[];
else
    cl=mc_results{mc_index};
    if (~isfield(cl,'freqs'))
        cl.freqs=1./(cl.wlparam.fourier_factor*cl.scale);
    end;
end;
