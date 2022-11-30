function mt_psp_hold(sp11,sp22,sp12,params,maxf);
%function mt_psp_hold(sp11,sp22,sp12,params,[maxf]);
%
% Spectral plotting routine for use with
%   mt_sp, mt_sp2, mt_sp2_pp
%
% Plot multitaper estimates with jackknife errors if available
%
%function mt_psp_hold(sp11,sp22,sp12,params,[maxf]);

% Plot conf. limits then asymptotic / null confidence lines then statistic
% in order to preserve statistic as the top level display.

% Determine number of cols
if (isfield(params,'q'))
    cols=3;
else
    cols=2;
end;

% Determine maximum plotting frequency
if (~exist('maxf'))
    maxf=params.rate/2;
end;

% Plot summary figure
%figure;
subplot(cols,2,1); hold('on'); mt_psp_a(sp11,params,['f' num2str(maxf)]);              % Auto-spectra ch.1
subplot(cols,2,2); hold('on'); mt_psp_a(sp22,params,['f' num2str(maxf)]);              % Auto-spectra ch.2
subplot(cols,2,3); hold('on'); mt_psp_ph(sp12,params,maxf,sp11,sp22);   % Phase
subplot(cols,2,4); hold('on'); mt_psp_coh(sp11,sp22,sp12,params,maxf,'p1');  % Coherence

% Cumulant density (if available)
if (isfield(params,'q'))
	subplot(3,2,5); mt_psp_q(params,sp11,sp22);
end;
