function h = mt_psp_q(params,sp11,sp22);
%function mt_psp_q(params,[sp11,sp22]);
%
% Cumulant plotting routine for use with
%   mt_sp, mt_sp2, mt_sp2_pp
%
% Plot multitaper estimates with jackknife errors if available
%
%function mt_psp_q(params,[sp11,sp22]);

% Determine common parameters
if (params.jackknife)
    % Spectral confidence limits
    p=0.05;
    ttest=0;%tq(1-p/2,params.jkcount-1);
end;
% Handles for conf. lines
h1=[];

% Check if label specified
if (~isfield(params,'what'))
    what='';
else
    what=[': ' params.what];
end;

% Calculate empirical 95% conf. limit (based on estimated spectra - Neurospec only)
q95=0;
T=params.duration*params.rate/1000; R=params.jkcount*T;
if (params.mtparams.spec_norm/(params.duration*params.rate/1000)==2*pi)
    qnorm=(2*pi)^2/(T*R);
elseif (params.mtparams.spec_norm==1)
    qnorm=1/(R*T);
else
    warning(' Unrecognised spectral norm for cumulant plotting');
    qnorm=0;
end;
if (params.jackknife)
    q95=1.96*sqrt(qnorm*sum(2*mean(params.jk11,2).*mean(params.jk22,2)));
end;

% Plot cumulant
h1 = [];
if (params.jackknife)
    warning(' Cumulant density jackknife not applicable.');
%     if (params.bootstrap)
%         % h1(end+1)=plot(params.qlags,params.q+1.96*sqrt(params.bsqv)); hold('on');
%         % h1(end+1)=plot(params.qlags,params.q-1.96*sqrt(params.bsqv));
%         h1(end+1)=plot(params.qlags,0-1.96*sqrt(params.bsqv)); hold('on');
%         h1(end+1)=plot(params.qlags,0+1.96*sqrt(params.bsqv));
%     else
%         % h1(end+1)=plot(params.qlags,params.q+1.96*sqrt(params.jkqv)); hold('on');
%         % h1(end+1)=plot(params.qlags,params.q-1.96*sqrt(params.jkqv));
%         h1(end+1)=plot(params.qlags,0+1.96*sqrt(params.jkqv)); hold('on');
%         h1(end+1)=plot(params.qlags,0-1.96*sqrt(params.jkqv));
%     end;
end;
% Empirical confidence limits
if (isempty(find(diff(params.L))) & exist('sp11') & exist('sp22')) % If all L the same
    qnorm=1/R/(R*params.L(1));
    c95=1.96*sqrt(qnorm*sum((real(params.fmax~=0)+1)*sp11.*sp22));
    h = plot(params.qlags([1 end]),-c95*[1 1],'k',params.qlags([1 end]),c95*[1 1],'k'); hold on;
end;
% Plot cumulant
plot(params.qlags([1 end]),[0 0],'k--'); hold on;
if (params.jackknife)
    h = [ h; plot(params.qlags([1 end]),-q95*[1 1],'k',params.qlags([1 end]),q95*[1 1],'k') ];
end;
h = [ plot(params.qlags,params.q,'k'); h ]; hold off;
axis('tight');
xlims=xlim; xlim(max(abs(xlims))*[-1 1]);
xlabel('OFFSET (MSECS)');
title(['q' what]);

% Format confidence limits
set(h1,'color',0.6*[1 1 1]);       % Solid gray line
