function cfplot_phamp_segs_summary( params, params0, label )

if (~exist('params0','var'))
    params0=[];
end;
if (~exist('label'))
    label = '';
end;

segs = size(params.rmag,1);

% Display segments
figure; clear('h'); maxylim = []; h1 = zeros(1,segs); h2 = h1;
for k = (1:segs)
    h1(k) = subplot(2,segs+2,k);
        cfplot_phamp_segs_amp_change( params.phaseangle, params.rmag(k,:), params.perm.rmag(:,k,:) );
        title([label ' ' num2str(k)]); maxylim = max(abs([maxylim ylim]));
    h2(k) = subplot(2,segs+2,k+segs+2);
        cf_phamp_rose( params.phaseangle, params.rmag(k,:) );
        title([label ' ' num2str(k)]);
end;
% Rescale all plots
for k = (1:segs), subplot(2,segs+2,k); ylim(maxylim*[-1 1]); end;
% Average
newh = subplot(2,segs+2,segs+1);
    averageplots( h1, newh );
    ylim(maxylim*[-1 1]); title('MEAN');
newh = subplot(2,segs+2,2*segs+3);
    cf_phamp_rose( params.phaseangle, mean(params.rmag,1) );
    %averageplots( h2, newh ); axis('equal'); axis('off');
    title('MEAN');
% Separately computed average
if (~isempty(params0))
    subplot(2,segs+2,segs+2);
        cfplot_phamp_amp_change( params0.phaseangle, params0.rmag, params0.perm.rmag );
        ylim(maxylim*[-1 1]); title('GLOBAL');
    subplot(2,segs+2,2*segs+4);
        cf_phamp_rose( params0.phaseangle, params0.rmag );
        title('GLOBAL');
end;
