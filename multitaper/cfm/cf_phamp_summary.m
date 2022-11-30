function cf_phamp_summary( params, pvalue )

% New figure
figure;

% Mean amplitude
subplot(321);
    cfplot_phamp_amp( params.phaseangle, params.rmag, params.perm.rmag, params.settings.amp_norm_secs, pvalue );

% Change in amplitude about the mean
subplot(222);
    cfplot_phamp_amp_change( params.phaseangle, params.rmag, params.perm.rmag, pvalue );

% Incidence
subplot(323);
    cfplot_phamp_amp_incidence( params.phaseangle, params.rcount, params.perm.rcount, pvalue )

% Rose plot of amplitude
subplot(224);
    polar( params.phaseangle([end 1:end]),        (params.rmag([end 1:end])-min(params.rmag))/(max(params.rmag)-min(params.rmag)), 'b' ); hold('on');
    polar( params.phaseangle([end 1:end]), smooth((params.rmag([end 1:end])-min(params.rmag))/(max(params.rmag)-min(params.rmag)),5)', 'r' );
    polar( params.phaseangle([end 1:end]),   (mean(params.rmag)-min(params.rmag))/(max(params.rmag)-min(params.rmag))*ones(1,length(params.phaseangle)+1), 'r' );

% Instantaneous frequency
subplot(325);
    cfplot_phamp_amp_instfreq( params.phaseangle, params.rif, params.perm.rif, pvalue )
