function analysis = freqtolerance_instfreq( instfreq, instmag, precision, deltaIF, varargin )
%function analysis = freqtolerance_instfreq( instfreq, instmag, precision, [deltaIF], varargin )
%
% Parameters
%       instfreq    [Nx1] vector of instantaneous frequency (zero crossing) entries
%       instmag     [Nx1] vector of instantaneous amplitude entries
%       precision   Frequency resolution of analysis (Hz), e.g. 0.1
%       deltaIF     (opt) [Nx1] Change in frequency associated with each instfreq entry
%                   This permits fast downsample and bootstrap estimates
%
% Parameter-value pairs:
%       criterion   Selection criterion, options include:
%                       'X2':   Chi-squared error (weighted sum of squares)
%                       'X2r':  Reduced chi-squared
%                       'AIC':  AIC (uncorrected)
%                       'AICc': (default) AIC small sample corrected
%                   Whichever criterion is selection is substituted into
%                   the aic_1 and aic_3 fields and used for model selection
%       dfrange     Allowable range of the (mean binned) deltaIF vector
%
%function analysis = freqtolerance_instfreq( instfreq, instmag, precision, [deltaIF], varargin )

% Input parse
p = inputParser;
addRequired(p,'instfreq',@isnumeric);
addRequired(p,'instmag',@isnumeric);
addRequired(p,'precision',@isnumeric);
addOptional(p,'deltaIF',[],@isnumeric);
addParamValue(p,'criterion',[]);
addParamValue(p,'dfrange',[],@isnumeric);
addParamValue(p,'instsnr',[],@isnumeric);
parse(p,instfreq,instmag,precision,deltaIF,varargin{:});
% Assign parsed values
criterion = p.Results.criterion;
dfrange = p.Results.dfrange;
instsnr = p.Results.instsnr;
if (isempty(criterion))
    criterion = 'AICc';
end

% Display some options for the user
%fprintf('Using %s selection criterion\n',criterion);
%fprintf('      +/- %g delta-instfreq range\n',dfrange);

% Default parameters on optional arguments
if (isempty(deltaIF))
    deltaIF = diff( instfreq );
end
if (isempty(instsnr))
    instsnr = NaN*instmag;
end
% Normalise vector lengths
instfreq = instfreq(1:length(deltaIF));
instmag = instmag(1:length(deltaIF));
instsnr = instsnr(1:length(deltaIF));

% Stratify frequency
ifstrat = round(instfreq/precision);
for fn = (1:((10)/precision))        % Max frequency (Hz)
    ifchangeMed(fn)   = nanmedian( deltaIF( ifstrat == fn ) );
    ifchangeMean(fn)  =   nanmean( deltaIF( ifstrat == fn ) );
    ifchangeStd(fn)   =    nanstd( deltaIF( ifstrat == fn ) );
    ifchangeCount(fn) =    nansum(        ( ifstrat == fn ) );
    ampMean(fn)       =   nanmean( instmag( ifstrat == fn ) );
    ampStd(fn)        =    nanstd( instmag( ifstrat == fn ) );
    ampSNR(fn)        =   nanmean( instsnr( ifstrat == fn ) );
    freqs(fn) = fn*precision;
end

if ( ~isempty(dfrange) )
    keep = ( ifchangeMean > -dfrange ) & ( ifchangeMean < dfrange );
    ifchangeMed( ~keep ) = NaN;
    ifchangeMean( ~keep ) = NaN;
    ifchangeStd( ~keep ) = NaN;
    ifchangeCount( ~keep ) = NaN;
    ampMean( ~keep ) = NaN;
    ampStd( ~keep ) = NaN;
    ampSNR( ~keep ) = NaN;
end

% Quantify -- std dev crossing points
fc = nanmean( instfreq );
fn = dsearchn( freqs.', fc );
f1 = freqs( find( (ifchangeMean(1:fn)  -0.5*ifchangeStd(1:fn)) > 0, 1, 'last' ) );
f2 = freqs( find( (ifchangeMean(fn:end)+0.5*ifchangeStd(fn:end)) < 0, 1, 'first' ) + fn - 1 );

% Optimisation constraints
min_gradient = 0.25;
max_gradient = 10.0;
separation_freq = fc;
min_freq = separation_freq - 4*nanstd(instfreq);
max_freq = separation_freq + 4*nanstd(instfreq);

% Quantify -- piecewise fit
x0 = [ 1 1 fc+nanstd(instfreq)*[-1 1] ];
options = optimoptions(@fmincon,'Algorithm','active-set','Display','off');
[piecefit,fval_3] = fmincon(@(x)trilinearstep_err(x,freqs,ifchangeMean,ifchangeStd),x0,[0 0 1 -1; -1 0 0 0; 0 -1 0 0],[0; 0; 0],[],[],[min_gradient min_gradient min_freq separation_freq],[max_gradient max_gradient separation_freq max_freq],[],options);
bestfit_3 = trilinearstep(piecefit,freqs);
Np = length(x0);                                                % Number of parameters
Nb = sum( (~isnan(ifchangeMean)) & (~isnan(ifchangeStd)) );     % Number of bins
switch ( criterion )
    case {'X2'}         % Chi-squared error (weighted sum of squares)
        aic_3 = fval_3;
        aic_0 = linearstep_err([0 0],freqs,ifchangeMean,ifchangeStd);
    case {'X2r'}        % Reduced chi-squared
        aic_3 = fval_3/(Nb-Np-1);                                       % Reduced chi-square
        aic_0 = linearstep_err([0 0],freqs,ifchangeMean,ifchangeStd)/(Nb-1);
    case {'AIC'}        % AIC (uncorrected)
        loglik = nansum( -0.5*(log(2*pi)+log(ifchangeStd.^2)+((bestfit_3-ifchangeMean).^2)./(ifchangeStd.^2)) );
        aic_3 = -2*loglik + 2*Np;
        loglik0 = nansum( -0.5*(log(2*pi)+log(ifchangeStd.^2)+((0*bestfit_3-ifchangeMean).^2)./(ifchangeStd.^2)) );
        aic_0 = -2*loglik0;
    case {'AICc'}       % AIC (small sample corrected)
        loglik = nansum( -0.5*(log(2*pi)+log(ifchangeStd.^2)+((bestfit_3-ifchangeMean).^2)./(ifchangeStd.^2)) );
        chi2r = -2*loglik/(Nb-Np-1);                                    % Reduced chi-square
        aic_3 = chi2r + 2*(Np*(Np+1)/(Nb*Np-1));                        % Small-sample corrected AIC
        loglik0 = nansum( -0.5*(log(2*pi)+log(ifchangeStd.^2)+((0*bestfit_3-ifchangeMean).^2)./(ifchangeStd.^2)) );
        aic_0 = -2*loglik0/(Nb-1);
    otherwise, error(' Unknown criterion selected: %s',criterion);
end

% Quantify -- linear fit
linearfittype = 2;
switch ( linearfittype )
    case 1  % Original method (polyfit does not account for std.dev in fitting the data)
        fval_1 = (( ifchangeMean(~isnan(ifchangeMean)) - polyval(polyfit(freqs(~isnan(ifchangeMean)),ifchangeMean(~isnan(ifchangeMean)),1),freqs(~isnan(ifchangeMean))) ).^2)./(ifchangeStd(~isnan(ifchangeMean)).^2);
        fval_1 = nansum( fval_1(~isinf(fval_1)) );
        aic_1 = fval_1 + 2*2;
        bestfit_1 = polyval(polyfit(freqs(~isnan(ifchangeMean)),ifchangeMean(~isnan(ifchangeMean)),1),freqs);
    case 2  % Chi-square error criteria
        x0 = [ 1 fc ];
        switch ( 2 )
            case 1      % Unconstrained fit
                options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');
                [linearfit,fval_1] = fminunc(@(x)linearstep_err(x,freqs,ifchangeMean,ifchangeStd),x0,options);
            case 2      % Same constraints as trilinear fit
                [linearfit,fval_1] = fmincon(@(x)linearstep_err(x,freqs,ifchangeMean,ifchangeStd),x0,[],[],[],[],[min_gradient min_freq],[max_gradient max_freq],[],options);
                bestfit_1 = linearstep(linearfit,freqs);
        end
        % AIC calculation
        Np = length(x0);                                                % Number of parameters
        Nb = sum( (~isnan(ifchangeMean)) & (~isnan(ifchangeStd)) );     % Number of bins
        switch ( criterion )
            case {'X2'}         % Chi-squared error (weighted sum of squares)
                aic_1 = fval_1;
            case {'X2r'}        % Reduced chi-squared
                aic_1 = fval_1/(Nb-Np-1);                                       % Reduced chi-square
            case {'AIC'}        % AIC (uncorrected)
                loglik = nansum( -0.5*(log(2*pi)+log(ifchangeStd.^2)+((bestfit_1-ifchangeMean).^2)./(ifchangeStd.^2)) );
                aic_1 = -2*loglik + 2*Np;
            case {'AICc'}       % AIC (small sample corrected)
                loglik = nansum( -0.5*(log(2*pi)+log(ifchangeStd.^2)+((bestfit_1-ifchangeMean).^2)./(ifchangeStd.^2)) );
                chi2r = -2*loglik/(Nb-Np-1);                                    % Reduced chi-square
                aic_1 = chi2r + 2*(Np*(Np+1)/(Nb*Np-1));                        % Small-sample corrected AIC
            otherwise, error(' Unknown criterion selected: %s',criterion);
        end
end

% Probability of broad-tolerance (under AIC)
rellik1 = exp((aic_0-aic_1)/2);
rellik3 = exp((aic_0-aic_3)/2);

% Best fit lines
if ( aic_3 < aic_1 )
    % Piecewise linear
    bestfit = bestfit_3;
else
    % Linear
    bestfit = bestfit_1;
end

if ( false )
    % Quantify -- direct measure of tolerance width by std.dev about mean -- longest chain of contiguous values
    span = 5;
    lowerlimit = smooth( ifchangeMed-ifchangeStd./sqrt(ifchangeCount), span ); lowerlimit( isnan(ifchangeMed) ) = NaN;
    upperlimit = smooth( ifchangeMed+ifchangeStd./sqrt(ifchangeCount), span ); upperlimit( isnan(ifchangeMed) ) = NaN;
    v = find( (lowerlimit<=0) & (upperlimit>=0) );
    x = [0; cumsum(diff(v)~=1)];     % Find continugous indices
    if (isempty(v))
        directwidth = 0;
    else
        v = v(x==mode(x));          % Indicies of most common (longest) chain
        directwidth = diff(freqs(v([1 end])));
    end
end

% Record analysis
analysis.instfreq           = instfreq(:);
analysis.instmag            = instmag(:);
analysis.instsnr            = instsnr(:);
analysis.deltaIF            = deltaIF(:);

analysis.freqs              = freqs;
analysis.ifchangeMed        = ifchangeMed;
analysis.ifchangeMean       = ifchangeMean;
analysis.ifchangeStd        = ifchangeStd;
analysis.ifchangeCount      = ifchangeCount;
analysis.ampMean            = ampMean;
analysis.ampStd             = ampStd;
analysis.ampSNR             = ampSNR;

analysis.fit.piecefit  = piecefit;
analysis.fit.freqtol   = (aic_1>aic_3)*(piecefit(4)-piecefit(3));
analysis.fit.aic_1     = aic_1;
analysis.fit.fval_1    = fval_1;
analysis.fit.linearfit = linearfit;
analysis.fit.linearfittype = linearfittype;
analysis.fit.aic_3     = aic_3;
analysis.fit.fval_3    = fval_3;
analysis.fit.aic_0     = aic_0;
analysis.fit.aic_rel_lik_1 = rellik1;
analysis.fit.aic_rel_lik_3 = rellik3;
analysis.fit.criterion = criterion;

analysis.fit.bestfit   = bestfit;
analysis.fit.bestfit_1 = bestfit_1;
analysis.fit.bestfit_3 = bestfit_3;

if (exist('directwidth','var'))
    analysis.directwidth = directwidth;
end
