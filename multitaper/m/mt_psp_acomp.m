function mt_psp_acomp(sp,params,maxf,logplot)
%
%
% Takes an cell-array of `sp' and `params' as input and overlays with
% errorbars.
%

% Check input parameters
if (~exist('maxf'))
    maxf=[];
end;
if (~exist('logplot'))
    logplot=true;
end;

% Default parameters
if (isempty(maxf))
    maxf=0;
    for ind=(1:length(params))
        maxf=max([maxf max(params{ind}.freqs)]);
    end;
end;

% Plot parameters
col=[100 100 255; 255 100 100; 100 255 100; 255 255 100]/255;
col2='brgy';
confscaling=1.96;
alphalevel=0.6;

% Prepare plot
hold('on');

% Recurse datasets, plot confidence limits
for ind=(length(sp):-1:1)
    switch (logplot)
        case true,
            conf_lower=10*log10(exp(log(sp{ind})-confscaling*sqrt(1./params{ind}.L)));
            conf_upper=10*log10(exp(log(sp{ind})+confscaling*sqrt(1./params{ind}.L)));
        case false,
            conf_lower=exp(log(sp{ind})-confscaling*sqrt(1./params{ind}.L));
            conf_upper=exp(log(sp{ind})+confscaling*sqrt(1./params{ind}.L));
    end;
    
    % Confidence intervals
    h=fill([params{ind}.freqs fliplr(params{ind}.freqs)],[conf_lower; flipud(conf_upper)]',col(ind,:));
    set(h,'linestyle','none'); alpha(h,alphalevel);
end;

% Recurse datasets, plot mean spectra
for ind=(length(sp):-1:1)
    switch (logplot)
        case true,
            mdat=10*log10(sp{ind});
        case false,
            mdat=sp{ind};
    end;    
    % Overlay mean ERP
    h=plot(params{ind}.freqs,mdat,col2(ind));
    set(h,'linewidth',2);
end;

% Close plot editing
xlim([0 maxf]);
hold('off');
