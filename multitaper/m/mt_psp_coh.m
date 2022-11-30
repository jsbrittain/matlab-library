function h=mt_psp_coh(sp11,sp22,sp12,params,maxf,opts);
%function h=mt_psp_coh(sp11,sp22,sp12,params,[maxf],[opts]);
%
% Coherence plotting routine for use with
%   mt_sp, mt_sp2, mt_sp2_pp
%
% Plot multitaper estimates with jackknife errors if available
%
% Input parameters
%       sp11        Auto-spectrum ch.1
%       sp22        Auto-spectrum ch.2
%       sp12        Cross-spectrum ch.1,2
%       params      Parameter structure
%       maxf        (opt) Upper plot frequency
%       opts        (opt) Options string
%                       c<n>    Crop y-limit (default: 0=N/A)
%                       p<0|1>  Pointwise confidence limits (default: 0)
%                       k<col>  Change default colour
%
%function h=mt_psp_coh(sp11,sp22,sp12,params,[maxf],[opts]);

% Default parameters
topcrop=0;              % Upper y-limm (0=N/A)
pointwise=logical(0);   % Pointwise confidence limits
col='k';                % Default output colour

% Check input parameters
if (~exist('opts'))
    opts='';
end;

% Parse options string
options=deblank(opts);
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'c'                    % Crop y-limit
            topcrop=str2num(optarg);
            if (isempty(topcrop))
                error(['Error in option argument -- c' optarg]);
            end;
        case 'p'                    % Pointwise confidence limits
            pointwise=logical(1);
        case 'k'                    % Colour output
            col=optarg;
        otherwise                   % Options for wavelet analysis
            error(['Unknown option -- ' opt]);
    end;
end;

% Determine common parameters
if (params.jackknife)
    % Spectral confidence limits
    p=0.05;
    ttest=tq(1-p/2,params.jkcount-1);
end;
% Handles for conf. lines
h=[];

% Check if label specified
if (~isfield(params,'what'))
    what='';
else
    what=[': ' params.what];
end;

% Determine maximum plot frequency
if (~exist('maxf'))
    maxf=params.rate/2;
end;
fmin=find(params.freqs>0); fmin=fmin(1);
fmax=find(min(abs(params.freqs-maxf))==abs(params.freqs-maxf));
frange=fmin:fmax;
freqs=params.freqs(fmin:fmax);

% Calculate empirical 95% conf. limit
L=double(params.L(fmin:fmax));
R95=(1-0.05.^(1./(L-1)));

% Plot coherence
if (isfield(params,'coh'))
    % Jack-knife or bootstrap confidence limits
	if (params.jackknife)
        if (params.bootstrap), cohv=params.bscohv(fmin:fmax);
        else                   cohv=params.jkcohv(fmin:fmax);
        end;
        if (pointwise)
            lowerc95=params.coh(fmin:fmax)-ttest*sqrt(cohv);
            lowerc95(lowerc95<0)=0;
            %h(end+1)=plot(freqs,tanh(lowerc95).^2); hold('on');                     % Pointwise
            %h(end+1)=plot(freqs,tanh(params.coh(fmin:fmax)+ttest*sqrt(cohv)).^2);   %  |
            fill([fliplr(freqs) freqs],[flipud(tanh(lowerc95).^2); tanh(params.coh(fmin:fmax)+ttest*sqrt(cohv)).^2],0.7*[1 1 1],'linestyle','none');
            hold('on');
        else
            h(end+1)=plot(freqs,tanh(ttest*sqrt(cohv)).^2); hold('on');             % Null limit
        end;
    end;
    % Empirical confidence limits
    if ((~params.jackknife) && (pointwise))
        %lowerc95=tanh(params.coh-1.96./sqrt(2*L)).^2;
        %upperc95=tanh(params.coh+1.96./sqrt(2*L)).^2;
        plot(freqs,tanh(params.coh-1.96./sqrt(2*L)).^2,[col '--']); hold('on');	% Pointwise
        plot(freqs,tanh(params.coh+1.96./sqrt(2*L)).^2,[col '--']);             %  |
        %fill([fliplr(freqs) freqs],[flipud(tanh(lowerc95).^2); tanh(params.coh(fmin:fmax)+ttest*sqrt(cohv)).^2],0.7*[1 1 1],'linestyle','none');
    else
        plot(freqs,R95,[col '--']); hold('on');                                          % Null limit
    end;
    % Plot transformed and squared magnitude coherence
    outh=plot(freqs,tanh(params.coh(frange)).^2,col);
	ylims=ylim;
    if (ylims(2)>1)
        ylim([0 1]);
    end;
else
    % Overlay confidence limits
    coh=abs(sp12(fmin:fmax)).^2./(sp11(fmin:fmax).*sp22(fmin:fmax));
    hold('on');
    if (pointwise)          % Pointwise
        lower95=real(max(0,tanh(atanh(sqrt(coh))-1.96./sqrt(2*L)))).^2;
        upper95=real(max(0,tanh(atanh(sqrt(coh))+1.96./sqrt(2*L)))).^2;
        %plot(freqs,lower95,[col '--'],freqs,upper95,[col '--']);
        fill([fliplr(freqs) freqs],[flipud(lower95); upper95],0.7*[1 1 1],'linestyle','none');
    end;
    % Display coherence
    outh=plot(freqs,coh,col);
    
    plot(freqs,R95,[col '--']);  % Null limit
end;
if (isfield(params,'bandwidth'))
    if (isfield(params,'coh')), coh=params.coh;
    else coh=abs(sp12(fmin:fmax)).^2./(sp11(fmin:fmax).*sp22(fmin:fmax));
    end;
    if (topcrop==0), topcrop=1; end;
    xpos=freqs(1)+0.95*(freqs(end)-freqs(1));
    ypos=0.95*topcrop;
    plot(xpos+[-1 1]*params.bandwidth/2,ypos*[1 1],col); hold('on');
elseif (topcrop>0)
    ylim([0 topcrop]);
end;
box('on');
xlabel('FREQ (Hz)');
xlim([0 freqs(end)]);
title(['coh: ' what]);

% Format confidence limits
set(h,'color',[128 128 128]/255);       % Solid gray line

% Determine output parameters
if (nargout>0)
    h=outh;
else
    clear('h');
end;
