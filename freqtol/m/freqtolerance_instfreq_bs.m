function analysis = freqtolerance_instfreq_bs( instfreq, instmag, precision, bscount, bssamples, deltaIF, varargin )
%function analysis = freqtolerance_instfreq_bs( instfreq, instmag, precision, bscount, bssamples, [deltaIF] )
%
% Bootstrap version
%
%function analysis = freqtolerance_instfreq_bs( instfreq, instmag, precision, bscount, bssamples, [deltaIF] )

% Parameter defaults
if (~exist('deltaIF','var')), deltaIF = []; end;
if (isempty(deltaIF))
    deltaIF = instfreq(2:end) - instfreq(1:end-1);
end;
if (~exist('bssamples','var')), bssamples = []; end;
if (isempty(bssamples))
    bssamples = length(deltaIF);
end;
if (~exist('bscount','var')), bscount = []; end;
if (isempty(bscount))
    bscount = 1;
end;
if (bscount==0), bscount = 1; end;

% Bootstrap procedure
fprintf('Bootstrap sampling (%g iters, %g/%g samples [%g%%]) \n<',bscount,bssamples,length(deltaIF),round(100*bssamples/length(deltaIF)));
progress = 0;
for k = (1:bscount)
    currentprogress = round( 100*k/bscount );
    if (currentprogress > progress )
        progress = currentprogress;
        if ((mod(progress,10)==0) && (progress~=100))
            fprintf('%g',progress);
        else
            fprintf('.');
        end;
    end;
    if ( bscount == 1 )
        kk = (1:length(deltaIF));
    else
        kk = randi(length(deltaIF),bssamples,1);
    end;
    analysis(k) = freqtolerance_instfreq( instfreq(kk), instmag(kk), precision, deltaIF(kk), varargin{:} );
end;
fprintf('>\n');
