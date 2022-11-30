function kf_psp_spm2w_grid(sp,params,opt_str,chlabels)
%function kf_psp_spm2w_grid(sp,params,opt_str,chlabels)
%
% Grid display spectra, coherence and phase (by default)
%
% Takes redundant cell array input, ie.
%   sp{1,1}, ... , sp{1,M},
%    sp{2,2}, ... , sp{2,M},
%     sp{3,3}, ... , sp{3,M},
%          |        |       |
%       sp{M-1,M-1}, sp{M-1,M},
%       sp{M,M}
%
% Input parameters
%   sp          Cell array of spectra    { M x M }
%   params      Cell array of parameters { M x M }
%   opt_str     (opt) Options string
%                 s<n,m>    Smooth images over <n> Hz and <m> trials
%                 f<n>      Plot frequencies to <n> Hz
%                 l         Linear spectra
%                 c         Upper limit on coherence plot
%   chlabels    (opt) Cell array of channel labels
%
%function kf_psp_spm2w_grid(sp,params,opt_str,chlabels)

% Check input data
if ((~iscell(sp)) || (~iscell(params)))
    error(' Input must be cell arrays.');
end
if (sum([diff(size(sp)) diff(size(params))])~=0)
    error(' Input must be square cell arrays.');
end

% Default parameters
smooth_n=1; smooth_m=1;
logspectra=true;
maxf=params{1,1}.rate/2;
maxcoh=1;

% Parse options string
options=deblank(opt_str); opt_str='';
while (any(options))                    % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));          % Determine option argument.
	switch (opt(1))
	    case 's'                        % Smooth output
            smooth_n = str2num(optarg);
            smooth_m=smooth_n(2);
            smooth_n=smooth_n(1);
        case 'l'                        % Linear-spectra
            logspectra=false;
        case 'f'                        % Maximum plot frequency
            maxf=str2num(optarg);
        case 'c'                        % Maximum coherence limit
            maxcoh=str2num(optarg);
    end
end

% Convert smooth parameter from Hz to (frequency) samples
smooth_n=round(smooth_n/diff(params{1,1}.freqs([1 2])));
if (smooth_n<1)
    smooth_n=1;     % Null smoother
end
smoothing_kernel = ones(smooth_n,smooth_m)/(smooth_n*smooth_m);

% Determine data parameters
M=length(sp);

% Check channel labels
if (~exist('chlabels','var'))
    for ind=(1:M)
        chlabels{ind}=num2str(ind);
    end
end

% Plot grid-array
figure;
for ch1=(1:M)
    for ch2=(1:M)
        if (ch1==ch2)       % Diagonal (spectra)
            
            subplot(M,M,(ch1-1)*M+ch2);
            %kf_psp_a(sp{ch1,ch2},params{ch1,ch2},'');
            
            flen=size(sp,1);
            faa=10*log10(squeeze(double(sp{ch1,ch2}(:,1,:))));
            imagesc([1 size(faa,2)],params{ch1,ch2}.freqs,convolve2(faa,smoothing_kernel,'same'));
            
            ylabel(''); xlabel(''); title([num2str(ch1) ': ' chlabels{ch1}]);
            
        elseif (ch1<ch2)    % Upper quadrant (Coherence)
            
            subplot(M,M,(ch1-1)*M+ch2);
            %kf_psp_coh(sp{ch1,ch2},params{ch1,ch2},'');
            coh=squeeze(double(sp{ch1,ch2}(:,4,:)));
            imagesc([1 size(coh,2)],params{ch1,ch2}.freqs([1 end]),convolve2(coh,smoothing_kernel,'same'));
            set(gca,'clim',[0 maxcoh]);
            ylabel(''); xlabel('');
            
        else                % Lower quandrant
            
            %mt_psp_ph(smooth(sp12{ch2,ch1},smoother),params{ch2,ch1},maxf,smooth(sp11{ch2,ch1},smoother),smooth(sp22{ch2,ch1},smoother));
            %title(''); xlabel('');
            %if (ch1==1)
            %    ylabel(num2str(ch2));
            %end;
            
        end;
    end;
end;
