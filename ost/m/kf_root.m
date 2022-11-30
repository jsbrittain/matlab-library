function [sp,params]=kf_root(dat1,dat2,trig,offset,duration,width,rate,seg_pwr,opt_str);
%function [sp,params]=kf_root(dat1,dat2,trig,offset,duration,width,rate,seg_pwr,opt_str);
%
% Calculates root-transform derived spectra, phase and coherence.  Also
% displays a summary plot with comparisons to their equivalent
% arithmetically determined quantities.
%
% Input parameters as kftf_spm2w
% Forces k0 (null smoothing)
%
% Additional input parameters
%       opt_str     Additional parameters
%                       v<n>        Take n-th root variance stabilising transform
%                                   Default: 4
%
% Output parameter
%       sp          Spectral matrix of geometrically determined spectra, phase
%                   and coherence.
%       params      Parameters structure.
%
%function [sp,params]=kf_geom(dat1,dat2,trig,offset,duration,width,rate,seg_pwr,opt_str);

% Default parameters
root=4;                             % Variance stabilising transform root

% Parse options string
options=deblank(opt_str);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'v'                    % Variance stabilising transform
            root=str2num(optarg);
            if isempty(root)
                
            end;
    	otherwise                   % Pass on to analysis
            opt_str=[opt_str ' ' opt];
    end;
end;

% Calculate periodograms (use KF routines with null smoothing)
opt_str=[opt_str ' k0'];
[sp2,params]=kftf_spm2w2(dat1,dat2,trig,0,1,offset,duration,width,rate,seg_pwr,opt_str);
freqs=params.freqs;
f_max=size(sp2,1);
M=size(sp2,3);

% Extract spectra
faa=squeeze(sp2(:,1,:));
fbb=squeeze(sp2(:,2,:));
fab=squeeze(sp2(:,3,:));

% Determine weighting scheme
w=(width.*duration)./sum(width.*duration);
w=w(:,ones(1,f_max))';

% Root spectra
sp(:,1)=(sum(w.*(faa.^(1/root)),2)).^(root);                                        % Autospectra ch.1
sp(:,2)=(sum(w.*(fbb.^(1/root)),2)).^(root);                                        % Autospectra ch.2
sp(:,3)=real((sum(w.*(real(fab).^(1/root)),2)).^(root)) + ...                       % Cross-spectra
        i*real((sum(w.*(imag(fab).^(1/root)),2)).^(root));                          %  |
sp(:,4)=abs(sp(:,3)).^2./(sp(:,1).*sp(:,2));                                        % Coherence
sp(:,5)=angle(sp(:,3));                                                             % Phase

% Calculate arithmetic and geometric mean
coh_arit=(abs(sum(w.*fab,2)).^2)./(sum(w.*faa,2).*sum(w.*fbb,2));

% Display results
figure;
subplot(2,2,1); plot(freqs,log10(sum(w.*faa,2)),freqs,log10(sp(:,1))); axis('tight'); xlim([0 freqs(end)]);
xlabel('FREQ (Hz)'); ylabel('x10 dB'); title('Autospectra ch.1');
subplot(2,2,2); plot(freqs,log10(sum(w.*fbb,2)),freqs,log10(sp(:,2))); axis('tight'); xlim([0 freqs(end)]);
xlabel('FREQ (Hz)'); ylabel('x10 dB'); title('Autospectra ch.2');
legend('Arithmetic',[int2str(root) '-root']);
subplot(2,2,3); plot(freqs,angle(sum(w.*fab,2)),freqs,sp(:,5)); xlim([0 freqs(end)]);
xlabel('FREQ (Hz)'); ylim([-pi pi]); title('Phase');
subplot(2,2,4); plot(freqs,coh_arit,freqs,sp(:,4),[0 freqs(end)],(1-0.05^(1/(M-1)))*[1 1],'k--');
xlim([0 freqs(end)]); xlabel('FREQ (Hz)'); title('Coherence');
