function [sp,params]=kf_geom(dat1,dat2,trig,offset,duration,width,rate,seg_pwr,opt_str);
%function [sp,params]=kf_geom(dat1,dat2,trig,offset,duration,width,rate,seg_pwr,opt_str);
%
% Calculates geometric mean derived spectra, phase and coherence.  Also
% displays a summary plot with comparisons to their equivalent
% arithmetically determined quantities.
%
% Input parameters as kftf_spm2w
% Forces options 'Q0' and 'k0' (zero process noise and null smoothing)
%
% Output parameter
%       sp          Spectral matrix of geometrically determined spectra, phase
%                   and coherence.
%       params      Parameters structure.
%
%function [sp,params]=kf_geom(dat1,dat2,trig,offset,duration,width,rate,seg_pwr,opt_str);

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
%w=(width.*duration)./sum(width.*duration);
%w=w(:,ones(1,f_max))';
w=1/size(sp2,3);

% Geometric spectra
sp(:,1)=exp(sum(w.*log(faa),2));                                                    % Autospectra ch.1
sp(:,2)=exp(sum(w.*log(fbb),2));                                                    % Autospectra ch.2
sp(:,3)=real(exp(sum(w.*log(real(fab)),2)))+i*real(exp(sum(w.*log(imag(fab)),2)));  % Cross-spectra
sp(:,4)=abs(sp(:,3)).^2./(sp(:,1).*sp(:,2));                                        % Coherence
sp(:,5)=angle(sp(:,3));                                                             % Phase

% Display results if output matrices not required
if (nargout)
    
    % Calculate arithmetic derived statistics
    faa=sum(w.*faa,2);
    fbb=sum(w.*fbb,2);
    fab=sum(w.*fab,2);
    coh_arit=(abs(fab).^2)./(faa.*fbb);
    phase_arit=angle(fab);
    
    % Plot results
    figure;
    subplot(2,2,1); plot(freqs,log10(faa),freqs,log10(sp(:,1))); axis('tight'); xlim([0 freqs(end)]);
    xlabel('FREQ (Hz)'); ylabel('x10 dB'); title('Autospectra ch.1');
    subplot(2,2,2); plot(freqs,log10(fbb),freqs,log10(sp(:,2))); axis('tight'); xlim([0 freqs(end)]);
    xlabel('FREQ (Hz)'); ylabel('x10 dB'); title('Autospectra ch.2');
    legend('Arithmetic','Geometric');
    subplot(2,2,3); plot(freqs,phase_arit,freqs,sp(:,5)); xlim([0 freqs(end)]);
    xlabel('FREQ (Hz)'); ylim([-pi pi]); title('Phase');
    subplot(2,2,4); plot(freqs,coh_arit,freqs,sp(:,4),[0 freqs(end)],(1-0.05^(1/(M-1)))*[1 1],'k--');
    xlim([0 freqs(end)]); xlabel('FREQ (Hz)'); title('Coherence');
    
end;
