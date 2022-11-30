function [f,t,cl,sp] = sp2_p_fn(d1,d2,samp_rate,seg_samp_min,seg_size,samp_tot,mains_flag,d3);
% function [f,t,cl,sp] = sp2_p_fn(d1,d2,samp_rate,seg_samp_min,seg_size,samp_tot,mains_flag,d3);
%
%

if (nargin==8)      % Check for a partial coherence predictor
    partial=1;
else
    partial=0;
end;

seg_tot=size(d1,2); % Calculate total number of segments
fd1=fft(d1);        % DFT across columns/segments ch 1, PBMB (4.1) or (4.2)
fd2=fft(d2);        % DFT across columns/segments ch 2, PBMB (4.1) or (4.2)

t_fac=2*pi*samp_tot;% Normalisation for periodogram spectral estimates.  samp_tot is equivalent to R=LT

f11=sum((abs(fd1.*fd1)/t_fac)')';   % Spectrum 1, PBMB (5.2), NB Mag squared for autospectra
f22=sum((abs(fd2.*fd2)/t_fac)')';   % Spectrum 2, PBMB (5.2), NB Mag squared for autospectra
f21=sum((fd2.*conj(fd1)/t_fac)')';  % Cross spectrum (complex valued), PBMB (5.2)

if (partial)
    % Calculate cross and auto spectra for predictor d3
    fd3=fft(d3);
    f13=sum((fd1.*conj(fd3)/t_fac)')';  % Cross spectrum of d1 and predictor
    f23=sum((fd2.*conj(fd3)/t_fac)')';  % Cross spectrum of d2 and predictor
    f33=sum((abs(fd3.*fd3)/t_fac)')';   % Predictor autospectrum
    
    % Find any zeros in f33 and add eps (minimum possible value) to them to
    % avoid 'divide by zero' errors
    f33=f33+(f33==0)*eps;
    
    % Partial autospectra
    f11=f11 - abs(f13.*f13)./f33;       % Partial autospectra (i.e. excluding effect from d3)
    f22=f22 - abs(f23.*f23)./f33;
    
    % Partial cross spectra
    f21=f21 - (f23.*conj(f13))./f33;    % Partial cross spectra (i.e. excluding effect from d3)
end;

deltaf=samp_rate/seg_size;  % Reolution - Spacing of Fourier frequencies in Hz

if (mains_flag)     % Suppression of mains - smooth out using adjacent values
    mains_ind=round(50.0/deltaf)+1;     % NB index 1 is DC
    f11(mains_ind)=0.5*(f11(mains_ind-2)+f11(mains_ind+2)); % Spectrum ch 1.
    f11(mains_ind-1)=0.5*(f11(mains_ind-2)+f11(mains_ind-3));
    f11(mains_ind+1)=0.5*(f11(mains_ind+2)+f11(mains_ind+3));
    f22(mains_ind)=0.5*(f22(mains_ind-2)+f22(mains_ind+2)); % Spectrum ch 2.
    f22(mains_ind-1)=0.5*(f22(mains_ind-2)+f22(mains_ind-3));
    f22(mains_ind+1)=0.5*(f22(mains_ind+2)+f22(mains_ind+3));
    f21(mains_ind)=0.5*(f21(mains_ind-2)+f21(mains_ind+2)); % Cross spectrum
    f21(mains_ind-1)=0.5*(f21(mains_ind-2)+f21(mains_ind-3));
    f21(mains_ind+1)=0.5*(f21(mains_ind+2)+f21(mains_ind+3));
    % Smooth elements in upper hermitian section of cross spectral
    % estimate.  This data used in ifft() to generate cumulant. NB Data is
    % complex conjugate.
    f21(seg_size-mains_ind+2)=conj(f21(main_ind));
    f21(seg_size-mains_ind+3)=conj(f21(main_ind-1));
    f21(seg_size-mains_ind+1)=conj(f21(main_ind+1));
end;

% Construct output spectral matrix f.
seg_size_2=(2:seg_size/2+1)';   % Indexing for output, DC component not output
f(:,1)=(seg_size_2-1)*deltaf;   % Column 1 - frequencies in Hz
f(:,2)=log10(f11(seg_size_2));  % Column 2 - Log spectrum ch 1.
f(:,3)=log10(f22(seg_size_2));  % Column 3 - Log spectrum ch 2.
f(:,4)=abs(f21(seg_size_2)).*abs(f21(seg_size_2))./(f11(seg_size_2).*f22(seg_size_2));
f(:,5)=angle(f21(seg_size_2));  % Column 5 - Phase, PBMB (5.7)

% Estimate cumulant density using inverse DFT of cross spectrum
deltat=1000.0/samp_rate;    % dt in msec

cov=ifft(f21);              % Inverse DFT

% Construct output time domain matrix t.
% Column 1 - time in msec. Range (-T/2)*dt to (T/2-1)*dt
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;

% Column 2 - Cumulant, shifted by T/2 so that time zero is centre.
% NB 2pi/T factor is 2*pi since ifft includes 1/T term
t([seg_size/2+1:seg_size,1:seg_size/2],2)=real(cov(1:seg_size))*2*pi;   % PBMB (5.9)

% Estimate variance of cumulant density estimate.
var_fac=4*pi*pi/(seg_size*samp_tot);        % Factor (2pi/T) (2pi/R)
q_var=var_fac*2*sum(f11(1:seg_size/2).*f22(1:seg_size/2));              % PBMB (6.10)

% Construct sp structure
sp.f11=f11(seg_size_2);
sp.f22=f22(seg_size_2);
sp.f21r=real(f21(seg_size_2));
sp.f21i=imag(f21(seg_size_2));

% Construct cp structure, confidence limits for parameter estimates
cl.seg_size=seg_size;       % Set values in cl structure
cl.seg_tot=seg_tot;
cl.samp_tot=samp_tot;
cl.samp_rate=samp_rate;
cl.dt=deltat;
cl.df=deltaf;
cl.q_c95=1.96*sqrt(q_var);  % Confidence limits for cumulant, PBMB (6.11)
