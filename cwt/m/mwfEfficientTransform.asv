function [S11,S22,S12,wlparam]=mwfEfficientTransform(dat1,dat2,dt,df,frange,opt,mother,wlopt);
%function [S11,S22,S12,wlparam]=mwfEfficientTransform(dat1,dat2,dt,df,frange,opt,mother,wlopt);
%
% Efficient Muliwavelet Transform
%
% Performs the usual multi-wavelet transform one scale at a time to reduce
% memory requirements.  The final auto- and cross-spectra are returned as
% 'single' precision matrices, also to reduce memory requirements.
%
% Input parameters
%   dat1        Time-series ch.1
%   dat2        Time-series ch.2
%   dt          Time resolution
%   df          Frequency resolution (1/(scales/octave))
%   frange      Frequency range [fmin fmax]
%   opt         Analysis options (as mwfTransform), (log-scale only)
%   mother      Mother wavelet
%   wlopt       Wavelet options
%
% Weighting schemes (internal variable)
%   0 - Uniform
%   1 - Eigenspectra (default)
%
% Supported wavelet
%   'morse'     {beta, gamma, A(area)}
%
%function [S11,S22,S12,wlparam]=mwfEfficientTransform(dat1,dat2,dt,df,frange,opt,mother,wlopt);

% Determine min and max scale
N=length(dat1);
f_min=frange(1);
f_max=frange(2);

% Calculate scales
wlparam=wlp_morse(wlopt{:});
s_max=1./(wlparam.fourier_factor*f_min);
s_min=1./(wlparam.fourier_factor*f_max);

% Convert min frequency to max scale
s0=s_min;
dj=df;
J=ceil((1/dj)*log2(s_max/s0));
wlparam.scale=s0.*2.^([0:J]*dj);
scales=wlparam.scale;

% Reserve variable space
sizeW=[length(wlparam.scale) N];
S11=single(zeros(sizeW));
S22=single(zeros(sizeW));
S12=single(zeros(sizeW));

% Perform multi-wavelet transform one scale at a time
for ind=1:length(scales)
    disp(['Scale ' int2str(ind) ' of ' int2str(length(scales))]);
    [S11a,S22a,S12a,wlparam]=mwTransform(dat1,dat2,dt,dj,scales(ind),opt,mother,wlopt);
    S11(ind,:)=single(S11a);
    S22(ind,:)=single(S22a);
    S12(ind,:)=single(S12a);
end;
wlparam.wlopt=wlopt;
wlparam.k=wlopt{end};
wlparam.scale=scales;
wlparam.linearscale=0;
wlparam.freqs=1./(wlparam.fourier_factor*scales);
