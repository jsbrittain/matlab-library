function [Wsp,wlparam]=wlsp2_gabor_epochs(dat1,dat2,rate,df,frange,opt_str,mother,wlopt);
%function [Wsp,wlparam]=wlsp2_gabor_epochs(dat1,dat2,rate,df,frange,opt_str,mother,wlopt);
%
% Wavelet spectrum analysis (Type 2)
%
% Calculate auto-spectra, cross-spectra, coherence and phase for
% two channels.
%
% Recursive implementation
%
% Input parameters
%   Data parameters
%       dat1        Time series 1
%       dat2        Time series 2
%       trig        Trigger
%       rate        Sampling rate
%   Analysis parameters
%       df          Frequency resolution
%       frange      Frequency range [fmin fmax] ([]=default)
%       opt_str     Options string
%                       r<0,1,2>    rectify (none, ch1, ch1&2)
%                       p           pad with zeros
%   Wavelet parameters
%       mother      Mother wavelet
%       wlopt       Wavelet options
%
% Ouput parameters
%       Wsp         Spectral coefficients (3rd dimension)
%                       1 - Auto-spectra ch.1
%                       2 - Auto-spectra ch.2
%                       3 - Cross-spectra
%                       4 - Coherence
%                       5 - Phase
%       wlparam     Wavelet parameters
%
%function [Wsp,wlparam]=wlsp2_gabor_epochs(dat1,dat2,rate,df,frange,opt_str,mother,wlopt);

% Remove singleton dimensions
dat1=squeeze(dat1);
dat2=squeeze(dat2);

% Extract epoch parameters
offset = 0;                         % Offset (null)
dur=size(dat1,1);                   % Duration of epoch (samples)
duration = dur/rate*1000;           %  | (msecs)
M = size(dat1,2);                   % Trial count

% Flatten time-series
dat1=dat1(:);
dat2=dat2(:);

% Determine trigger locations
trig = (1:dur:M*dur);

% Perform sp2 analysis
[Wsp,wlparam]=wlsp2_gabor(dat1,dat2,trig,offset,duration,rate,df,frange,opt_str,mother,wlopt);
