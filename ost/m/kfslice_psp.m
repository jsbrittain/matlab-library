function kfslice_psp(sp,params,freq);
%function kfslice_psp(sp,params,freq);
%
% Kalman-Fourier frequency slice summary plotting routine
% for auto-spectra, coherence and phase
%
% Input parameters
%       sp          Spectral matrix
%       params      Parameters structure
%       freq        Slice frequency
%
%function kfslice_psp(sp,params,freq);

% Determine label
label=kf_label(params);

% Plot spectral summary of adaptive method
figure;
subplot(2,2,1); kfslice_psp_a(sp,params,freq,label);
subplot(2,2,2); kfslice_psp_b(sp,params,freq,label);
subplot(2,2,3); kfslice_psp_ph(sp,params,freq,label);
subplot(2,2,4); kfslice_psp_coh(sp,params,freq,label);
