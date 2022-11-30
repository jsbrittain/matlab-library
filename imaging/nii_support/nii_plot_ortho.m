function nii_plot_ortho(nii,threshold)
%function nii_plot_ortho(nii,threshold)
%
% Supplementary plot routine for Nifti Matlab toolbox
%
%function nii_plot_ortho(nii,threshold)

if (~exist('threshold','var'))
    threshold=-Inf;
end;

% Orthographic representation
figure;
subplot(2,2,1); nii_plot_coronal(nii,threshold);
subplot(2,2,2); nii_plot_sagittal(nii,threshold); 
subplot(2,2,3); nii_plot_axial(nii,threshold);
