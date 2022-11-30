function nii_plot_ortho_mip(nii,threshold)
%function nii_plot_ortho_mip(nii,threshold)
%
% Supplementary plot routine for Nifti Matlab toolbox
%
%function nii_plot_ortho_mip(nii,threshold)

if (~exist('threshold','var'))
    threshold=-Inf;
end;

% Orthographic representation
figure;
subplot(2,2,1); nii_plot_coronal_mip(nii,threshold);
subplot(2,2,2); nii_plot_sagittal_mip(nii,threshold); 
subplot(2,2,3); nii_plot_axial_mip(nii,threshold);
