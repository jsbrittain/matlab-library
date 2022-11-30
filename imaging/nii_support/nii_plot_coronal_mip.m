function nii_plot_coronal_mip(nii,threshold)
%function nii_plot_coronal(nii,threshold)
%
% Supplementary plot routine for Nifti Matlab toolbox
%
%function nii_plot_coronal(nii,threshold)

% Threshold image at zero
img=nii.img;
if (exist('threshold','var'))
    img(img<threshold)=threshold;
end;

% Display X-ray composition
imagesc(fliplr(squeeze(max(img,[],2)))');
colormap(gray(256));
set(gca,'DataAspectRatio',nii.hdr.dime.pixdim([4 2 3]));
