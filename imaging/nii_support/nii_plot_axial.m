function nii_plot_axial(nii,threshold)
%function nii_plot_axial(nii,threshold)
%
% Supplementary plot routine for Nifti Matlab toolbox
%
%function nii_plot_axial(nii,threshold)

% Threshold image at zero
img=nii.img;
if (exist('threshold','var'))
    img(img<threshold)=threshold;
end;

% Display X-ray composition
imagesc(fliplr(squeeze(mean(img,3)))');
colormap(gray(256));
set(gca,'DataAspectRatio',nii.hdr.dime.pixdim([2 3 4]));
