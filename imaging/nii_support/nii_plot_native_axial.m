function nii_plot_native_axial(nii,threshold)
%function nii_plot_native_axial(nii,threshold)
%
% Supplementary plot routine for Nifti Matlab toolbox
%
% Scale to native space (MNI coordinate system)
%
%function nii_plot_native_axial(nii,threshold)

% Threshold image at zero
img=nii.img;
if (exist('threshold'))
    img(img<threshold)=threshold;
end;

% Display X-ray composition
imagesc(fliplr(squeeze(mean(img,3)))');
colormap(gray(256));
set(gca,'DataAspectRatio',nii.hdr.dime.pixdim([2 3 4]));
