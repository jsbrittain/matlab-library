function nii=nii_fuse_ctmr(ctfilename,mrfilename,ctmrfilename,reslice,th)
%function [nii]=nii_fuse_ctmr(ctfilename,mrfilename,ctmrfilename,[reslice,[th]])
%
% Assumes CT and MR are in the same space
%
% Requires third-party Nifti & Analyze Matlab toolbox
%
% Input parameters
%   ctfilename
%
% Output parameters
%
% 
%
%function [nii]=nii_fuse_ctmr(ctfilename,mrfilename,ctmrfilename,[reslice,[th]])

% Check input arguments
if ((nargin<3) || (nargin>5))
    error(' Incorrect number of input parameters.');
end;
if (~exist('reslice'))
    reslice=false;
end;
if (~exist('th'))
    th=0.001;       % 0.1%
end;

% Reslice images if requested
if (reslice)
    disp('Reslicing MR image');
    reslice_nii(mrfilename,[mrfilename(1:end-4) '_resliced.nii']);
    mrfilename=[mrfilename(1:end-4) '_resliced.nii'];
    disp('Reslicing CT image');
    reslice_nii(ctfilename,[ctfilename(1:end-4) '_resliced.nii']);
    ctfilename=[ctfilename(1:end-4) '_resliced.nii'];
else
    warning(' Not reslicing has flipped the CT in the past!');
end;

% Load MR and CT files
% niimr=load_nii(mrfilename);
% niict=load_nii(ctfilename);
niimr=load_untouch_nii(mrfilename);
niict=load_untouch_nii(ctfilename);

% Remove intensities less than zero
niimr.img(niimr.img<0)=0;
niict.img(niict.img<0)=0;
% Reject intensities larger than `th' of total integral
[y,x]=hist(single(niimr.img(:)),1000);
maxintmr=x(find((1-cumsum(y)/sum(y))<th,1,'first'));
niimr.img(niimr.img>maxintmr)=maxintmr;
[y,x]=hist(single(niict.img(:)),1000);
maxintct=x(find((1-cumsum(y)/sum(y))<th,1,'first'));
niict.img(niict.img>maxintct)=maxintct;
% Normalise intensity range (0-1)
niimr.img=single(niimr.img)/single(maxintmr);
niict.img=single(niict.img)/single(maxintct);

% Create new (composite) image (scale image to enable 3D rendering in FSLview)
nii=make_nii((niimr.img+niict.img)*1024);%,niimr.hdr.dime.pixdim(2:4));

% Save image file
save_nii(nii,ctmrfilename);

% Check output arguments
if (nargout<1)
    clear('nii');
end;
