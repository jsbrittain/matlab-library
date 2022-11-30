function [jk,m]=kfroot_jk(sp,w,root);
%function [jk,m]=kfroot_jk(sp,w,root);
%
% Construct jackknife estimates for kfroot analysis
%
% Input parameters
%       sp      Matrix of periodograms (root-domain)
%       w       Weight vector
%       root    Root-transform
%
% Output parameters
%       jk      Matrix of (jackknifed) variances
%
%function [jk,m]=kfroot_jk(sp,w,root);

% Check input parameters
if (nargin~=3)
    error(' Incorrect number of input parameters');
end;
if (size(sp,2)~=4)
    error(' Incorrect spectral matrix dimensions');
end;

% Determine data parameters
fcount=size(sp,1);
M=size(sp,3);

% Construct delete-one estimates
d=zeros(fcount,4,M);
wmsp=sum(double(w).*sp,3);
if (M==1)               % Single trial (No variance)
    d(:,:,1)=wmsp;
else                    % Multiple trials
    for ind=1:M
        currentsp=double(w(:,:,ind)).*sp(:,:,ind);
        d(:,:,ind)=(wmsp-currentsp)./(sum(double(w),3)-double(w(:,:,ind)));     % Auto- and (real/imag) cross- spectra
	end;
end;
d(:,3,:)=(real((d(:,3,:)).^root)+i*real((d(:,4,:)).^root)).^(1/root);           % Complex cross-spectra
d(:,4,:)=atanh(abs(d(:,3,:).^root)./sqrt((d(:,1,:).^root).*(d(:,2,:).^root)));  % Mag. coherency (with normalising transform)
d(:,5,:)=angle(d(:,3,:).^root);                                                 % Phase

% Calculate mean parameters
m(:,1:2)=wmsp(:,1:2,:);
m(:,3)=(real((sum(double(w(:,3,:)).*sp(:,3,:),3)).^root)+i*real((sum(double(w(:,4,:)).*sp(:,4,:),3)).^root)).^(1/root);
m(:,4)=atanh(abs(m(:,3).^root)./sqrt((m(:,1).^root).*(m(:,2).^root)));
m(:,5)=angle(m(:,3).^root);

% Estimate variance
jk=(M-1)*mean(abs(d-m(:,:,ones(1,1,M))).^2,3);

% Ensure no NaN values occur (Generated by atanh(1)=Inf)
nanvalues=find(isnan(jk(:,4)));
jk(nanvalues,4)=atanh(1-eps);
