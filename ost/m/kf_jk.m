function jk=kf_jk(sp,w,jk);
%function jk=kf_jk(sp,w,jk);
%
% Construct jackknife estimates for kf analysis
%
% Input parameters
%       sp      Matrix of periodograms (time domain)
%       w       Weight matrix
%
% Output parameters
%       jk      Matrix of (jackknifed) variances
%
%function jk=kf_jk(sp,w,jk);

warning(' Jackknife estimates for time-domain kf analysis are not currently available');
jk=zeros(size(sp,1),5);

return;

% Check input parameters
if (size(sp,2)~=3)
    error(' Incorrect spectral matrix dimensions');
end;

% Determine data parameters
fcount=size(sp,1);
M=size(sp,3);

% Rotate weight vector
w=reshape(transpose(w),fcount,3,M);

% Construct delete-one estimates
d=zeros(fcount,3,M);
wsp=sum(w.*sp,3);
if (M==1)               % Single trial (No variance)
    d(:,:,1)=wsp;
else                    % Multiple trials
    for ind=1:M
        currentsp=w(:,:,ind).*sp(:,:,ind);
        d(:,:,ind)=(wsp-currentsp)./(sum(w,3)-w(:,:,ind));      % Auto- and cross- spectra
	end;
end;
d(:,4,:)=atanh(abs(d(:,3,:))./sqrt(d(:,1,:).*d(:,2,:)));        % Mag. coherency (with normalising transform)
d(:,5,:)=angle(d(:,3,:));                                       % Phase

% Calculate mean parameters
m(:,1:3)=wsp;
m(:,4)=atanh(abs(m(:,3))./sqrt(m(:,1).*m(:,2)));
m(:,5)=angle(m(:,3));

% Estimate variance
jk=(M-1)*mean(abs(d-m(:,:,ones(1,1,M))).^2,3);
