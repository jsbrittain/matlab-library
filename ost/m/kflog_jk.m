function [jk,m]=kflog_jk(sp,w);
%function [jk,m]=kflog_jk(sp,w);
%
% Construct jackknife estimates for kflog analysis
%
% Input parameters
%       sp      Matrix of periodograms (ln domain)
%       w       Weight vector
%
% Output parameters
%       jk      Matrix of (jackknifed) variances
%
%function [jk,m]=kflog_jk(sp,w);

% Check input parameters
if (size(sp,2)~=4)
    error(' Incorrect spectral matrix dimensions');
end;
if ((size(w,1)~=1) & (size(w,2)~=1))
    error(' Spectral weights must be in vector form');
end;

% Determine data parameters
fcount=size(sp,1);
M=size(sp,3);

% Rotate weight vector
if (size(w,2)>size(w,1))
    w=transpose(w);
end;

% Construct delete-one estimates
d=zeros(fcount,4,M);
wm=permute(w(:,ones(1,fcount),ones(1,1,4)),[2 3 1]);
wmsp=sum(wm.*sp,3);
if (M==1)               % Single trial (No variance)
    d(:,:,1)=wmsp;
else                    % Multiple trials
    for ind=1:M
        currentsp=w(ind)*sp(:,:,ind);
        d(:,:,ind)=(wmsp-currentsp)/(sum(w)-w(ind));                    % Auto- and (real/imag) cross- spectra
	end;
end;
d(:,3,:)=log(real(exp(d(:,3,:)))+i*real(exp(d(:,4,:))));           % Complex cross-spectra
d(:,4,:)=atanh(abs( exp(d(:,3,:)) )./sqrt(exp(d(:,1,:)).*exp(d(:,2,:)))); % Mag. coherency (with normalising transform)
d(:,5,:)=angle(exp(d(:,3,:)));                                          % Phase

% Calculate mean parameters
wm=wm(:,1,:);
m(:,1:2)=wmsp(:,1:2,:);
%m(:,3)=log(real(exp(sum(wm.*sp(:,3,:),3))))+i*log(real(exp(sum(wm.*sp(:,4,:),3))));
m(:,3)=log(real(exp(sum(wm.*sp(:,3,:),3)))+i*real(exp(sum(wm.*sp(:,4,:),3))));
%m(:,4)=atanh(abs( real(exp(real(m(:,3))))+i*real(exp(imag(m(:,3)))) )./sqrt(exp(m(:,1)).*exp(m(:,2))));
m(:,4)=atanh(abs( exp(m(:,3)) )./sqrt(exp(m(:,1)).*exp(m(:,2))));
m(:,5)=angle(exp(m(:,3)));

% Estimate variance
jk=(M-1)*mean(abs(d-m(:,:,ones(1,1,M))).^2,3);

% Ensure no NaN values occur (Generated by atanh(1)=Inf)
nanvalues=find(isnan(jk(:,4)));
jk(nanvalues,4)=atanh(1-eps);