function [coh1,coh2,dz,dzv]=jk_comp_coh(params1,params2);
%
% Function to compare coherence
%
% Applicable to the multi-sample case of unequal sample sizes.
% Estimated quantity: magnitude coherency (sqrt coherence)
%
% References
%   H.Bokil, D.J.Thomson & P.P.Mitra (2005), Comparing Spectra and Coherences
%       for Samples of Unequal Size.  Conference paper (689.22).
%

% Coherence comparison from jackknife estimates

% Determine DOF (2m,2n) for data (m,n = no. jackknife estimates)
m=size(params1.jk11,2);
n=size(params2.jk11,2);

% Calculate transformed coherence minus bias term (all samples) - from pseudovalues!
%       --- CHECK SPECTRA - USING MEAN(JK11) NOT SP11 OR MEAN(PV11); is closer! ---
coh1=atanh(abs(mean(params1.jk12,2))./sqrt(mean(params1.jk11,2).*mean(params1.jk22,2)))-1/(2*m-2);
coh2=atanh(abs(mean(params2.jk12,2))./sqrt(mean(params2.jk11,2).*mean(params2.jk22,2)))-1/(2*n-2);

% Calculate transformed coherence minus bias term (delete-one estimates)
z1=atanh(abs(params1.jk12)./sqrt(params1.jk11.*params1.jk22))-1/(2*m-2);
z2=atanh(abs(params2.jk12)./sqrt(params2.jk11.*params2.jk22))-1/(2*n-2);

% Determine dz (scaled diff of transformed coherence)
dz=(coh1-coh2)/sqrt(1/(2*m-2)+1/(2*n-2));
dz1=(z1-coh2(:,ones(1,n)))/sqrt(1/(2*(m-1)-2)+1/(2*n-2)); % delete-one dz (ch.1)
dz1=m*dz(:,ones(1,m))-(m-1)*dz1;                          %  => pseudovalues
dz2=(coh1(:,ones(1,m))-z2)/sqrt(1/(2*m-2)+1/(2*(n-1)-2)); % delete-one dy (ch.2)
dz2=n*dz(:,ones(1,n))-(n-1)*dz2;                          %  => pseudovalues

% Calculate variance
meandz1=mean(dz1,2); meandz1=meandz1(:,ones(1,m));
meandz2=mean(dz2,2); meandz2=meandz2(:,ones(1,n));
dzv=mean((dz1-meandz1).^2,2)/(m-1)+mean((dz2-meandz2).^2,2)/(n-1);

%disp([num2str(length(find((abs(dz))>1.96))/length(coh1)*100) '% outside 95%']);
