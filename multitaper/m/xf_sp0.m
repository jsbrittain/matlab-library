function [xf11,xf22,xf12,params] = xf_sp0(dat1,dat2,duration,rate,opt_str)
%function [xf11,xf22,xf12,params] = xf_sp0(dat1,dat2,duration,rate,opt_str)
%
% Cross-frequency (xf) correlation analysis
%
%function [xf11,xf22,xf12,params] = xf_sp0(dat1,dat2,duration,rate,opt_str)

%warning(' I have no idea if the cross-frequency output is valid!');

opt_str = [ opt_str ' e'];      % Return eigenspectra (!!! IMPLEMENTED PER EPOCH NOT PER TAPER IN MT_SP !!!)
[~,~,~,params] = mt_sp0(dat1,dat2,duration,rate,opt_str);

disp(' Computing correlations...');
xf11 = zeros(length(params.freqs));
xf22 = xf11; xf12 = xf11;
for f1 = (1:length(params.freqs))
    disp([ num2str(f1) ' of ' num2str(length(params.freqs)) ]);
    for f2 = (1:length(params.freqs))
        cc = corrcoef( params.p11(f1,:), params.p11(f2,:) ); xf11(f1,f2) = cc(1,2);
        cc = corrcoef( params.p22(f1,:), params.p22(f2,:) ); xf22(f1,f2) = cc(1,2);
        cc = corrcoef( params.p11(f1,:), params.p22(f2,:) ); xf12(f1,f2) = cc(1,2);
    end;
end;
disp('done.');
