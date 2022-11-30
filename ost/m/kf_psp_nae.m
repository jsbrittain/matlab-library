function kf_psp_nae(time,q,c95);
%function kf_psp_nae(time,q,c95);
%
% Plot Normalised Activity Envelope (NAE)
%
% Input parameters
%       time        Vector of NAE relative-time
%       q           NAE
%       c95         95% conf. interval based on assumption of independence
%
%function kf_psp_nae(time,q,c95);

plot(time,q,'k');
if (exist('c95'))
    if (~isempty(c95))
        hold('on');
        plot(time([1 end]),[0 0],'k--',time([1 end]),-c95*[1 1],'k',time([1 end]),c95*[1 1],'k');
    end;
end;
axis('tight');
title('NAE');
