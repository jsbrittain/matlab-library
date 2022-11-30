function layered_coherence(experiments, ch1, ch2)

% Draw layered coherence diagram
%
% See fig.15 of Walk1 locomotion report (panum)

%
% Usage:
%   layered_coherence( [Exp405k1010 Exp406k1004 Exp406k1006], [2 2 2], [3 3 3] );
%

freq = 60;
ch_max = 1;
offset=(-500:200:-100)';
duration = 200;
seg_pwr_t = 10;
[dumb, num_exp] = size(experiments);

figure
for n=1:num_exp
    [f0,t0,cl0] = tf_sp2a2(experiments(:,n).d(:,ch1(n)), experiments(:,n).d(:,ch2(n)), ...
                experiments(:,n).t*(experiments(:,n).rate/experiments(:,n).rate_t), ...
                experiments(:,n).rate, offset, duration, seg_pwr_t, 't2 r2 a21 a21', ...
                experiments(:,n).name);
    for k=1:3
        subplot(2,2,k);
        hold on;
        psp_ch1(f0(:,:,k), cl0(:,k), freq, ch_max);
    end;
end;
    
hold off;