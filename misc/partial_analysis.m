function [fp1_out, tp1_out, clp1_out, sp1_out] = partial_analysis(experiment, ch1, ch2, ch3, draw_graphs)
% Script to perform partial spectra analysis on
% experiment EMG1(ch2) EMG2(ch3) with EEG(ch1) as predictor ( f23/1 )
%
% Input:
%   experiment
%   ch1, ch2    Input and output channels
%   ch3         Predictor channel
%   draw_graphs (Optional)  If set (=1) then graphs drawn one for each
%                           segment
%
% Output:
%   [fp1, tp1, clp1, sp1]
%   correspond to partial spectra analysis over n segments
%   therefore fp1(:,:,1) is the partial analysis of segment 1

% debug (allows scripting)
%experiment = Exp405k1010;
%ch1=2;
%ch2=3;
%ch3=1;
% end debug

if nargin < 5
    draw_graphs = 0;
end;

offset = (-500:200:-100)';
tot_seg=3;
duration=200;
seg_pwr_t=12;
t_fac = experiment.rate_t / experiment.rate;

[f_1,t_1,cl_1,sp_1] = tf_sp2a21(experiment.d(:,ch3), experiment.d(:,ch2), experiment.t/t_fac, experiment.rate, offset, duration, seg_pwr_t, 't2 r2 a0 a21', experiment.name);
[f_2,t_2,cl_2,sp_2] = tf_sp2a21(experiment.d(:,ch3), experiment.d(:,ch1), experiment.t/t_fac, experiment.rate, offset, duration, seg_pwr_t, 't2 r2 a0 a21', experiment.name);
[f_3,t_3,cl_3,sp_3] = tf_sp2a21(experiment.d(:,ch1), experiment.d(:,ch2), experiment.t/t_fac, experiment.rate, offset, duration, seg_pwr_t, 't2 r2 a21 a21', experiment.name);

% Partial spectral analysis on segments
for segment=1:tot_seg
    sp(:,:,1)=sp_1(:,:,segment);
    sp(:,:,2)=sp_2(:,:,segment);
    sp(:,:,3)=sp_3(:,:,segment);

    cl(1)=cl_1(segment);
    cl(2)=cl_2(segment);
    cl(3)=cl_3(segment);

    [fp1, tp1, clp1, sp1] = f_partial1(sp, cl);
    
    if draw_graphs == 1
        figure
        psp2(fp1, tp1, clp1, 60, 500, 400, 1);
    end;
    
    fp1_out(:,:,segment) = fp1;
    tp1_out(:,:,segment) = tp1;
    clp1_out(segment) = clp1;
    sp1_out(:,:,segment) = sp1;
end;

