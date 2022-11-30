function [f,t,cl,sp] = pc_xyz(dat1,dat2,sp0,rate,seg_pwr,opt_str);
% function [f,t,cl,sp] = pc_xyz(dat1,dat2,sp0,rate,seg_pwr,opt_str);
%
% Function to calculate spectral coefficient matrix for 1st order partial spectra
% for three time series: f_xy/z.
%
% Inputs dat1    Input time series
%        dat2    Output time series
%        sp0     Predictor time series
%        rate    Sampling rate
%        seg_pwr Segment power for periodogram length
%        opt_str Options string for sp2a2_m1
%
% Outputs 
%
% function [f,t,cl,sp] = pc_xyz(dat1,dat2,sp0,rate,seg_pwr,opt_str);

[f_dummy,t_dummy,cl_1(1),sp_1(:,:,1)]=sp2a2_m1(sp0,dat1,rate,seg_pwr,opt_str);
[f_dummy,t_dummy,cl_1(2),sp_1(:,:,2)]=sp2a2_m1(sp0,dat2,rate,seg_pwr,opt_str);
[f_dummy,t_dummy,cl_1(3),sp_1(:,:,3)]=sp2a2_m1(dat1,dat2,rate,seg_pwr,opt_str);
[f,t,cl,sp]=f_partial1(sp_1,cl_1);

