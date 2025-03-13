function [base_pts_mean,base_pts_sd,noNorm_p_to_p,Norm_base_pts_sd,base_voltages,base_norm_voltages,iP1,iP2,center_time] = baseline_normalize_center(npts,voltages,time)
%adjust baseline and center the raw voltages recorded

%baseline substraction: subtract mean amplitude of the first and last 512 points (1024 points total, half of the total points per recording (npts = 2048)). 
n_base_pts = 512;
base_pts = voltages([1:n_base_pts npts-n_base_pts+1:npts]);

base_pts_mean = mean(base_pts);
base_pts_sd = std(base_pts);

base_voltages = voltages - base_pts_mean;

%normalize peak-to-peak amplitude:

%extract voltage values and indexes of P1 and P2
[P1, iP1] = max(base_voltages);
[P2, iP2] = min(base_voltages);

%peak-to-peak amplitude
noNorm_p_to_p = P1 - P2;

%normalization
base_norm_voltages = base_voltages/noNorm_p_to_p;

% the standard dev of the baseline should also be normalized (this is
% equivalent to calculating the sd dev on the normalized voltages)
Norm_base_pts_sd = base_pts_sd/noNorm_p_to_p;

%Center. Set time = 0 at P1
tP1 = time(iP1);
center_time = time - tP1;

end