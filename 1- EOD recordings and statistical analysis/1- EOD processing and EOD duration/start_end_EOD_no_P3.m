function [p_to_p,threshold,thres_criterion,iT1,tT1,vT1,iT2,tT2,vT2,P3] = start_end_EOD_no_P3(hard_threshold,base_pts_sd,npts,base_norm_voltages,center_time,Norm_base_pts_sd)
%% calculate start and endpoints of EOD

%log that P3 is not assessed
P3 = NaN;

%extract voltage values of P1 and P2
[P1, ~] = max(base_norm_voltages);
[P2, iP2] = min(base_norm_voltages);

%peak-to-peak amplitude: (if normalized, this is 1)
p_to_p = P1 - P2;

%set a hard threshold: a % of peak-to-peak amplitude
threshold = hard_threshold*p_to_p/100;
thres_criterion = sprintf('%.1f %%', hard_threshold);

%if wave is too noisy, set threshold from std_dev of baseline
if Norm_base_pts_sd > threshold	
     threshold = 3*Norm_base_pts_sd;
     thres_criterion = '3 std_dev';
end


%Find start and end points.
%% start point: 
% the first point of the first three consecutive points that go over the threshold
for j = 1:npts
    if abs(base_norm_voltages(j)) > threshold && abs(base_norm_voltages(j+1)) > threshold && abs(base_norm_voltages(j+2)) > threshold
        iT1 = j;
        tT1 = center_time(j);
        vT1 = base_norm_voltages(j);
        break
    end
end

%% end point:

% I want the end point of the EOD to be when the waveform levels off after
% climbing from P2 (i.e., when it looks "horizontal")

% start at iP2, find where the negative threshold is met (this will ignore
% P3s)
k = iP2;
while base_norm_voltages(k) < -threshold 
    k = k + 1;
end

% I will consider the voltage returned to baseline when the average
% voltage of the previous 11 points (including k) is larger than the
% negative threshold (this will treat P3s as normal returns to baseline)

potT2 = mean(base_norm_voltages(k-10:k));
while potT2 <= -threshold
    k = k + 1;
    potT2 = mean(base_norm_voltages(k-10:k));
end

iT2 = k;
tT2 = center_time(k);
vT2 = base_norm_voltages(k);

end
