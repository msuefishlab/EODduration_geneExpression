function [measurement_data] = standard_eod_measurement(wave,s_rate,samp_name,period,treatment,treatdate,individual)
%%%%%%%%%%
%STANDARD_EOD_MEASUREMENT
%performs EOD measurements on previously normalized EOD waveform, based on
%measurescript.m
%modified 3 November 2008, for JRG's study of signal variation in P.
%kingsleyae, based on MEA's code for analyzing variation in P. magnostipes
%reproductive isolation.

%%--Settings: Should be Paramaterized...---
firstforty = wave(1:40);
stddev = std(firstforty);
hard_threshold = 0.5;

p_to_p=max(wave)-min(wave);

threshold = hard_threshold*p_to_p/100;

npts=length(wave);
n_base_pts = 512;
base_pts = wave([1:n_base_pts npts-n_base_pts+1:npts]);
baseline_std = std(base_pts);


if baseline_std > threshold				%if wave too noisy
        threshold = 3*Norm_base_pts_sd;
end
[~,n_pts]=size(wave);
measurement_data.period=period;
measurement_data.treatement=treatment;
measurement_data.treatdate=treatdate;
measurement_data.individual=individual;
measurement_data.npts=n_pts;
time=linspace(0,1000*n_pts/s_rate,n_pts);

%% --calc first derivative used for calculating slopes---------
ddtwave=ddteod(wave,s_rate);									%returns vector 'ddtwave'

%% --Describe Peak P1, always present, by definition ------ %%
[~,measurement_data.iP1]=max(wave);
time=time-time(measurement_data.iP1);
measurement_data.tP1 = time(measurement_data.iP1);                                          %time of P1
measurement_data.vP1 = wave(measurement_data.iP1);                                          %voltage at P1
[measurement_data.sS1,measurement_data.iS1] = max(ddtwave(1:measurement_data.iP1));         %index (iS1) and Slope (sS1) of Maximum Slope of P1 (S1)
measurement_data.tS1 = time(measurement_data.iS1);                                          %Time of S1
measurement_data.vS1 = wave(measurement_data.iS1);                                          %voltage at S1

%% --Describe Peak P2, assumed always present ------ %%
[measurement_data.vP2, measurement_data.iP2] = min(wave(measurement_data.iP1:n_pts)); % Voltage and Index of P2
measurement_data.iP2 = measurement_data.iP1 + measurement_data.iP2 - 1;               % Adjust index of P2
measurement_data.tP2 = time(measurement_data.iP2);                                   % Time of P2

[measurement_data.sS2, measurement_data.iS2] = min(ddtwave(measurement_data.iP1:measurement_data.iP2)); % Index (iS2) and Slope (sS2) of Maximum Slope of P1-P2 (S2)
measurement_data.iS2 = measurement_data.iP1 + measurement_data.iS2 - 1;              % Adjust index of S2
measurement_data.tS2 = time(measurement_data.iS2);                                   % Time of S2
measurement_data.vS2 = wave(measurement_data.iS2);                                   % Voltage of S2

%% --Describe Beginning and End of Waveform, always present ---- %%
[measurement_data.vT1,measurement_data.tT1,measurement_data.iT1,measurement_data.vT2,measurement_data.tT2,measurement_data.iT2]=eod_start_end(wave,time,threshold,measurement_data.iP2);

%% --Describe P0.  Always Assumed to be present -------- %%
for k = measurement_data.iP1:-1:1				
	if wave(k) < 0	       
        measurement_data.iZC1 = k;				
		break;
    else
        measurement_data.iZC1 = measurement_data.iP1;
    end
end

measurement_data.tZC1 = time(measurement_data.iZC1);				%time of ZC1
measurement_data.vZC1 = wave(measurement_data.iZC1);				%voltage at ZC1
measurement_data.sZC1 = ddtwave(measurement_data.iZC1);             %slope at ZC1 

width=fix(s_rate*.0005);
markpoint=measurement_data.iZC1-width;
measurement_data.aP0=sum(wave(markpoint:measurement_data.iZC1))/s_rate;
measurement_data.aP0=measurement_data.aP0/1e-6;

[measurement_data.vP0,measurement_data.iP0] = min(wave(markpoint:measurement_data.iZC1)); %attempt to find the minimum
measurement_data.tP0 = time(measurement_data.iP0+markpoint);     %attempt to find the time of minimum


%% --Fit Decay of P2 to y = a * exp(-τ * t) + b ------ %%
%% --Fit Decay of P2 to y = a * exp(-τ * t) + b ------ %%
% Define the region of interest (decay phase after P2, up to iT2)
decay_indices = measurement_data.iP2:measurement_data.iT2;
time_decay = time(decay_indices) - measurement_data.tP2; % Time relative to P2
voltage_decay = wave(decay_indices);

% Use MATLAB's curve fitting toolbox to fit the decay
ft = fittype('a * exp(-tau * t) + b', 'independent', 't', 'coefficients', {'a', 'tau', 'b'});

% Initial guesses for parameters: [a, τ, b]
initial_a = measurement_data.vP2 - mean(voltage_decay(end-5:end)); % Peak voltage minus baseline
initial_tau = 0.001; % Small initial guess for the time constant
initial_b = mean(voltage_decay(end-5:end)); % Baseline voltage

% Fit the decay model using robust fitting
fit_opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                      'StartPoint', [initial_a, initial_tau, initial_b], ...
                      'Lower', [-Inf, 0, -Inf], ...
                      'Upper', [Inf, Inf, Inf], ...
                      'Robust', 'LAR'); % LAR handles outliers better

% Ensure time_decay and voltage_decay are column vectors
time_decay = time_decay(:);
voltage_decay = voltage_decay(:);

[fit_result, gof] = fit(time_decay, voltage_decay, ft, fit_opts);

% Extract optimized parameters
a = fit_result.a;
tau = fit_result.tau;
b = fit_result.b;

% Store results
measurement_data.decay_a = a;
measurement_data.decay_tau = tau;
measurement_data.decay_b = b;
measurement_data.time_decay=time_decay;
measurement_data.voltage_decay=voltage_decay;
measurement_data.fit_result=fit_result;
measurement_data.gof=gof;


%% --Zero Crossings ------ %%
dumpos1 = n_pts-1;						%default near end of wave
dumpos2 = n_pts;						%default at end of wave
for i = measurement_data.iP1:+1:n_pts					%start at P1 and read forward
	if wave(i) <= threshold		
   	dumpos1 = i - 1;				%index just before wave crosses +threshold
		break;     						%break out of for loop and leave dumpos1 = corrected i
  	end
end
for i = dumpos1:+1:n_pts				%start at dumpos1 and read forward
	if wave(i) <= (-threshold)	
     	dumpos2 = i;					%index where wave crosses -threshold
     	dumpos2_present = 1;			%if so, set cond for phase2 presence to true      
     	break;    						%break out of for loop and leave dumpos2 = new i
    end
end

if dumpos2_present					%if condition true, find ZC2 landmarks
   j = round((dumpos1 + dumpos2)/2);
   if abs(wave(dumpos1))<abs(wave(j))
      measurement_data.iZC2=dumpos1;
   elseif abs(wave(dumpos2))<abs(wave(j))
     	measurement_data.iZC2=dumpos2;					%these statements check for rounding error
   else									%and find the actual voltage nearest zero-xing
      measurement_data.iZC2=j;							%in this (usually) steep voltage transition
   end
  	measurement_data.tZC2 = time(measurement_data.iZC2);				%time of ZC2
  	measurement_data.vZC2 = wave(measurement_data.iZC2);				%voltage at ZC2
  	measurement_data.sZC2 = ddtwave(measurement_data.iZC2);			%slope at ZC2
end

%% --calculate phase areas---------
measurement_data.aP1 = (sum(wave(measurement_data.iT1:measurement_data.iZC2)))/s_rate;
measurement_data.aP2 = (sum(wave(measurement_data.iZC2:measurement_data.iT2)))/s_rate;

%% --calculate the power spectrum-------
PS=power_spectrum_eod(wave,s_rate);	%call to func that plots power spectr & calcs fftstats

measurement_data.Fmax = PS.Fmax;
measurement_data.Pmax = PS.Pmax;
measurement_data.Flow = PS.Flow;
measurement_data.Plow = PS.Plow;
measurement_data.Fhigh = PS.Fhigh;
measurement_data.Phigh = PS.Phigh;
measurement_data.Bndwidth = PS.Fhigh - PS.Flow;
measurement_data.sample_name=string(samp_name);

