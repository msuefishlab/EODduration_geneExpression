%% This code is intended to analyze EOD recordings from testosterone experiment
% It will correct for gain and normalize voltages, although some
% non-normalized data is collected

%% Set hard threshold for baseline: a percentage of the peak-to-peak amplitude
% %(write down the desired percentage value, without the percentage
%sign, and without dividing by 100. E.g. if you want 2%, set har
%%
% _ITALIC TEXT_ d_threshold
%= 2)
hard_threshold = 0.5;

%% set a Gain value. This is the gain set on the amplifier. 
% Note that this will affect all recordings. If different recordings have
% different gain values they must be corrected manually
gain = 20;

%% get an array with all the folders with recordings, sorted by date

% Get a list of all files and folders in the working directory
files = dir;
% Get a logical vector that tells which are directories
dirFlags = [files.isdir];
% Extract only those that are directories
folders = files(dirFlags);
% get rid of the '.' and '..'
folders = folders(~ismember({folders.name},{'.','..'}));

% sort by date. It looks like the structure array must be converted to a
% table first.
folders = struct2table(folders);
sortedFolders = sortrows(folders, 'datenum');

% extract an array with only the sorted directory names
dates = table2cell(sortedFolders(:,1));

dates_total = length(dates);

%% Make a list of all samples
%Each of these will have a final structure array with the results of its recordings

%create an empty string array
AllSamples = string.empty;

%Look for samples (= .mat files) in all recording dates
for i = 1:dates_total
    dummy = dir(fullfile(dates{i},'*.mat'));
    dummy = struct2table(dummy);
    dummy = dummy.name;
    %remove mat files from a timed recording (as opposed to single mode)
    iTimed = find(contains(dummy, 'leod.mat'));
    dummy(iTimed) = [];
    for j = 1:length(dummy)
        AllSamples(j,i) = dummy{j};
    end
end

%keep only the unique samples, without the .mat extension
AllSamples_u = unique(AllSamples);
AllSamples_u = rmmissing(AllSamples_u); % this function removes entire rows or 
%columns. Use it after the samples are converted from matrix to vector 
%form (previous step).
AllSamples_mod = extractBefore(AllSamples_u, ".mat");
%the following is necessary because variable names can't start with a
%number
AllSamples_mod = strcat("f", AllSamples_mod);

%% Create a table that will hold the results of every sample in structure arrays

Nsamples = length(AllSamples_mod);

AllResults = table('Size',[Nsamples, 1],'VariableTypes',{'cell'}, ...
    'RowNames', AllSamples_mod, 'VariableNames', {'Data'});

% Inside each sample's structure arrray, create necessary tables and
% structures. Notice that recording dates may differ between samples (each
% recording date is a row)
ColumnNames = {'fish', 'date', 'duration_mean', 'duration_std',...
    'N_recordings', 'recordings', 'fish_orientation', 'noisy_recordings',...
    'P3', 'time_averaged_EOD', 'voltage_averaged_EOD', 'joint_recordings',...
    'averagedEOD__duration', 'norm_vP0_mean', 'norm_vP0_sd', ... 
    'norm_vP1_mean', 'norm_vP1_sd', 'norm_vP2_mean', 'norm_vP2_sd', ...
    'ratio_vP0_vP1_mean', 'ratio_vP0_vP1_sd', 'ratio_vP2_vP1_mean', 'ratio_vP2_vP1_sd', ...
    'No_norm_vP0_mean', 'No_norm_vP0_sd', 'No_norm_vP1_mean', 'No_norm_vP1_sd',...
    'No_norm_vP2_mean', 'No_norm_vP2_sd', 'No_norm_Vrange_mean',... 
    'No_norm_Vrange_sd', 'PeakFreq_mean_Hz', 'PeakFreq_sd_Hz'};

VarTypes = {'string', 'string', 'double', 'double',...
    'double', 'cell', 'string', 'string',...
    'string', 'cell', 'cell', 'cell',...
    'double', 'double','double', ... 
    'double', 'double', 'double', 'double',...
    'double', 'double', 'double', 'double',...
    'double', 'double', 'double', 'double',...
    'double', 'double', 'double',... 
    'double', 'double', 'double'};
Ncols = length(ColumnNames);

for i = 1:Nsamples
    
    sample = AllSamples_mod{i};
      
    %count on how many dates the sample was recorded
    dates_rec = sum(count(AllSamples(:), AllSamples_u{i}));
    
    %create necessary tables and structures
    %AllResults{sample, 'Data'} = {struct()};    
    AllResults{sample, 'Data'}{1,1} = table('Size',[dates_rec,...
        Ncols], 'VariableTypes', VarTypes, 'VariableNames', ColumnNames);
    
    %make an empty struct inside the cell array under recordings. This is
    %necessary to please Matlab (when adding rows to that structure)
    AllResults{sample, 'Data'}{1,1}{:, 'recordings'} = {struct()};
    
    %Matlab adds zeros to cells in numeric variables. I rather not have
    %zeros assigned beforehand, The only option I found was to change them
    %to NaN:
    AllResults{sample, 'Data'}{1,1}{:, vartype('double')} = NaN;   
end

    
% Create an additional table to save additional info about each recording
recordings_extra_details = table('Size', [0, 8], 'VariableNames', {'Fish',...
    'date', 'recording', 'basepts_mean', 'basepts_std', 'iT1', 'iT2',...
    'eod_pts'}, 'VariableTypes', {'string', 'string', 'int64', 'double',...
    'double', 'int64', 'int64', 'int64'});


%% loop through the folder structure and analyze each recording

%for each date of recording, list the recorded samples. This list is a
%structure array, the column .name has the samples
for i = 1:dates_total
     
    %save the date of recording as a variable to use in the analysis
    date = convertCharsToStrings(dates{i});    
    
    % list the samples recorded on date i
    samples = dir(fullfile(dates{i},'*.mat'));
    
    %% select each sample (j) for the chosen date (i)
    for j = 1:length(samples)
        
        %extract the name of the sample
        sample = strcat("f", extractBefore(samples(j).name, ".mat"));
        
        %% The following piece is necessary to adjust for analyses where not 
        % all samples were recorded on every recording date. It will
        % determine which row in each sample's table should receive the
        % data for this date. I call this row "ii".
        
        % I need to count how many times the fish (= "samples(j).name") was
        % recorded in the dates prior to and including "date". For example,
        % if the fish was recorded for the fourth time in the date "date",
        % then the data from "date" should be stored in the fourth row of
        % that fish's table (again, I call this row "ii").
        
        % Subset which samples were recorded on each date prior to and
        % including "date". This would be the columns 1 : column of
        % recording in the string array AllSamples. Column of recording is
        % the index of "date" in the cell array "dates", which is i
        dates_subset = AllSamples(:,1:i);
        
        % Find and count how many times the sample of interest is present
        % in that subset (Note: this code assumes that each fish is never
        % recorded more than once per recording date)
        dummy0 = count(dates_subset, samples(j).name);
        ii = sum(dummy0(:)); 
        
        %% row ii will will hold the data for this sample and date            
        % log sample j and date i
        AllResults{sample, 'Data'}{1,1}{ii, 'fish'} = sample;
        AllResults{sample, 'Data'}{1,1}{ii, 'date'} = date;      
       
        %the following line loads the structed array named as 'eod' by the
        %recording software. It contains all the recordings for sample j on
        %date i
        load(fullfile(date, samples(j).name));
                  
        %number of recordings in sample j
        [~, Nrecordings] = size(eod);
        
        %% select and analyze each recording (k) from sample (j) on date (i)
        for k = 1:Nrecordings
            
            % extract sampling rate
            fs = eod(k).Rate;
            
            %select the voltage values for recording k
            voltages = eod(k).wave/gain;
            
            %check that the waveform is head positive. If it isn't,
            %multiply by -1 and print a warning 
            [~, iP1_raw] = max(voltages);
            [~, iP2_raw] = min(voltages);
            
            if iP1_raw > iP2_raw
                warning('This waveform is head negative. I multiplied all votages by (-1). Date = %s, fish = %s, recording = %d', date, sample, k)
                warn_flag = 1;
                voltages = -1*voltages;
            else
                warn_flag = 0;
            end
            
            %extract the number of points present in recording k
            [npts, ~] = size(voltages);

            %make a vector with the times in ms at which each voltage point was
            %recorded
            time = linspace(0,1000*npts/fs, npts)'; 

            %% For recording k, adjust baseline, normalize by peak-to-peak amplitude and center time at t = 0 at P1
            [base_pts_mean,base_pts_sd,noNorm_p_to_p,Norm_base_pts_sd,base_voltages,base_norm_voltages,iP1,iP2,center_time] = baseline_normalize_center(npts,voltages,time);

            %% calculate start and endpoints of the EOD present in recording
            % use this code to include P3 in analysis
          %  [p_to_p,threshold,thres_criterion,iT1,tT1,vT1,iT2,tT2,vT2,P3] = start_end_EOD_P3(hard_threshold,base_pts_sd,npts,base_norm_voltages,center_time);
            
            % use this code to exclude P3 from analysis
            [p_to_p,threshold,thres_criterion,iT1,tT1,vT1,iT2,tT2,vT2,P3] = start_end_EOD_no_P3(hard_threshold,base_pts_sd,npts,base_norm_voltages,center_time,Norm_base_pts_sd);
            
            %% Calculate EOD duration, save
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).threshold = threshold;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).threshold_crit = thres_criterion;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).duration = tT2 - tT1;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).voltages_all = base_norm_voltages;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).time_all = center_time;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).voltages_eod = base_norm_voltages(iT1:iT2);
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).time_eod = center_time(iT1:iT2);

            %% log fish orientation
            if warn_flag == 1
                AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).fish_orientation = 'head negative. I multiplied all votages by (-1)';
            else
                AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).fish_orientation = 'head positive';
            end
            
            %% log presence/absence of P3
            if P3 == 0
                AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).P3 = 'absent';
            elseif P3 == 1
                AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).P3 = 'present';
            elseif isnan(P3)
                AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).P3 = 'not_assessed';
            end
            
            %% Calculate normalized & non-normalized voltage peaks, ratios P0/P1 & P2/P1, & freq. of max power
            % Note that this code assumes P0 exists. If it doesn't exist,
            % this code will still designate a P0.
            
            % Determine normalized P0, P1, P2
            vP2 = min(AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).voltages_eod);
            [vP1, iP1] = max(AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).voltages_eod);
            vP0 = min(AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).voltages_eod(1:iP1));
            
            % Determine non-normalized values of P0, P1, P2. 
            vNoNorm_P2 = min(base_voltages);
            [vNoNorm_P1, iNoNormP1] = max(base_voltages);
            vNoNorm_P0 = min(base_voltages(1:iNoNormP1));
            
            % Determine ratios P0/P1, P2/P1, voltage range, and Log all peak values
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).norm_vP0 = vP0;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).norm_vP1 = vP1;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).norm_vP2 = vP2;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).ratio_vP0_vP1 = abs(vP0/vP1);
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).ratio_vP2_vP1 = abs(vP2/vP1);
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).No_norm_vP0 = vNoNorm_P0;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).No_norm_vP1 = vNoNorm_P1;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).No_norm_vP2 = vNoNorm_P2;
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).No_norm_Vrange = abs(vNoNorm_P2) + abs(vNoNorm_P1);     
            
            % Calculate and log the frequency of maximum power
            [PeakFreq] = power_spectrum(base_norm_voltages, fs);
            AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}(k).PeakFreq_Hz = PeakFreq;
            
            %% save additional info about each recording
            eod_pts = iT2 - iT1 + 1;
            ptsdata = {sample, date, k, base_pts_mean, base_pts_sd, iT1, iT2, eod_pts};
            recordings_extra_details = [recordings_extra_details; ptsdata];
            
        end

        
        %% calculate mean and std of duration for the (k) recordings of sample (j) on date (i)
        AllResults{sample, 'Data'}{1,1}{ii, 'duration_mean'} = mean([AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.duration]);
        AllResults{sample, 'Data'}{1,1}{ii, 'duration_std'} = std([AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.duration]);
        AllResults{sample, 'Data'}{1,1}{ii, 'N_recordings'} = Nrecordings;
        
        %% calculate the averaged EOD for each sample and date
        [AllResults] = averaged_EODs_testosterone(AllResults, sample, date, Nrecordings);
        
        
        %% Calculate mean & sd of voltage peaks & ratios, & freq. of max power
        
        % Log normalized voltage peaks & ratios
        AllResults{sample, 'Data'}{1,1}.norm_vP0_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.norm_vP0}));
        AllResults{sample, 'Data'}{1,1}.norm_vP0_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.norm_vP0}));
        AllResults{sample, 'Data'}{1,1}.norm_vP1_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.norm_vP1}));
        AllResults{sample, 'Data'}{1,1}.norm_vP1_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.norm_vP1}));
        AllResults{sample, 'Data'}{1,1}.norm_vP2_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.norm_vP2}));
        AllResults{sample, 'Data'}{1,1}.norm_vP2_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.norm_vP2}));
        AllResults{sample, 'Data'}{1,1}.ratio_vP0_vP1_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.ratio_vP0_vP1}));
        AllResults{sample, 'Data'}{1,1}.ratio_vP0_vP1_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.ratio_vP0_vP1}));
        AllResults{sample, 'Data'}{1,1}.ratio_vP2_vP1_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.ratio_vP2_vP1}));
        AllResults{sample, 'Data'}{1,1}.ratio_vP2_vP1_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.ratio_vP2_vP1}));
       
        % Log non-normalized voltage peaks & Voltage range        
        AllResults{sample, 'Data'}{1,1}.No_norm_vP0_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_vP0}));
        AllResults{sample, 'Data'}{1,1}.No_norm_vP0_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_vP0}));
        AllResults{sample, 'Data'}{1,1}.No_norm_vP1_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_vP1}));
        AllResults{sample, 'Data'}{1,1}.No_norm_vP1_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_vP1}));        
        AllResults{sample, 'Data'}{1,1}.No_norm_vP2_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_vP2}));
        AllResults{sample, 'Data'}{1,1}.No_norm_vP2_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_vP2}));   
        AllResults{sample, 'Data'}{1,1}.No_norm_Vrange_mean(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_Vrange}));
        AllResults{sample, 'Data'}{1,1}.No_norm_Vrange_sd(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.No_norm_Vrange}));   
        
        % Log means and sd of frequency of maximum power
        AllResults{sample, 'Data'}{1,1}.PeakFreq_mean_Hz(ii) = mean(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.PeakFreq_Hz}));
        AllResults{sample, 'Data'}{1,1}.PeakFreq_sd_Hz(ii) = std(cell2mat({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.PeakFreq_Hz}));
        
        %% if any recording for this date and sample was head negative, log it 
        %make a cell array with all the unique messages in this column
        dummy1 = unique({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.fish_orientation});
        %are there both head positive and head negative recordings?
        if length(dummy1) > 1
            AllResults{sample, 'Data'}{1,1}.fish_orientation(ii) = 'both head positive & head negative';        
        %if there is only one orientation, find which one it was:
        %convert to matrix and find the positition of the character where the word 'negative'
        %begins (if it is not present, no number is given)
        elseif strfind(cell2mat(dummy1), 'negative') > 0
            AllResults{sample, 'Data'}{1,1}.fish_orientation(ii) = 'head negative';
        else
            AllResults{sample, 'Data'}{1,1}.fish_orientation(ii) = 'head positive';
        end

        %% if any recording for this date and sample was noisy, log it
        %make a cell array with all the unique messages in this column
        dummy2 = unique({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.threshold_crit});
        %are there both noisy and clean recordings?
        if length(dummy2) > 1
            AllResults{sample, 'Data'}{1,1}.noisy_recordings(ii) = 'some';
        %if only one threshold was used, find which one it was:
        %convert to matrix and find the positition of the character where the word '3 std_dev'
        %begins (if it is not present, no number is given)
        elseif strfind(cell2mat(dummy2), '3 std_dev') > 0
            AllResults{sample, 'Data'}{1,1}.noisy_recordings(ii) = 'yes';
        else
            AllResults{sample, 'Data'}{1,1}.noisy_recordings(ii) = 'no';
        end
        
        %% log the presence or absence of P3 for this date and sample
        %make a cell array with all the unique messages in this column
        dummy3 = unique({AllResults{sample, 'Data'}{1,1}{ii, 'recordings'}{1,1}.P3});
        %is there more than one state for P3?
        if length(dummy3) > 1
            AllResults{sample, 'Data'}{1,1}.P3(ii) = 'present & absent';
        %if there is only one state for P3
        else 
            AllResults{sample, 'Data'}{1,1}.P3(ii) = cell2mat(dummy3);
        end 
   
    end
 
end

%% If desired, uncomment below to export the table to csv

% matlab iterates through columns of the vector, hence the transposition
for m = AllSamples_mod'
        
    % Remove the following variables, they don't export properly (they
    % could be exported individually if desired)
    export_data = removevars(AllResults{m, 'Data'}{1,1}, {'recordings', ...
        'time_averaged_EOD', 'voltage_averaged_EOD', 'joint_recordings'});
    
    filename = sprintf("../../postMatlab_EOD/%s_EOD_data.csv", m);

    % export
    writetable(export_data, filename)
    
end

%export the table with the additional info about each recording
writetable(recordings_extra_details, ...
    '../../postMatlab_EOD/recordings_extra_details.csv')

%% old code
% for each sample, if there are any warnings, log this fact in the main table 
% for i = 1:length(AllSamples_mod)
%     if isfield(AllResults{AllSamples_mod{i}, 'Data'}{1,1}, 'warnings') == 1
%         AllResults{AllSamples_mod{i}, 'warnings'} = "yes";
%     end
% end 

%% Plot and save plots

% to compare EOD recordings from the same sample & date (=how variable are
% my recordings)
run Overlay_all_recordings_testosterone
% 
% to compare how the averaged EOD for a sample changes through time
% normalized voltages
run Plot_EOD_through_time_testosterone
% non-normalized voltages
run Plot_EOD_through_time_testosterone_Nonorm

%% consolidate the last day's recording's duration in an easy,
% quick-n-dirty, copy-paste manner: 
today_date = sprintf("%s_duration", dates{end});
today = table('Size', [Nsamples, 2], 'VariableNames', ["Fish", today_date], ...
    'VariableTypes', {'string', 'double'});

for s = 1:Nsamples
    
    today.Fish(s) = AllSamples_mod(s);
    
    % display only today's recordings
    if AllResults{AllSamples_mod(s), 'Data'}{1,1}.date(end) == dates{end}
        today{s, today_date} = AllResults{AllSamples_mod(s), 'Data'}{1,1}.duration_mean(end);
    else 
       today{s, today_date} = NaN;
    end
        
    %today.today_duration(s) = AllResults{AllSamples_mod(s), 'Data'}{1,1}.duration_mean(end);
end

rows2vars(today)
