function [AllResults] = averaged_EODs_testosterone(AllResults, sample, date, Nrecordings)
% *This code will average the voltage values for those time points with a
% recorded voltage in at least half the recordings.*

%get the index of the chosen date
idate = find(strcmp(date, AllResults{sample, 'Data'}{1,1}.date) == 1);

%% make a table with one column for the EOD_time and one column for each 
%EOD_voltage with the voltage values for a given time on the same row (some
%time values at the begining and end will have voltage values for only some
%recordings

%create empty table to hold the desired data
AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1} = ...
    table('Size', [1, 1],  'VariableTypes', {'double'}, 'VariableNames', {'time_all_eods'});

%iterate through all the recordings
for i = 1:Nrecordings
    
    % this is the column name for the voltages for EOD i
    volt_col = sprintf('voltages_eod_%i', i);
    
    %this is a dummy variable with the times and voltages of the EOD i. 
    %They must be arranged in a table for outerjoin to work. 
    %Rounding the time values is necessary because some values differ on
    %the ~15 decimal (I chose 10 digits to preserve very small numbers as
    %such
    dummy = table(round(AllResults{sample, 'Data'}{1,1}{idate, 'recordings'}{1,1}(i).time_eod, 10, 'significant'), ...
        AllResults{sample, 'Data'}{1,1}{idate, 'recordings'}{1,1}(i).voltages_eod, 'VariableNames', {'time_all_eods', volt_col});
    
    %this code outerjoins the tables of the already joint recordings (empty
    %at first) and EOD i. This means to arrange on the same row all voltage
    %values that share the same time value. NaN is added where no voltage
    %values were recorded for that time value (at the beginning or end of
    %the EOD).
    AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1} = ...
        outerjoin(AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1}, ...
        dummy, 'Keys', 'time_all_eods', 'MergeKeys', true);

end


%% I decided to only consider a time point if it has a voltage value for
%more than half of the EODs recorded

%the following line adds a column with the counts of these ocurrences
AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1}.total_volt_pts = ...
    Nrecordings - sum(ismissing(AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1}), 2);

%get the time values that meet the criterion above
AllResults{sample, 'Data'}{1,1}{idate, 'time_averaged_EOD'}{1,1} =  ...
    table2array(AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1}...
    (AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1}.total_volt_pts...
    >= Nrecordings/2, 'time_all_eods'));

%get the average voltage for the selected time values. I broke the code in
%pieces to make it more readable

%This flags the rows that meet the criterion above 
keep_rows = AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1}.total_volt_pts >= Nrecordings/2;

%This lists the columns that have the voltages values for the eods
keep_columns = contains(AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1}.Properties.VariableNames, 'voltages_eod');

%By using the previous 2 variables, this extracts the cells that need be
%averaged
keep_data = AllResults{sample, 'Data'}{1,1}{idate, 'joint_recordings'}{1,1} (keep_rows, keep_columns);

%this averages the voltage value for each of the selected time points and
%omits NaN values. Conversion to array is required by the mean function
AllResults{sample, 'Data'}{1,1}{idate, 'voltage_averaged_EOD'}{1,1} = mean(table2array(keep_data), 2, 'omitnan');


%% Calculate the duration of the averaged EOD
AllResults{sample, 'Data'}{1,1}.averagedEOD__duration(idate) = ...
    AllResults{sample, 'Data'}{1,1}{idate, 'time_averaged_EOD'}{1,1}(end)...
    - AllResults{sample, 'Data'}{1,1}{idate, 'time_averaged_EOD'}{1,1}(1);

