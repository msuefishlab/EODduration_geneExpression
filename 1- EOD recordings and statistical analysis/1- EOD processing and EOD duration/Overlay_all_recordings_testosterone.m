%% This code will plot every recording of every date for each sample, with its averaged EOD

% It requires the AllResults file
% If only some samples or dates are desired, adjust parameters accordingly

%% Choose samples 

% list all samples in the results file
fish = AllResults.Properties.RowNames;

% If only some samples are desired, build the array here:
% fish = [fish(1), fish(3)];

%% Choose recording dates for each sample

for i = 1:length(fish)
    dates_fish = AllResults{fish{i}, 'Data'}{1,1}.date;
    
% If only some recording dates are desired, build the array here 
    %dates_fish = [dates_fish(1), dates_fish(7), dates_fish(13)];

%(Note: if very specific sample/sate combinations are desired those will
%have to be specified individually)
    
%% Define how many subplots are needed
% this function comes from
% https://www.mathworks.com/matlabcentral/fileexchange/26310-numsubplots-neatly-arrange-subplots
[dims, ~] = numSubplots(length(dates_fish));

%% For fish i, iterate through its recording dates
    for j = 1:length(dates_fish)
        
        %find the row in the table for fish i that contains the date at
        %hand
        jdate = find(strcmp(AllResults{fish{i}, 'Data'}{1,1}.date, dates_fish(j)));
        %find how many recordings were made for that fish/date
        Nrecordings = AllResults{fish{i}, 'Data'}{1,1}.N_recordings(jdate);

%% Overlay plots of the k recordings and EODs for the chosen fish/date
       
        % choose a subplot
        subplot(dims(1), dims(2), j)
        
        % plot averaged_EOD
        time_av = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'time_averaged_EOD'}{1,1};
        volt_av = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'voltage_averaged_EOD'}{1,1};
        plot(time_av, volt_av, 'k', 'LineWidth', 5, 'DisplayName', 'averaged_EOD')  
        
        hold on
        
        for k = 1:Nrecordings
            % choose recording and EOD times & voltages
            time_all = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'recordings'}{1,1}(k).time_all;
            volt_all = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'recordings'}{1,1}(k).voltages_all;
            time_eod = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'recordings'}{1,1}(k).time_eod;
            volt_eod = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'recordings'}{1,1}(k).voltages_eod;
            % plot           
            plot(time_all, volt_all, 'LineWidth', 0.8, 'HandleVisibility', 'off')
            plot(time_eod, volt_eod, 'LineWidth', 1.5, 'HandleVisibility', 'off')
            
                      
        end
        
        %add lines at x = 0, y = 0, and y = +/- 0.05
        xline(0, 'HandleVisibility', 'off');
        yline(0, 'LineWidth', 0.4, 'HandleVisibility', 'off');
        yline(-0.005, ':', 'HandleVisibility', 'off');
        yline(0.005, ':', 'HandleVisibility', 'off');
               
        % add axis labels
        xlabel('time (ms)')
        ylabel('normalized Voltage (V)')
        
        % add title
        title(sprintf("%s - %s", fish{i}, dates_fish{j}), 'Interpreter','none');
        
        % add legend
        leg = legend('Interpreter','none', 'Location', 'southwest', 'FontSize', 8);
        leg.ItemTokenSize = [5, 18];
        
        hold off
    end
    
    %% Save as figure
    
    fig = gcf;
    
    % save as matlab fig
    filename = sprintf("../../EOD_plots/within_rec_date/%s.fig", fish{i});
    savefig(fig, filename, 'compact');
    
    % export as picture
    filename = sprintf("../../EOD_plots/within_rec_date/%s.png", fish{i});
    width = 2000;
    height = 1000;
    fig.Position = [0 0 width height];
    print(filename, '-dpng', '-r0')
       
    close(fig);

end  
    