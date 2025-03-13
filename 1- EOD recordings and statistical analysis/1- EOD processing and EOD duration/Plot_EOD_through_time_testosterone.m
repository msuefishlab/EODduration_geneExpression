%% This code will plot the averaged EOD for all recording dates, for each sample.

% this code helps visualize how the EOD for a given sample changes through
% time

% It requires the AllResults file
% If only some samples or dates are desired, adjust parameters accordingly

%% Choose samples 

% list all samples in the results file
fish = AllResults.Properties.RowNames;

% If only some samples are desired, build the array here:
% fish = [fish(1), fish(3)];

%% Define how many subplots are needed
% this function comes from
% https://www.mathworks.com/matlabcentral/fileexchange/26310-numsubplots-neatly-arrange-subplots
% [dims, ~] = numSubplots(length(fish)); 

%% Choose recording dates for each sample

for i = 1:length(fish)
    dates_fish = AllResults{fish{i}, 'Data'}{1,1}.date;
    
% If only some recording dates are desired, build the array here 
    % dates_fish = [dates_fish(1), dates_fish(7), dates_fish(13)];

%(Note: if very specific sample/sate combinations are desired those will
%have to be specified individually)
    
%% For fish i, iterate through its recording dates and overlay plots of averaged EODs
    for j = 1:length(dates_fish)
        
        %find the row in the table for fish i that contains the date at
        %hand
        jdate = find(strcmp(AllResults{fish{i}, 'Data'}{1,1}.date, dates_fish(j)));
        
        % choose a subplot
        % subplot(dims(1), dims(2), i)
        hold on
                
        %plot averaged_EOD
        time_av = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'time_averaged_EOD'}{1,1};
        volt_av = AllResults{fish{i}, 'Data'}{1,1}{jdate, 'voltage_averaged_EOD'}{1,1};
        plot(time_av, volt_av, 'LineWidth', 2, 'DisplayName', sprintf("%s", dates_fish{j}));
        
        %add lines at x = 0, y = 0, and y = +/- 0.05
        xline(0, 'HandleVisibility', 'off');
        yline(0, 'LineWidth', 0.4, 'HandleVisibility', 'off');
        yline(-0.005, ':', 'HandleVisibility', 'off');
        yline(0.005, ':', 'HandleVisibility', 'off');
        
        % add axis labels
        xlabel('time (ms)')
        ylabel('normalized Voltage (V)')
        
        % add title
        title(sprintf("%s", fish{i}), 'Interpreter','none');
        
        % add legend
        leg = legend('Interpreter','none', 'Location', 'southwest', 'FontSize', 8);
        leg.ItemTokenSize = [5, 18];
        
        %extend x and y axes, for aesthetics
        xlim([-2, 2]);
        ylim([-1, 1]);
    
        hold off
    end
    
        %% Save as figure
    
    fig = gcf;
    
    % save as matlab fig
    filename = sprintf("../../EOD_plots/between_rec_dates/Normalized/%s.fig", fish{i});
    savefig(fig, filename, 'compact');
    
    % export as picture
    filename = sprintf("../../EOD_plots/between_rec_dates/Normalized/%s.png", fish{i});
    width = 2000;
    height = 1000;
    fig.Position = [0 0 width height];
    print(filename, '-dpng', '-r0')
       
    close(fig);
    
end