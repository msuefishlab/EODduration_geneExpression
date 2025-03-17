
for i=1:length(averaged_eods)
norm_measurement_data(i)=standard_eod_measurement(normalized_eods(i).wave',normalized_eods(i).sampRate,normalized_eods(i).sample_name,normalized_eods(i).period,normalized_eods(i).treatment,normalized_eods(i).treatdate,normalized_eods(i).individual);
avg_measurement_data(i)=standard_eod_measurement(averaged_eods(i).wave',averaged_eods(i).sampRate,averaged_eods(i).sample_name,averaged_eods(i).period,averaged_eods(i).treatment,averaged_eods(i).treatdate,averaged_eods(i).individual);

    
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'units','normalized')

%% Plot All Normalized EODs
subplot(2,2,1)
        hold on;
            [~,n_waves]=size(normalized_eods(i).wave);
            for j=1:n_waves

                if isempty(normalized_eods(i).wave(:,j))
                    continue
                end
                time=linspace(0,1000*length(normalized_eods(i).wave(:,j))/normalized_eods(i).sampRate,length(normalized_eods(i).wave(:,j)));
                [~,ivmax]=max(normalized_eods(i).wave(:,j));
                time = time - (time(ivmax));
                plot(time,normalized_eods(i).wave(:,j)); %,'DisplayName',normalized_eods(i).info_text{j});
            end

            title(['N= ', num2str(n_waves) , ' EODs from ',string(normalized_eods(i).sample_name)],'Interpreter', 'none')
            hold off;

%% Plot All Normalized EODs, Focusing on P0    
    subplot(2,2,2)
        %% --Plot Decay Fit for Visualization ------ %%
    hold on;
    plot(avg_measurement_data(i).time_decay, avg_measurement_data(i).voltage_decay, 'o', 'DisplayName', 'Measured Data'); % Measured data
    plot(avg_measurement_data(i).time_decay, avg_measurement_data(i).fit_result(avg_measurement_data(i).time_decay), '-r', 'DisplayName', 'Decay Fit'); % Fitted curve

    plot(norm_measurement_data(i).time_decay, norm_measurement_data(i).voltage_decay, 'o', 'DisplayName', 'Measured Data'); % Measured data
    plot(norm_measurement_data(i).time_decay, norm_measurement_data(i).fit_result(norm_measurement_data(i).time_decay), '-r', 'DisplayName', 'Decay Fit'); % Fitted curve
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title('Decay Fit for P2');
    legend('show');
    grid on;
    hold off;

%% Plot the Averaged EOD with Landmarks
subplot(2,2,3)
    [~,ivmax]=max(averaged_eods(i).wave);
    time = time - (time(ivmax));
    plot(time,averaged_eods(i).wave);
    hold on;
    plot(avg_measurement_data(i).tT1,avg_measurement_data(i).vT1,'k+'); %begin and end time landmarks are always present
    text(avg_measurement_data(i).tT1+.05,avg_measurement_data(i).vT1,'\leftarrow EOD start', 'HorizontalAlignment','left')
    plot(avg_measurement_data(i).tT2,avg_measurement_data(i).vT2,'k+');
    text(avg_measurement_data(i).tT2+.05,avg_measurement_data(i).vT2,'\leftarrow EOD end', 'HorizontalAlignment','left')
    plot(avg_measurement_data(i).tP1,avg_measurement_data(i).vP1,'r+');					%phase1 landmarks (P1, S1, S2) always present
    text(avg_measurement_data(i).tP1+.05,avg_measurement_data(i).vP1,'\leftarrow P1', 'HorizontalAlignment','left')
    plot(avg_measurement_data(i).tS1,avg_measurement_data(i).vS1,'r+');
    text(avg_measurement_data(i).tS1+.05,avg_measurement_data(i).vS1,'\leftarrow Inflection 1', 'HorizontalAlignment','left')
    plot(avg_measurement_data(i).tS2,avg_measurement_data(i).vS2,'r+');
    text(avg_measurement_data(i).tS2+.05,avg_measurement_data(i).vS2,'\leftarrow Inflection 2', 'HorizontalAlignment','left')

    plot(avg_measurement_data(i).tP0,avg_measurement_data(i).vP0,'r+');
    text(avg_measurement_data(i).tP0+.05,avg_measurement_data(i).vP0,'\leftarrow P0', 'HorizontalAlignment','left','Rotation', -30)

    plot(avg_measurement_data(i).tZC1,avg_measurement_data(i).vZC1,'k+');
    text(avg_measurement_data(i).tZC1+.05,avg_measurement_data(i).vZC1,'\leftarrow Zero Crossing 1', 'HorizontalAlignment','left','Rotation', -30)
    plot(avg_measurement_data(i).tZC2,avg_measurement_data(i).vZC2,'k+');
    text(avg_measurement_data(i).tZC2+.05,avg_measurement_data(i).vZC2,'\leftarrow Zero Crossing 2', 'HorizontalAlignment','left', 'Rotation', -30)
    plot(avg_measurement_data(i).tP2,avg_measurement_data(i).vP2,'r+');
    text(avg_measurement_data(i).tP2+.05,avg_measurement_data(i).vP2,'\leftarrow P2', 'HorizontalAlignment','left')
    
    title(averaged_eods(i).sample_name,'Interpreter', 'none')
        
%% Plot Averaged EOD, with Focus on P0
subplot(2,2,4);
    [~,ivmax]=max(normalized_eods(i).wave);
    time = time - (time(ivmax));
    plot(time,normalized_eods(i).wave);
    hold on;
    plot(norm_measurement_data(i).tT1,norm_measurement_data(i).vT1,'k+'); %begin and end time landmarks are always present
    text(norm_measurement_data(i).tT1+.05,norm_measurement_data(i).vT1,'\leftarrow EOD start', 'HorizontalAlignment','left')
    plot(norm_measurement_data(i).tT2,norm_measurement_data(i).vT2,'k+');
    text(norm_measurement_data(i).tT2+.05,norm_measurement_data(i).vT2,'\leftarrow EOD end', 'HorizontalAlignment','left')
    plot(norm_measurement_data(i).tP1,norm_measurement_data(i).vP1,'r+');					%phase1 landmarks (P1, S1, S2) always present
    text(norm_measurement_data(i).tP1+.05,norm_measurement_data(i).vP1,'\leftarrow P1', 'HorizontalAlignment','left')
    plot(norm_measurement_data(i).tS1,norm_measurement_data(i).vS1,'r+');
    text(norm_measurement_data(i).tS1+.05,norm_measurement_data(i).vS1,'\leftarrow Inflection 1', 'HorizontalAlignment','left')
    plot(norm_measurement_data(i).tS2,norm_measurement_data(i).vS2,'r+');
    text(norm_measurement_data(i).tS2+.05,norm_measurement_data(i).vS2,'\leftarrow Inflection 2', 'HorizontalAlignment','left')

    plot(norm_measurement_data(i).tP0,norm_measurement_data(i).vP0,'r+');
    text(norm_measurement_data(i).tP0+.05,norm_measurement_data(i).vP0,'\leftarrow P0', 'HorizontalAlignment','left','Rotation', -30)

    plot(norm_measurement_data(i).tZC1,norm_measurement_data(i).vZC1,'k+');
    text(norm_measurement_data(i).tZC1+.05,norm_measurement_data(i).vZC1,'\leftarrow Zero Crossing 1', 'HorizontalAlignment','left','Rotation', -30)
    plot(norm_measurement_data(i).tZC2,norm_measurement_data(i).vZC2,'k+');
    text(norm_measurement_data(i).tZC2+.05,norm_measurement_data(i).vZC2,'\leftarrow Zero Crossing 2', 'HorizontalAlignment','left', 'Rotation', -30)
    plot(norm_measurement_data(i).tP2,norm_measurement_data(i).vP2,'r+');
    text(norm_measurement_data(i).tP2+.05,norm_measurement_data(i).vP2,'\leftarrow P2', 'HorizontalAlignment','left')
    
    title(averaged_eods(i).sample_name,'Interpreter', 'none')

%% Plot Power Spectrum
% subplot(2,3,6);
% a = semilogx(measurment_data.F,measurement_data.PdB,'b');
% b = get(a,'parent');
% set(b,'box','off', 'Fontname', 'Arial', 'Fontsize', 10);
% set(a,'LineWidth',1);
% xlabel ('frequency');
% title(['Power spectrum of current wave'], 'fontsize', 10);
% axis([0 10000 -60 0]);
% hold on;
% plot(measurment_data.Fmax,measurment_data.Pmax,'r+');
% plot(measurment_data.Flow,measurment_data.Plow,'r+');
% plot(measurment_data.Fhigh,measurment_data.Phigh,'r+');
    
%% Export PDF
pdfprinfig=gcf;
 pdffilename=['output_data' averaged_eods(i).sample_name];
 print(pdfprinfig, '-dpng', strjoin(pdffilename,"/")+".png"); % Printing tpdfprinfige figure to file
 
close(pdfprinfig)

end

%% Write Average Measurement File

% Get column names
colnames = fieldnames(avg_measurement_data);

% Dynamically filter out columns that contain numeric arrays
valid_fields = {};
for i = 1:numel(colnames)
    field_data = {avg_measurement_data.(colnames{i})}; % Extract all values for the field
    % Check if field contains numeric arrays (exclude them)
    if all(cellfun(@(x) ~(isstruct(x) || (isnumeric(x) && ~isscalar(x)) || isobject(x)), field_data))
        valid_fields{end+1} = colnames{i}; %#ok<AGROW>
    end
end

% Update column names to only valid fields
colnames = valid_fields;

% Extract rows for valid fields
filtered_data = arrayfun(@(x) cellfun(@(field) x.(field), colnames, 'UniformOutput', false), ...
    avg_measurement_data, 'UniformOutput', false);

% Convert the cell array of rows into a 2D cell matrix
filtered_data = vertcat(filtered_data{:});

% Add headers to the filtered data
ds_f = vertcat(colnames, filtered_data);

% Open file for writing
fid = fopen('output_data/avg_measurement_data.csv', 'wt');

if fid > 0
    for k = 1:size(ds_f, 1)
        if k == 1
            % Write header row
            fprintf(fid, '%s,', ds_f{k, 1:end-1});
            fprintf(fid, '%s\n', ds_f{k, end});
        else
            % Write data rows
            for j = 1:size(ds_f, 2)
                if isnumeric(ds_f{k, j}) || islogical(ds_f{k, j})
                    fprintf(fid, '%f,', ds_f{k, j});
                else
                    fprintf(fid, '%s,', string(ds_f{k, j}));
                end
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end


%% Write Normalized Measurement File

% Get column names
colnames = fieldnames(norm_measurement_data);

% Dynamically filter out columns that contain numeric arrays
valid_fields = {};
for i = 1:numel(colnames)
    field_data = {norm_measurement_data.(colnames{i})}; % Extract all values for the field
    % Check if field contains numeric arrays (exclude them)
    if all(cellfun(@(x) ~(isstruct(x) || (isnumeric(x) && ~isscalar(x)) || isobject(x)), field_data))
        valid_fields{end+1} = colnames{i}; %#ok<AGROW>
    end
end

% Update column names to only valid fields
colnames = valid_fields;

% Extract rows for valid fields
filtered_data = arrayfun(@(x) cellfun(@(field) x.(field), colnames, 'UniformOutput', false), ...
    norm_measurement_data, 'UniformOutput', false);

% Convert the cell array of rows into a 2D cell matrix
filtered_data = vertcat(filtered_data{:});

% Add headers to the filtered data
ds_f = vertcat(colnames, filtered_data);

% Open file for writing
fid = fopen('output_data/norm_measurement_data.csv', 'wt');

if fid > 0
    for k = 1:size(ds_f, 1)
        if k == 1
            % Write header row
            fprintf(fid, '%s,', ds_f{k, 1:end-1});
            fprintf(fid, '%s\n', ds_f{k, end});
        else
            % Write data rows
            for j = 1:size(ds_f, 2)
                if isnumeric(ds_f{k, j}) || islogical(ds_f{k, j})
                    fprintf(fid, '%f,', ds_f{k, j});
                else
                    fprintf(fid, '%s,', string(ds_f{k, j}));
                end
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end

%% Clean Up
clearvars -except normalized_eods averaged_eods subjectlist norm_measurement_data avg_measurement_data