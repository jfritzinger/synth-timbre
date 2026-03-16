

%% Load 

[base, ~, ~, ppi] = getPaths();
datapath = fullfile(base, 'data', '2025-manuscript-sub2');
tables_VS = readtable(fullfile(datapath, "peak_picking_w_thresholds_VS.xlsx"));

datapath = fullfile(base, 'data', '2025-manuscript');
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name),...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% 

% Create has_data
data_ind(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);


has_data = data_ind(:,1);
data_index = find(has_data);
num_neurons = sum(has_data);

% Plot each neuron
for isesh = progress(1:num_neurons)
    ineuron = data_index(isesh);

    % Load in session
    putative = sessions.Putative_Units{ineuron};
    CF = sessions.CF(ineuron);
    MTF_shape = sessions.MTF{ineuron};
    load(fullfile(datapath, 'neural_data', [putative '.mat']))

    if data_ind(ineuron,1)==1

        % Load in proper dataset for each idata
        param_ST = data(5+2, 2);

        % Analyze synthetic timbre
        data_ST = analyzeST(param_ST, CF);
        data_ST = data_ST{1};
        [~, peak_ind] = max(data_ST.rate);
        peak_f = data_ST.fpeaks(peak_ind);
        [~, peaksm_ind] = max(data_ST.rates_sm);
        peak_fsm = data_ST.fpeaks(peaksm_ind);

        temporal = analyzeST_Temporal(param_ST{1}, data_ST);
        VS(isesh, 1) = max(temporal.VS);

        % Get significant
        p_vals = temporal.VS_p_harms;
        is_sig = p_vals < 0.05;
        sig_harms = sum(is_sig,1);

        if sig_harms(1) > 20
            sig_200(isesh) = 1;
        else
            sig_200(isesh) = 0;
        end

        sig_all_harms(isesh,:) = sig_harms(2:end) > 20;
        if any(sig_harms(2:end) > 20)
            sig_harm(isesh) = 1;
        else
            sig_harm(isesh) = 0;
        end

        % Changes over spectral peak frequencies 
        %smoothed_max = smooth(max(temporal.VS_harms, [], 2));
        %[max_VS, max_ind] = max(smoothed_max);
        smooth_VS = smooth(temporal.VS);
        [~, max_ind] = max(smooth(temporal.VS));
        max_freq = param_ST{1}.fpeaks(max_ind);
        if max(smooth_VS) - min(smooth_VS) < 0.2
            near_peak(isesh) = 0;
        elseif max_freq < CF-600 || max_freq > CF+600
            near_peak(isesh) = 1;
        else
            near_peak(isesh) = 2;
        end

        % figure
        % plot(param_ST{1}.fpeaks, smooth_VS)
        % hold on
        % xline(CF)
        % title(near_peak(isesh))

        % Find peaks & prominence values
        [peaks, dips, type{isesh}, prom, width, lim, ~, ~, freq] = peakFinding(...
            data_ST, CF, 'Temporal', param_ST{1});
    end
    %fprintf('%s done, %d percent done\n', putative, round(isesh/num_neurons*100))
end

%% Analyze 

fprintf('Significant phase locking to 200 Hz: %0.1f%% (%d/%d)\n', ...
    sum(sig_200)/length(sig_200)*100, sum(sig_200), length(sig_200))
fprintf('Significant phase locking to harmonics: %0.1f%% (%d/%d)\n', ...
    sum(sig_harm)/length(sig_harm)*100, sum(sig_harm), length(sig_harm))

fprintf('Significant phase locking to both: %0.1f%% (%d/%d)\n', ...
    sum(sig_harm&sig_200)/length(sig_harm)*100, sum(sig_harm&sig_200), length(sig_harm))
sig_only_200 = sig_200 & ~sig_harm;
fprintf('Significant phase locking to 200Hz NOT harmonics: %0.1f%% (%d/%d)\n', ...
    sum(sig_only_200)/length(sig_harm)*100, sum(sig_only_200), length(sig_harm))
sig_only_harm = ~sig_200 & sig_harm;
fprintf('Significant phase locking to harmonics NOT 200Hz: %0.1f%% (%d/%d)\n', ...
    sum(sig_only_harm)/length(sig_harm)*100, sum(sig_only_harm), length(sig_harm))
not_sig = ~(sig_200|sig_harm);
fprintf('No phase locking: %0.1f%% (%d/%d)\n', ...
    sum(not_sig)/length(sig_harm)*100, sum(not_sig), length(sig_harm))

%% Per harm

nharm_sig = sum(sig_harm);
for ii = 1:9

    num_sig = sum(sig_all_harms(:,ii));
    fprintf('Significant phase locking to %dHz: %0.1f%% (%d/%d)\n', ...
    ii*200+200, num_sig/nharm_sig*100, num_sig, nharm_sig)

end

%% 

num_sig = sum(near_peak==0);
fprintf('No changes: %0.1f%% (%d/%d)\n', ...
    num_sig/length(sig_200)*100, num_sig, length(sig_200))
num_sig = sum(near_peak==2);
fprintf('Near CF: %0.1f%% (%d/%d)\n', ...
    num_sig/length(sig_200)*100, num_sig, length(sig_200))
num_sig = sum(near_peak==1);
fprintf('Away from CF: %0.1f%% (%d/%d)\n', ...
    num_sig/length(sig_200)*100, num_sig, length(sig_200))