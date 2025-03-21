%% plot_example_temporal
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Load in example

putative = 'R27_TT2_P8_N02'; %'R24_TT2_P13_N05'; %
[base, datapath, savepath, ppi] = getPaths();
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

%% Analysis

% Synthetic timbre analysis
params = data(7, 2);
params = params(~cellfun(@isempty, params));
data_ST  = analyzeST(params, CF);
data_ST = data_ST{1};
param = params{1};
temporal = analyzeST_Temporal(param, data_ST);

%% Plot PSTH and smoothed PSTH

% fpeaks_ind = find(mod(param.fpeaks, param.Delta_F)==0);
% fpeaks_harm = param.fpeaks(fpeaks_ind);
% freq_values = round(fpeaks_harm);
max_rate = max(temporal.PSTH, [], 'all');
freq_values = round(param.fpeaks);

figure('Position',[200,22,403,1233])
tiledlayout(1, 2, 'Padding','compact')

% Plot PSTH full response
nexttile
num_fpeaks = length(freq_values);
hold on
for j = 1:num_fpeaks

	% Plot PSTHs
	counts = temporal.PSTH(j,:);
	edges = temporal.t;
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');

	% Plot smoothed PSTH
	plot(temporal.t(1:end-1), temporal.PSTH_smooth(j,:)+offset,'k', 'LineWidth',1.5);

end
ylim([0 max_rate*num_fpeaks])
xlabel('Time (ms)')
box on
yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
yticklabels(freq_values)
title('PSTH')
set(gca, 'fontsize', 12)

% Plot period histogram
nexttile
max_rate = max(temporal.p_hist, [], 'all');
for j = 1:num_fpeaks

	% Plot PSTHs
	counts = temporal.p_hist(j,:);
	edges = temporal.t_hist;
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
end
ylim([0 max_rate*num_fpeaks])
xlabel('Time (ms)')
box on
yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
yticklabels(freq_values)
xlim([0 5])
xticks(1:5)
grid on
title('Period Histogram')
set(gca, 'fontsize', 12)

%% Calculate VS sync to 200 Hz

figure
hold on
plot(param.fpeaks, temporal.VS)
plot(param.fpeaks, smooth_rates(temporal.VS,zeros(num_fpeaks, 1),...
	ones(num_fpeaks, 1), CF), 'k')
xline(CF, '--')
xlabel('Spectral Peak Freq. (Hz)')
ylabel('Vector Strength')
legend('VS', 'Smoothed VS')
set(gca, 'fontsize', 12)
title('Example neuron VS, BS')

