%% plot_mode_locking_analysis
clear 

%% Load in spreadsheet

[base, ~, ~, ~] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Load in example 

% putative = 'R25_TT1_P8_N04'; % mode locking
% putative = 'R29_TT3_P2_N22'; % mode locking
% putative = 'R29_TT1_P2_N09'; % mode locking
% putative = 'R25_TT2_P9_N06'; % mode locking
putative = 'R25_TT1_P8_N12'; % mode locking
% putative = 'R29_TT1_P3_N01'; % no mode locking
% putative = 'R27_TT3_P1_N08'; % no phase locking

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

figure('Position',[200,22,403,1233])
tiledlayout(1, 2, 'Padding','compact')
num_fpeaks = length(data_ST.fpeaks);

% Plot period histogram
nexttile
max_rate = max(temporal.p_hist, [], 'all');
freq_values = round(param.fpeaks);
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


%% Plot ISI histogram

% Calculate ISI for each rep
nreps = param.nrep;
ISI_all = cell(1, num_fpeaks);
nbins = 30;
edges = linspace(0, 20, nbins+1); % Bin edges
counts_all = zeros(num_fpeaks, nbins); % Pre-allocate counts_all
for jj = 1:num_fpeaks
    x = temporal.x{jj} / 1000; % ms
    y = temporal.y{jj};

    ISI = arrayfun(@(ii) diff(x(y == ii)), 1:nreps, 'UniformOutput', false);
    ISI_all{jj} = vertcat(ISI{:});

    counts_all(jj, :) = histcounts(ISI_all{jj}, edges);
end


% PLot ISI histograms
nexttile
max_rate = max(counts_all, [], 'all');
freq_values = round(param.fpeaks);
for j = 1:num_fpeaks
	counts = counts_all(j,:);

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
grid on
title('ISI Histogram')
set(gca, 'fontsize', 12)


%% Calculate ISI scatter plot 

% Calculate ISI for each rep
ISI_all = cell(1, num_fpeaks);
counts_all = zeros(num_fpeaks, 50); % Pre-allocate counts_all
for jj = 1:num_fpeaks
	ISI = [];
	x_isi = [];
	y_isi = [];
	x = temporal.x{jj}/1000; % ms
	y = temporal.y{jj};
	for ii = 1:30
		valid = y==ii;
		isi = diff(x(valid));
		ISI = [ISI isi'];
		x_isi = [x_isi isi(1:end-1)'];
		y_isi = [y_isi isi(2:end)'];
	end
	ISI_all{jj} = ISI;
	x_isi_all{jj} = x_isi;
	y_isi_all{jj} = y_isi;
end


figure('Position',[560,42,1050,806])
tiledlayout(6, 7, 'Padding','compact')
for jj = 1:num_fpeaks
	T = 5; % period
	nexttile
	hold on
	line(repmat([0;10*T],1,10),[1;1]*(1:10)*T,'Color','r', 'linewidth', 0.3)
	line([1;1]*(1:10)*T,repmat([0;10*T],1,10),'Color','r', 'linewidth', 0.3)
	axis square
	scatter(x_isi_all{jj}, y_isi_all{jj},15, 'filled', 'k')

	if ismember(jj, 37:41)
		xlabel('ISI n (ms)')
	end

	if ismember(jj, 1:7:41)
		ylabel('ISI n+1 (ms)')
	end
	title(sprintf('%d Hz', freq_values(jj)))
	xlim([0 30])
	ylim([0 30])
end

%% Significance 
% An Ris score of 13.8 (Rayleigh criterion) and a z score of 2 is 
% required to meet our criterion for "mode-locking."

% Ris score 



% Z-score 

%% Calculate VS sync to 200 Hz
% 
% figure
% hold on
% plot(param.fpeaks, temporal.VS)
% plot(param.fpeaks, smooth_rates(temporal.VS,zeros(num_fpeaks, 1),...
% 	ones(num_fpeaks, 1), CF), 'k')
% xline(CF, '--')
% xlabel('Spectral Peak Freq. (Hz)')
% ylabel('Vector Strength')
% legend('VS', 'Smoothed VS')
% set(gca, 'fontsize', 12)
% title('Example neuron VS, BS')
 
%% Calculate FFT on PSTHs
temporal = analyzeST_Temporal(param, data_ST);

figure('Position',[200,22,403,1233])
tiledlayout(1, 1, 'Padding','compact')
num_fpeaks = length(data_ST.fpeaks);

spike_hist = temporal.PSTH;
edges = temporal.t;
t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
bin_width = (edges(2) - edges(1))/1000; % seconds 
fs = 1/bin_width;

[VS, fft_output, f] = calcFFT(spike_hist, param, fs);

max_rate = max(fft_output, [], 'all');
freq_values = round(param.fpeaks);
hold on
for j = 1:num_fpeaks % Plot FFTs
	offset = (j-1)*max_rate; % Adjust offset amount
	plot(f, fft_output(j,:) + offset);
end
ylim([0 max_rate*num_fpeaks])
xlabel('Frequency (Hz)')
box on
yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
yticklabels(freq_values)
grid on
title('FFT(PSTH)')
xlim([0 fs/2])
set(gca, 'fontsize', 12)

%% FUNCTIONS 

function [VS, fft_output, f] = calcFFT(spike_hist, params, fs)

freq = 200;
for ii = 1:length(params.fpeaks)

	% Cut to integer # of cycles (50 cycles, excluding onset)
	onsetwin = 0.05; % ms
	spike_wo_onset = spike_hist(ii, round(onsetwin*fs):end);

	% Find frequency index for 200 Hz
	L = length(spike_wo_onset);

	% Take FFT 
	Y = fft(spike_wo_onset);
	P2 = abs(Y)/L; % Normalize by signal length
	P1 = P2(1:L/2+1);
	P1(2:end-1) = 2*P1(2:end-1);
	f = fs/L*(0:(L/2));

	% Get FFT at 200 Hz 
	ind = f==freq;
	%VS(ii) = P1(ind);
	fft_output(ii,:) = P1;
	VS = 0;

end
end



