%% single_unit_examples_2.m
%
%
%
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Set up figure

linewidth = 1.5;
figure('Position',[4,295,796,612])
%tiledlayout(3, 3)
data_colors = {'#03882F', '#82BB95'};
legsize = 14;

%% Plot example 

examples = {'R24_TT2_P13_N05', 'R27_TT2_P8_N02', 'R27_TT2_P8_N05', ...
	'R25_TT2_P9_N01', 'R27_TT3_P1_N08', 'R27_TT2_P7_N01', ...
	'R29_TT4_P5_N15', 'R25_TT2_P8_N02', 'R29_TT1_P3_N05'};

ineuron = 1;

% Load in data
putative = examples{ineuron};
[base, datapath, savepath, ppi] = getPaths();
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

% RM to get spont
params_RM = data{2, 2};
data_RM = analyzeRM(params_RM);
spont = data_RM.spont;

% Synthetic timbre analysis
params = data(7, 2);
params = params(~cellfun(@isempty, params));
data_ST  = analyzeST(params, CF);
data_ST = data_ST{1};
rate = data_ST.rate;
rate_std = data_ST.rate_std;
rlb = data_ST.rlb;
rub = data_ST.rub;
fpeaks = data_ST.fpeaks;
spl = data_ST.spl;
rate_sm = data_ST.rates_sm;
max_rate = max(rate);

% Plot
fpeaks_re_CF = log2(fpeaks/CF);
hold on
rates_sm = smooth_rates(rate, rlb, rub, CF);
errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
	'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{1})
plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{1})
plot(fpeaks./1000, rates_sm, 'linestyle', '-', 'linewidth', linewidth, 'color', 'k')
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)

% Figure parameters
plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
set(gca, 'Fontsize', 14)
xlim(plot_range);
grid on
ylim([0 max_rate+5])
ylabel('Avg. Rate (sp/s)')
xlabel('Spectral Peak Freq. (Hz)')

if ineuron == 1 || ineuron == 4 || ineuron == 7
	xlim([0.4 2.4])
elseif ineuron == 2 || ineuron == 5 || ineuron == 8
	xlim([0.8 3.2])
else
	xlim([2.8 9.2])
end

%% Calculate thresholds 

[threshold_percent, threshold_freq, slope_rate] = calculateThresholds(fpeaks, rate, rate_std, CF);
plot(threshold_freq/1000, slope_rate, 'r', 'LineWidth',3)




