%% plot_temporal_examples 

%% Load in spreadsheet

[base, ~, ~, ~] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Load in example 

putative = 'R25_TT1_P8_N04'; % mode locking
%putative = 'R29_TT3_P2_N22'; % mode locking
%putative = 'R29_TT1_P2_N09'; % mode locking
%putative = 'R25_TT2_P9_N06'; % mode locking
%putative = 'R25_TT1_P8_N12'; % mode locking
% putative = 'R29_TT1_P3_N01'; % no mode locking
%putative = 'R27_TT3_P1_N08'; % no phase locking

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

%% Plot as heatmap 

% Plot period histogram
num_fpeaks = length(data_ST.fpeaks);
max_rate = max(temporal.p_hist, [], 'all');
freq_values = round(param.fpeaks);
p_hist = temporal.p_hist;
edges = temporal.t_hist;
t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers

grayMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', linspace(0, 1, 256)'];
grayMap = flipud(grayMap);


figure;
h = pcolor(t_bin, data_ST.fpeaks, p_hist);
set(h, 'EdgeColor', 'none');
hold on
yline(CF, 'r', 'LineWidth',3)
title('Example');
colorbar;
%shading interp
axis square;
colormap(grayMap);
%xlim([0 5])
clim([0 70])
ylabel('Spectral Peak Freq. (Hz)')
xlabel('Period (ms)')
