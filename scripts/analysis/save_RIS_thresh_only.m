%% save_peak_finding_results.m
%
% Script that...
%
%
% Author: J. Fritzinger
% Created: 2026-03-11
%
% -------------------------------------------------------------------------
clear% --- Replacement for setupUmichClusters ---

%% Add paths 

addpath(genpath('/home/jofritzi/projects/synth-timbre/scripts/SPIKY'))
addpath(genpath('/home/jofritzi/projects/synth-timbre/scripts/helper-functions'))
addpath(genpath('/home/jofritzi/projects/shared-physio'), '-end')

%%

% Load in spreadsheet
[base, datapath, ~, ~] = getPaths();
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Create table
varNames = ["Putative", "CF", "MTF","MTF_at200", "MTF_str", ...
    "SPL", "binmode", "F0", ...
    "Type", "Prom", "Width", "Lim", "Freq", "Q", "Q_log"];
varTypes = ["string", "double", "string", "string", "double", ...
    "double", "double", "double", ...
    "string", "double", "double", "double", "double", "double", "double"];
est_num_rows = 160; % set to number larger than
num_cols = length(varNames);
table_size = [est_num_rows num_cols];
tables = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames);

% Create has_data
data_ind(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
data_ind(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
has_data = any(data_ind(:,2),2);
data_index = find(has_data);
num_neurons = sum(has_data);

%% Plot each neuron

row_idx = 1;
for isesh = 1:num_neurons
    ineuron = data_index(isesh);

    % Load in data
    putative = sessions.Putative_Units{ineuron};
    CF = sessions.CF(ineuron);
    MTF_shape = sessions.MTF{ineuron};
    at200 = sessions.MTF_at200{ineuron};
    load(fullfile(datapath, 'neural_data', [putative '.mat']))

    % Load 
    try
        load(fullfile(datapath, 'RIS_thresholds', [putative '_RIS.mat']))
    catch
        continue
    end
    nrep = param_ST{1}.nrep;

    %plot_SPIKE_Similarity(param_ST{1}.fpeaks, RI_S_dist, nrep)

    [threshold_percent, threshold_freq, d_prime] = calculate_SPIKE_Thresholds(param_ST{1}.fpeaks, RI_S_dist, nrep);

    % [threshold_percent, threshold_freq, d_prime] = dprime_from_spike_distance(RI_S_dist, nrep, param_ST{1}.fpeaks);
    % disp(max(max(d_prime)))
    % disp(threshold_percent)

    tables.Putative{row_idx} = sessions.Putative_Units{ineuron};
    tables.CF(row_idx) = CF;
    tables.D_prime(row_idx) = max(max(d_prime));
    tables.Threshold(row_idx) = threshold_percent;
    tables.Thresh_Freq{row_idx} = threshold_freq;
    row_idx = row_idx+1;
end
writetable(tables,fullfile(datapath, 'peak_picking_w_thresholds_RIS_real.xlsx'))


%% Analysis
figure
scattersize = 30;
fontsize = 10;

% Set up figure 
h(1) = subplot(1, 2, 1);

for ibin = 2
	for ispl = 2

		% Data
		CFs = tables.CF;
		Qs = tables.Threshold;

		% Add in units without thresholds
		Qs(isnan(Qs)) = 100;
		Qs(Qs>100) = 100;

		% Plot
		%gfig = gscatter(CFs/1000, Qs,peaks, 'filled');
		%set(gfig, 'MarkerEdgeColor', 'k');
		scatter(CFs/1000, Qs, scattersize, 'filled', 'MarkerEdgeColor','k', ...
			'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
		hold on

		% Plot minimum thresholds based on stimuli
		%plot(CFs_array/1000, min_threshold, 'r', 'LineWidth',2)

		% Plot human threshold
		scatter(1200/1000, 4,scattersize, 'r', 'filled','markeredgecolor', 'k') %'filled')
		yline(4, 'r')

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		%title('Thresholds vs CF')
		xlabel('CF (kHz)')
		if ispl == 1
			ylabel('Q')
		end
		ylim([0.35 100])
		xlim([0.3 10])
		set(gca, 'fontsize', fontsize)
		set(gca, 'XScale', 'log')
		set(gca, 'YScale', 'log')
		ylabel('Threshold from VS (%)')
		xticks([0 200 500 1000 2000 5000 10000]/1000)
		yticks([0.2 0.5 1 2 5 10 20 50 100])
		yticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '50', '>100'})
		grid on
		box off
		hleg = legend('Neural','Human', '', 'Location','southwest', 'box',...
			'off');
		hleg.ItemTokenSize = [8, 8];
	end
end

%% 

h(2) = subplot(1, 2, 2);

spls = [43, 63, 73, 83];
is200 = tables.F0==200;
isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');

for ibin = 2
	isbin = tables.binmode == ibin;
	for ispl = 2

		% Get data
		islevel = tables.SPL == spls(ispl);
		index = islevel & isbin & is200; % & isMTF;

		% Data
		Qs = tables.Threshold(index);
        Qs_VS = tables_VS.Threshold(index);

		% Add in units without thresholds
		Qs(isnan(Qs)) = 50;
		Qs(Qs>50) = 50;
        Qs_VS(isnan(Qs_VS)) = 50;
        Qs_VS(Qs_VS>50) = 50;

		% Plot
		scatter(Qs, Qs_VS, scattersize, 'filled', 'MarkerEdgeColor','k', ...
			'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
		hold on
        plot([0.35 50], [0.35 50], 'k')

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		%title('Thresholds vs CF')
		xlabel('Threshold from Rate (%)')
		if ispl == 1
			ylabel('Q')
		end
		ylim([0.35 50])
		xlim([0.35 50])
		set(gca, 'fontsize', fontsize)
		set(gca, 'XScale', 'log')
		set(gca, 'YScale', 'log')
		ylabel('Threshold from VS (%)')
		xticks([0.2 0.5 1 2 5 10 20 50 70])
		yticks([0.2 0.5 1 2 5 10 20 50 70])
		yticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '>50'})
        xticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '>50'})
		grid on
		box off
		hleg.ItemTokenSize = [8, 8];
	end
end
