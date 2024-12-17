%% plot_binmode_example

%% Load in dataset and session

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

% Session of interest
putative = 'R24_TT2_P13_N05';


%% Plot

bin_colors = {'#3F985C', '#3F985C', '#3F985C', '#3F985C'};
contra_colors = {'#1267B9', '#1267B9', '#1267B9', '#1267B9'};

figure('Position',[70,607,951,250])
tiledlayout(1, 4)
linewidth = 1.5;
for ind = [1, 2, 4] %num_ds

	% Load in data
	[base, datapath, savepath, ppi] = getPaths();
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,filename)), 'data';
	index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
	CF = sessions.CF(index);

	% Plot 
	nexttile
	hold on
	for ibin = 1:2
		if ibin == 1
			colors = contra_colors;
		else
			colors = bin_colors;
		end

		% Analysis
		params = data(5+ind, 1:2);
		data_ST  = analyzeST(params);

		rate = data_ST{ibin}.rate;
		rate_std = data_ST{ibin}.rate_std;
		rlb = data_ST{ibin}.rlb;
		rub = data_ST{ibin}.rub;
		fpeaks = data_ST{ibin}.fpeaks;
		spl = data_ST{ibin}.spl;
		rate_sm = data_ST{ibin}.rates_sm;

		% Plot
		rates_sm = smooth_rates(rate, rlb, rub);
		errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
			'linestyle', 'none', 'linewidth', 0.8, 'color', colors{ind})
		plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',colors{ind})
		%label{1} = [num2str(params{order(ind)}.spl) ' dB SPL'];

		[~, peak_ind] = max(data_ST{ibin}.rate);
		peak_f(ibin) = data_ST{ibin}.fpeaks(peak_ind);
	end

	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
	xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	set(gca, 'Fontsize', 14, 'XTick', plot_range(1)+0.200:0.400:plot_range(2)-0.200);
	xlim(plot_range);
	grid on
	ylimit = ylim;
	ylim([0 ylimit(2)])
	ylabel('Avg. rate (sp/s)')
	xlabel('Spectral Peak  Frequency (Hz)')
	title(sprintf('%d dB SPL', data_ST{1}.spl))

	q_factor(ind,1:2) = [data_ST{1}.width/peak_f(1) data_ST{2}.width/peak_f(2)];
	q_factor2(ind,1:2) = [peak_f(1)/data_ST{1}.width peak_f(2)/data_ST{2}.width];


end

% Q-value
% q_factor(3,:) = [];
% nexttile
% plot(q_factor')
% title('Q')
% ylabel('Freq_p_e_a_k / HHBW')
% xticks([1, 2])
% xlim([0 3])
% xticklabels({'Contra', 'Bin'})
% xlabel('SPL')
% set(gca, 'Fontsize', 14)

% Q-value
q_factor2(3,:) = [];
nexttile
plot(q_factor2', 'o-')
title('Q_H_H_B_W')
ylabel('Freq_p_e_a_k / HHBW')
xticks([1, 2])
xlim([0 3])
xticklabels({'Contra', 'Bin'})
xlabel('SPL')
set(gca, 'Fontsize', 14)
legend('43', '63', '83', 'Location','best')
grid on


