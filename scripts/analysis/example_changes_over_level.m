% Fig___ How responses change/don't change over level
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Set up figure
figure('Position',[77,128,1080,608])
tiledlayout(2, 3, 'TileSpacing','compact')
linewidth = 1.5;

%% Plot
for ineuron = 1:6
	switch ineuron
		case 1
			putative = 'R24_TT2_P13_N03'; % BS
			%putative = 'R24_TT1_P12_N01'; % Low CF
		case 2
			putative = 'R27_TT2_P8_N02'; % BS
			%putative = 'R27_TT3_P1_N01'; % Low CF
		case 3
			putative = 'R24_TT2_P13_N05'; % BS
			%putative = 'R27_TT2_P8_N02'; % Med CF
		case 4
			putative = 'R27_TT3_P7_N08'; % BE
			%putative = 'R27_TT2_P8_N03'; % Med CF
		case 5
			putative = 'R27_TT3_P7_N14'; % BE
			%putative = 'R27_TT4_P8_N10'; % High CF
		case 6
			putative = 'R29_TT3_P5_N10'; % BE
			%putative = 'R25_TT4_P8_N05'; % High CF
	end

	% Load in data
	[base, datapath, savepath, ppi] = getPaths();
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename)), 'data';
	index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
	CF = sessions.CF(index);

	params = data(6:9, 2);
	params = params(~cellfun(@isempty, params));
	data_ST  = analyzeST(params, CF);

	% RM to get spont
	params_RM = data{2, 2};
	data_RM = analyzeRM(params_RM);
	spont = data_RM.spont;

	num_ds = size(data_ST, 2);
	if num_ds == 4
		%data_colors = {'#1b9e77', '#d95f02', '#7570b3', '#e7298a'};
		%data_colors = {'#ffffcc', '#a1dab4', '#41b6c4', '#225ea8'};
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
	else
		%data_colors = {'#1b9e77', '#d95f02', '#e7298a'};
		%data_colors = {'#a1dab4', '#41b6c4', '#225ea8'};
		data_colors = {'#034E1C', '#03882F', '#82BB95'};
	end

	% Sort
	spls = cell2mat(cellfun(@(p) p.spl, data_ST, 'UniformOutput',false));
	[~, order] = sort(spls);
	order = fliplr(order);
	max_rate = max(cellfun(@(d) max(d.rate), data_ST));
	nexttile
	hold on
	label_ind = 1;
	for ind = 1:num_ds

		rate = data_ST{order(ind)}.rate;
		rate_std = data_ST{order(ind)}.rate_std;
		rlb = data_ST{order(ind)}.rlb;
		rub = data_ST{order(ind)}.rub;
		fpeaks = data_ST{order(ind)}.fpeaks;
		spl = data_ST{order(ind)}.spl;
		rate_sm = data_ST{order(ind)}.rates_sm;

		% Plot
		rates_sm = smooth_rates(rate, rlb, rub, CF);
		errorbar(fpeaks, rate, rate_std/sqrt(30), 'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{ind})
		plot(fpeaks, rate, 'LineWidth',linewidth, 'Color',data_colors{:,ind})
		%plot(fpeaks, rate_sm, 'linewidth', linewidth, 'color', data_colors{ind})
		label{1} = [num2str(params{order(ind)}.spl) ' dB SPL'];
		%label_ind = label_ind+1;

		plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)];
		xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
		%label{2} = 'Estimated CF';
		yline(spont, 'k', LineWidth=linewidth)
		%yline(data_RM.spont, 'Color','k', 'LineWidth',2)
		%label(label_ind+1) = {'Spont'};
	end
	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)];
	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	label{label_ind} = 'Estimated CF';
	%yline(data_RM.spont, 'Color','k', 'LineWidth',2)
	%label(label_ind+1) = {'Spont'};
	xlabel('Spectral Peak  Frequency (Hz)')
	ylabel('Avg. rate (sp/s)')
	set(gca, 'Fontsize', 14, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
	xlim(plot_range);
	grid on
	%title('Synthetic Timbre', 'FontSize',18);
	%legend(label)

end

%% Export figure

%exportgraphics(gcf, fullfile(savepath, 'Fig2_SingleUnitExamples.png'), 'Resolution', 600)
