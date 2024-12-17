%% single_unit_examples_2.m
%
%
%


%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Set up figure
linewidth = 1.5;
figure('Position',[4,189,1079,718])
tiledlayout(3, 3)
data_colors = {'#03882F', '#82BB95'};


%% Plot

examples = {'R24_TT2_P13_N05', 'R27_TT2_P8_N02', 'R27_TT2_P8_N05', ...
	'R25_TT2_P9_N01', 'R27_TT3_P1_N08', 'R27_TT2_P7_N01', ...
	'R29_TT4_P5_N15', 'R25_TT2_P8_N02', 'R29_TT1_P3_N05'};

for ineuron = 1:9

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
	data_ST  = analyzeST(params);
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
	nexttile
	hold on
	rates_sm = smooth_rates(rate, rlb, rub);
	errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
		'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{1})
	plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{1})
	plot(fpeaks./1000, rates_sm, 'linewidth', linewidth, 'color', 'k')

	% 
	xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)
	yline(0.1, 'k', LineWidth=linewidth)

	% Figure parameters 
	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
	set(gca, 'Fontsize', 14, 'XTick', plot_range(1)+0.200:0.400:plot_range(2)-0.200);
	xlim(plot_range);
	grid on
	ylim([0 max_rate+5])
	xlabel('Spectral Peak Frequency (Hz)')
	ylabel('Avg. Rate (sp/s)')

	% Legend
	if ineuron == 9
		legend({'', 'Data', 'Estimated CF', 'Spont Rate'}, 'Location','southwest')
	end

	% Labeling MTF Type
	if ineuron == 9
		text(0.05, 0.85, MTF_shape, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',16)
	else
		text(0.05, 0.95, MTF_shape, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',16)
	end

	% Titles
	titles_x = {'Low CF', 'Medium CF', 'High CF'};
	locs = linspace(0.16, 0.772, 3);
	for ii = 1:3
		annotation('textbox',[locs(ii) 0.933 0.126 0.0459],...
			'String',titles_x{ii},'FontWeight','bold',...
			'FontSize',20,'EdgeColor','none');
	end

	titles_y = {'Flat', 'Dip', 'Peak'};
	locs = linspace(0.19, 0.78, 3);
	for ii = 1:3
		annotation('textbox',[0.045 locs(ii) 0.126 0.0459],...
			'String',titles_y{ii},'FontWeight','bold',...
			'FontSize',20,'EdgeColor','none', 'Rotation',90);
	end
end

%% Export figure

%exportgraphics(gcf, fullfile(savepath, 'manuscript', 'examples-peaksanddips.png'), 'Resolution', 600)
