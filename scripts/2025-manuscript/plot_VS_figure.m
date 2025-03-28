%% plot_temporal_examples

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Load in example

figure('Position',[1119,186,327,703])
tiledlayout(5, 2, "TileSpacing","compact", 'Padding','compact')
legsize = 9;
fontsize = 11;
labelsize = 18;

%%
for ii = 1:4

	switch ii
		case 1 
			%putative = 'R25_TT1_P8_N02'; % blurred at CF
			putative = 'R25_TT1_P8_N01'; % blurred at CF
			%putative = 'R25_TT4_P7_N01'; % blurred at CF
			%putative = 'R29_TT4_P2_N03'; % blurred at CF
		case 2
			putative = 'R29_TT4_P2_N09'; % 2 peaks 	
			%putative = 'R29_TT1_P3_N02'; % 2 peaks 
			%putative = 'R29_TT4_P2_N10'; % 2 peaks 
			%putative = 'R29_TT1_P3_N05';
		case 3
			%putative = 'R29_TT1_P2_N03'; % 2 peaks, odd
			putative = 'R25_TT1_P8_N15'; %'R25_TT1_P8_N03'; % 2 peaks, odd 
			%putative = 'R25_TT1_P8_N04'; %'R29_TT4_P2_N16'; % 2 peaks, odd 
			%putative = 'R25_TT2_P9_N02'; % % 2 peaks, odd
		case 4
			%putative = 'R29_TT3_P2_N05'; % multiple
			putative = 'R29_TT4_P3_N04'; % multiple
			%putative = 'R29_TT3_P2_N04'; % multiple
			%putative = 'R25_TT3_P8_N06'; % multiple
	end
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

	% Plot as heatmap
	num_fpeaks = length(data_ST.fpeaks);
	max_rate = max(temporal.p_hist, [], 'all');
	freq_values = round(param.fpeaks);
	p_hist = temporal.p_hist;
	edges = temporal.t_hist;
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	t = linspace(0, 5, size(p_hist,2));

	grayMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', linspace(0, 1, 256)'];
	grayMap = flipud(grayMap);

	nexttile
	hh = pcolor(t, data_ST.fpeaks./1000, p_hist);
	set(hh, 'EdgeColor', 'none');
	hold on
	yline(CF/1000, 'r', 'LineWidth',2)
	colormap(grayMap);
	xlim([0 5])
	max_rate = max(p_hist, [], "all");
	clim([0 max_rate-max_rate*0.3])
	if ismember(ii, 3)
		ylabel('                                Spectral Peak Freq. (kHz)')
	end
	xticks(0:5)
	if ismember(ii, [4, 8, 12, 16, 20])
		xlabel('Period (ms)')
	end
	box off

	set(gca, 'fontsize', fontsize)
	if ii == 17
		hleg = legend('', 'CF', 'Location','northwest', 'fontsize', ...
			legsize, 'box', 'off');
		hleg.ItemTokenSize = [8, 8];
	end
	msg = sprintf('%d sp/s', max_rate);
	text(0.58, 1, msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize, 'Color','r')

	% Plot vector strength calculation
	nexttile
	plot(param.fpeaks/1000, temporal.VS)
	hold on
	VS_smooth = smooth_rates(temporal.VS,zeros(num_fpeaks, 1),...
	 	ones(num_fpeaks, 1), CF);
	plot(param.fpeaks/1000, VS_smooth, 'k', 'linewidth', 2)
	xline(CF/1000, '--')
	xlim([param.fpeaks(1) param.fpeaks(end)]/1000)
	ylim([0 1])
	if ii == 4
		xlabel('Spectral Peak Freq. (kHz)')
	end
	if ismember(ii, 3)
		ylabel('                                Vector Strength')
	end
end

%% Load in spreadsheet
nexttile([1 2])

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name),...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Load in spreadsheet with peak information
spreadsheet_name = 'peak_picking_VS.xlsx';
table = readtable(fullfile(datapath, spreadsheet_name));

% Histogram overall 
fontsize = 12;

% Plot histogram  
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = table.binmode == 2;
for ispl = 1:4
	isSPL = table.SPL == spl(ispl);


	index = isSPL  & isBin;

	num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), table.Type(index)));
	num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), table.Type(index)));
	num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), table.Type(index)));
	all = sum([num_peak num_dip num_flat]);

	percent_peak(ispl) = num_peak;
	percent_dip(ispl) = num_dip;
	percent_flat(ispl) = num_flat;
	percent_all = [percent_peak; percent_dip; percent_flat]';

	bar(percent_all)
	%title([num2str(spl(ispl)) ' dB SPL'])
	xticklabels(spl)
	legend('Peak', 'Dip', 'Flat', 'Location','north', 'NumColumns', 3)
	ylabel('# Neurons')
	xlabel('Level (dB SPL)')
	ylim([0 150])
	set(gca, 'fontsize', fontsize)
	title('Peak/Dip Picking in VS Profiles')
end