%% methods_peak_quantification 

%% Load in spreadsheet 

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Set up figure 

figure('Position',[830,549,800,263])
fontsize = 14;

%% Peak/Dip/Sloping Examples 
 
examples = {'R25_TT3_P9_N01', 'R27_TT3_P1_N08', 'R29_TT1_P2_N04'};
CF_color = [0.3 0.3 0.3];
% Or R27_TT2_P8_N02 for peak

for ineuron = 1:3

	% Load in examples
	putative = examples{ineuron};
	[base, datapath, savepath, ppi] = getPaths();
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename)), 'data';
	index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
	CF = sessions.CF(index);
	MTF_shape = sessions.MTF{index};

	% Analysis
	param_ST = data(7, 2);
	data_ST = analyzeST(param_ST, CF);
	data_ST = data_ST{1};

	% Z-score
	rate = zscore(data_ST.rate);
	rate_sm = zscore(data_ST.rates_sm);

	% Cut down to +/- one octave range
	hi_limit = CF*2;
	lo_limit = CF/2;
	indices = data_ST.fpeaks > lo_limit & data_ST.fpeaks < hi_limit;
	
	% Calculate the peak/dip/flat
	[peaks, dips, type, prom, ~, lim, bounds_freq, halfheight] = peakFinding(data_ST, CF);

	% Plots
	h(ineuron) = subplot(2, 3, ineuron);
	patch([lo_limit lo_limit hi_limit hi_limit]./1000,[-4 4 4 -4], 'r', 'FaceAlpha',0.05, 'EdgeColor', 'none');
	hold on
	plot(data_ST.fpeaks./1000,rate, 'linewidth', 0.9, 'Color',"#0072BD");
	errorbar(data_ST.fpeaks./1000,rate, zscore(data_ST.rate_std)/sqrt(30), 'linewidth', 0.9, 'Color',"#0072BD");
	plot(data_ST.fpeaks./1000,rate_sm, 'linewidth', 1.5,'Color','k');
	ylim([-4 4])

	scatter(peaks.locs./1000, peaks.pks,50,  'filled', 'r')
	num_peaks = length(peaks.pks);
	scatter(dips.locs./1000, -1*dips.pks, 50, 'filled', 'r')
	num_dips = length(dips.pks);
	if ineuron == 1
		line([bounds_freq(1)/1000, bounds_freq(2)/1000], [halfheight, halfheight], 'Color', 'g', 'LineWidth', 1.5);
		line([peaks.locs peaks.locs]./1000, [peaks.pks-0.75 peaks.pks], 'Color', 'r', 'LineWidth', 1.5);
		xline(CF./1000, '--', 'Color',CF_color, 'linewidth', 1.5)
		hleg = legend('ROI', '', 'Data', 'Smoothed', 'Ref. Value', '', '+/- 0.75', ...
 			'Bandwidth', 'CF', 'Location','northeastoutside');
	elseif ineuron == 2
		line([dips.locs dips.locs]./1000, -1*[dips.pks-0.75 dips.pks], 'Color', 'r', 'LineWidth', 1.5);
		line([bounds_freq(1)/1000, bounds_freq(2)/1000], [halfheight, halfheight], 'Color', 'g', 'LineWidth', 1.5);
		xline(CF./1000, '--', 'Color',CF_color, 'linewidth', 1.5)
	else
		xline(CF./1000, '--', 'Color',CF_color, 'linewidth', 1.5)
	end

	plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)]./1000;
	
	xlim(plot_range)
	if ineuron == 1
		ylabel('Z-score')
		title('Peak')
	elseif ineuron == 2
		yticklabels([])
		title('Dip')
		xlabel('Spectral Peak Freq. (kHz)')
	else 
		yticklabels([])
		title('Sloping')
	end
	grid on
	set(gca, 'fontsize', fontsize)
end

%% Arrange and annotate 

left = repmat(linspace(0.07, 0.6, 3), 1, 2);
bottom = 0.17;
height = 0.7;

for ii = 1:3
	set(h(ii), 'position', [left(ii) bottom 0.23 height])
end
set(hleg, 'Position', [0.844510314761978,0.461977186311787,0.125,0.41254752851711])

% Add labels 
labelsize = 24;
annotation('textbox',[left(1)-0.03 0.95 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2)-0.03 0.95 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(3)-0.03 0.95 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
