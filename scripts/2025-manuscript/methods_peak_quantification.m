%% methods_peak_quantification 

%% Load in spreadsheet 

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Set up figure 

figure('Position',[830,392,596,420])
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
	errorbar(data_ST.fpeaks./1000,rate, 1/sqrt(30), 'linewidth', 0.9, 'Color',"#0072BD");
	plot(data_ST.fpeaks./1000,rate_sm, 'linewidth', 1.5,'Color','k');
	ylim([-4 4])
	scatter(peaks.locs./1000, peaks.pks,50,  'filled', 'r')
	num_peaks = length(peaks.pks);
	% for ip = 1:num_peaks
	% 	plot([peaks.locs(ip) peaks.locs(ip)]./1000, [peaks.pks(ip)-peaks.p(ip) peaks.pks(ip)], 'Color',"#D95319")
	% end
	scatter(dips.locs./1000, -1*dips.pks, 50, 'filled', 'g')
	num_dips = length(dips.pks);
	% for ip = 1:num_dips
	% 	plot([dips.locs(ip) dips.locs(ip)]./1000, -1*([dips.pks(ip)-dips.p(ip) dips.pks(ip)]), 'Color',"#7E2F8E")
	% end
	%plot(bounds_freq, [halfheight halfheight], 'g')

	plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)]./1000;
	xline(CF./1000, '--', 'Color',CF_color, 'linewidth', 1.5)
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

%% Q-value examples 
data_colors = {'#03882F', '#82BB95'};
linewidth = 1.5;
examples = {'R27_TT2_P8_N02', 'R27_TT3_P1_N08'};
for ineuron = 1:2

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

	[peaks, dips, type, prom, ~, lim, bounds_freq, halfheight] = peakFinding(data_ST, CF);
	rate = zscore(data_ST.rate);
	rate_sm = zscore(data_ST.rates_sm);

	% Plot
	if ineuron == 1
		h(3+ineuron) = subplot(2, 3, 3+ineuron);
		hold on
		yline(0)
		plot(fpeaks./1000,rate,'LineWidth', 1.5)
		plot(fpeaks./1000, rate_sm, 'linewidth', linewidth, 'color', 'k')
		line([bounds_freq(1)/1000, bounds_freq(2)/1000], [halfheight, halfheight], 'Color', 'g', 'LineWidth', 1.5);
		scatter(peaks.locs/1000, peaks.pks, 'r', 'filled', 'LineWidth',1.5)
		line([peaks.locs peaks.locs]./1000, [peaks.pks-0.75 peaks.pks], 'Color', 'r', 'LineWidth', 1.5);
		xline(CF/1000, '--')
		xlabel('Spectral Peak Freq. (kHz)')
		ylabel('Z-score')
		xlim([fpeaks(1)/1000 fpeaks(end)/1000])
		box on
		grid on
	else
		h(3+ineuron) = subplot(2, 3, 3+ineuron);
		hold on
		yline(0)
		plot(fpeaks./1000,rate,'LineWidth', 1.5)
		plot(fpeaks./1000, rate_sm, 'linewidth', linewidth, 'color', 'k')
		scatter(dips.locs/1000, -1*dips.pks, 'r', 'filled', 'LineWidth',1.5)
		line([dips.locs dips.locs]./1000, -1*[dips.pks-0.75 dips.pks], 'Color', 'r', 'LineWidth', 1.5);
		line([bounds_freq(1)/1000, bounds_freq(2)/1000], [halfheight, halfheight], 'Color', 'g', 'LineWidth', 1.5);
		xline(CF/1000, '--')
		xlabel('Spectral Peak Freq. (kHz)')
		%ylabel('Z-score')
		yticklabels([])
		xlim([fpeaks(1)/1000 fpeaks(end)/1000])
		box on
		grid on
		legend('', 'Data', 'Smoothed', 'Ref. Value', '+/- 0.75', ...
			'Bandwidth', 'CF', 'Location','northeastoutside')
	end
	set(gca, 'fontsize', fontsize)
end

%% Arrange and annotate 

left = repmat(linspace(0.1, 0.7, 3), 1, 2);
left(5) = 0.445;
bottom = [0.6 0.6 0.6 0.11 0.11];
height = 0.32;

for ii = 1:3
	set(h(ii), 'position', [left(ii) bottom(ii) 0.28 height])
end

for ii = 4:5
	set(h(ii), 'position', [left(ii) bottom(ii) 0.33 height])
end


% Add labels 
labelsize = 24;
annotation('textbox',[0.01 0.95 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.49 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

