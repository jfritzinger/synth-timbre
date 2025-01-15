%% methods_peak_quantification 

%% Set up figure 

figure('Position',[103,736,671,420])
fontsize = 14;
%tiledlayout(2, 3)

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
	scatter(peaks.locs./1000, peaks.pks, 'filled', 'r')
	num_peaks = length(peaks.pks);
	for ip = 1:num_peaks
		plot([peaks.locs(ip) peaks.locs(ip)]./1000, [peaks.pks(ip)-peaks.p(ip) peaks.pks(ip)], 'Color',"#D95319")
	end
	scatter(dips.locs./1000, -1*dips.pks, 'filled', 'b')
	num_dips = length(dips.pks);
	for ip = 1:num_dips
		plot([dips.locs(ip) dips.locs(ip)]./1000, -1*([dips.pks(ip)-dips.p(ip) dips.pks(ip)]), 'Color',"#7E2F8E")
	end
	%plot(bounds_freq, [halfheight halfheight], 'g')

	% % Annotate prominence & width
	% message = sprintf('Prom: %0.2f', prom);
	% text(0.05, 0.98, message, 'Units', 'normalized', ...
	% 	'VerticalAlignment', 'top', 'FontSize',fontsize)
	% message = sprintf('Width: %0.0f Hz', width);
	% text(0.05, 0.94, message, 'Units', 'normalized', ...
	% 	'VerticalAlignment', 'top', 'FontSize',fontsize)

	plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)]./1000;
	xline(CF./1000, '--', 'Color',CF_color, 'linewidth', 1.5)
	xlim(plot_range)
	xlabel('Spectral Peak Freq. (kHz)')
	if ineuron == 1
		ylabel('Z-score')
		title('Peak')
	elseif ineuron == 2
		yticklabels([])
		title('Dip')
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

	% Plot
	h(3+ineuron) = subplot(2, 3, 3+ineuron);
	hold on
	rates_sm = smooth_rates(rate, rlb, rub, CF);
	errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
		'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{1})
	plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{1})
	plot(fpeaks./1000, rates_sm, 'linewidth', linewidth, 'color', 'k')
	xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)

	% Figure parameters 
	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
	set(gca, 'Fontsize', fontsize, 'XTick', plot_range(1)+0.200:0.400:plot_range(2)-0.200);
	xlim(plot_range);
	grid on
	ylim([0 max_rate+5])
	ylabel('Avg. Rate (sp/s)')
	xlabel('Spectral Peak Freq. (kHz)')

	% if ineuron == 2
	% 	legend({'', 'Data', 'CF', 'Spont.'}, 'Location',...
		% 	'northeastoutside')
	% end

end

%% Arrange and annotate 

left = repmat(linspace(0.1, 0.7, 3), 1, 2);
left(5) = 0.54;
bottom = [0.6 0.6 0.6 0.13 0.13];
width = 0.27;
height = 0.32;

for ii = 1:5
	set(h(ii), 'position', [left(ii) bottom(ii) width height])
end

% Add labels 
labelsize = 24;
annotation('textbox',[0.01 0.95 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.49 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

