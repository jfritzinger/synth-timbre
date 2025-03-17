%% changes_over_level
clear 

%% Load in data

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath,"LMM", "peak_picking_excludeflat.xlsx"));
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Set up figure 

figure('Position',[560,214,9*ppi,7*ppi])
fontsize = 12;
linewidth = 1;
legsize = 8;
labelsize = 22;
titlesize = 16; 

%% Set up and plot examples 

plot_ind = [1, 2, 4, 5, 7, 8];
for ineuron = 1:6

	h(plot_ind(ineuron)) = subplot(3, 3, plot_ind(ineuron));

	switch ineuron
		case 1 % Sharpening
			putative = 'R27_TT3_P7_N08';
		case 2 % Sharpening
			putative = 'R27_TT3_P7_N14';
		case 3 % No Change 
			putative = 'R29_TT4_P5_N02';
		case 4 % No Change 
			putative = 'R27_TT2_P8_N02';
			% okay: 'R24_TT2_P13_N02'; 'R24_TT2_P13_N06';
		case 5 % Broadening
			putative = 'R24_TT2_P13_N03';
		case 6 % Broadening
			putative = 'R24_TT2_P13_N05';
	end

	% Load in data
	[base, datapath, savepath, ppi] = getPaths();
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename)), 'data';
	index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
	CF = sessions.CF(index);

	params = data([6,7,9], 2);
	params = params(~cellfun(@isempty, params));
	data_ST  = analyzeST(params, CF);

	% RM to get spont
	params_RM = data{2, 2};
	data_RM = analyzeRM(params_RM);
	spont = data_RM.spont;

	num_ds = size(data_ST, 2);
	data_colors = {'#034E1C', '#03882F', '#82BB95'};

	% Sort
	spls = cell2mat(cellfun(@(p) p.spl, data_ST, 'UniformOutput',false));
	[~, order] = sort(spls);
	order = fliplr(order);
	max_rate = max(cellfun(@(d) max(d.rate), data_ST));
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
		errorbar(fpeaks./1000, rate, rate_std/sqrt(30), 'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{ind})
		%plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{:,ind})
		plot(fpeaks./1000, rates_sm, 'LineWidth',linewidth, 'Color',data_colors{:,ind})
		[peaks, dips, type, prom, width, lim, ~, ~, freq] = peakFinding(data_ST{order(ind)}, CF);
		Q(ind) = freq/width;
	end
	yline(spont, 'k', LineWidth=linewidth)
	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
	xline(CF./1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	xlabel('Spectral Peak Freq. (kHz)')
	ylabel('Avg. rate (sp/s)')
	set(gca, 'Fontsize', fontsize);
	xlim(plot_range);
	grid on

	hLeg = legend('', sprintf('83, Q=%0.1f', Q(1)), '', sprintf('63, Q=%0.1f', Q(2)), ...
		'', sprintf('43, Q=%0.1f', Q(3)));
	hLeg.ItemTokenSize = [6,6];
	hLeg.FontSize = legsize;
	hLeg.Box = 'off';
	hLeg.Location = "northeast";

	% Add Q:
	% tbl = table([43, 63,83]', Q', 'VariableNames', {'X', 'Q'});
	% mdl = fitlm(tbl, 'Q ~ X');
	% slope = mdl.Coefficients.Estimate(2); % slope
	% 
	% plot([43, 63,83], Q)
	% hold on
	% y_fit = predict(mdl, table([43, 63,83]', 'VariableNames', {'X'}));
	% plot([43, 63,83], y_fit, 'r-', 'LineWidth', 2);
	% 
	% criteria = 0.02; % backwards because starting from high level
	% if slope<criteria && slope > -1*criteria
	% 	changing = 'no change';
	% elseif slope<-1*criteria
	% 	changing = 'sharpening';
	% else
	% 	changing = 'broadning';
	% end
	% 
	% label = sprintf('Q43=%0.2f, Q63=%0.2f, Q83=%0.2f, %s', Q(1), ...
	% 	Q(2), Q(3), changing);
	% text(0.05, 0.95, label, 'Units', 'normalized', ...
		% 	'VerticalAlignment', 'top', 'FontSize',16)
end


%% Plot overall

% Find sessions for target synthetic timbre response
all_neurons = tables.Putative;
neurons = unique(all_neurons);
num_units = size(neurons, 1);
isbin = tables.binmode == 2;
is200 = tables.F0 == 200;

SPLs = [43, 63, 73, 83];
qs = NaN(num_units, 4);
for isesh = 1:num_units

	putative = neurons{isesh};
	isput = cellfun(@(s) strcmp(s, putative), tables.Putative);

	for ispl = 1:4
		ind = isput & isbin & is200 & tables.SPL==SPLs(ispl);
		if any(ind)
			qs(isesh, ispl) = tables.Q(ind);
			qs_log(isesh, ispl) = tables.Q_log(ind);
			CF_group(isesh) = tables.CF_Group(ind);
		end
	end
end


% Get matrix of units with 43, 63, 83 dB 
qs2 = qs(:,[1,2,4]);
rows_with_nan = any(isnan(qs2),2);
qs2(rows_with_nan,:) = [];
x = [43, 63, 83];

% Try 1: Criteria using slope
for ii = 1:length(qs2)

	y = qs2(ii, :)';
    tbl = table(x', y, 'VariableNames', {'X', 'Q'});
    mdl = fitlm(tbl, 'Q ~ X');
    slopes(ii) = mdl.Coefficients.Estimate(2); % slope

	% % Plot
	% figure
	% plot(x, y)
	% hold on
	% y_fit = predict(mdl, table(x', 'VariableNames', {'X'}));
	% plot(x, y_fit, 'r-', 'LineWidth', 2);
end
criteria = 0.03;
same = slopes<criteria & slopes > -1*criteria;
decrease = slopes<-1*criteria;
increase = slopes>criteria;
spls = [43, 63, 83];

indices = [3, 6, 9];
for ii = 1:3
	if ii == 1
		values = increase;
		color = [27,158,119]/256;
	elseif ii == 2
		values = same;
		color = [217,95,2]/256;
	else
		values = decrease;
		color = [117,112,179]/256;
	end

	% Increase
	h(indices(ii)) = subplot(3, 3, indices(ii));
	hold on
	plot(spls, qs2(values,:)', 'color',color , 'LineWidth',linewidth)
	xticks(spls)
	ylabel('Q')
	xlim([40 86])
	plot(spls, mean(qs2(values,:), 'omitnan'), 'k', 'LineWidth',2)
	plot(spls, median(qs2(values,:), 'omitnan'), ':k', 'LineWidth',2)
	set(gca, 'fontsize', fontsize)
	xlabel('Level (dB SPL)')

	label = ['n=' num2str(sum(values))];
	text(0.05, 0.95, label, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)
	hLeg = legend;
	num_lines = size(hLeg.String,2);
	for iii = 1:num_lines
		if iii==num_lines
			leg{iii} = 'Mean';
		elseif iii == num_lines-1
			leg{iii} = 'Median';
		else
			leg{iii} = '';
		end
	end
	hLeg = legend(leg, 'FontSize',legsize);
	hLeg.ItemTokenSize = [15,6];
	hLeg.Box = 'off';
	ylim([0 15])

end


%% Arrange figure 

left = [0.13 0.43 0.74];
bottom = linspace(0.07, 0.73, 3);
height = 0.23;
width = 0.23;

left = repmat(left, 1, 3);
bottom = repmat(bottom, 3, 1);
bottom = fliplr(reshape(bottom, 1, 9));

for ii = 1:9
	set(h(ii), 'Position', [left(ii) bottom(ii) width height])
end

% Annotations
titles_y = { 'Broadening','No Change','Sharpening', };
locs = linspace(0.095, 0.75, 3);
for ii = 1:3
	annotation('textbox',[0.06 locs(ii) 0.15 0.0459],...
		'String',titles_y{ii},...
		'FontSize',titlesize,'EdgeColor','none', 'Rotation',90);
end

bottom = linspace(0.32, 0.96, 3);
annotation('textbox',[0.01 bottom(3) 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 bottom(2) 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 bottom(1) 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');


