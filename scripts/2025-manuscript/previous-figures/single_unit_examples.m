%% Single-Unit Examples of ST responses
% J. Fritzinger
%
% Three example units for BE neurons and three example units for BS neurons
clear



%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Set up figure
linewidth = 1.5;
figure('Position',[4,427,1313,480])
%tiledlayout(3, 6, 'TileIndexing','columnmajor')

%% Plot
loc = [1, 7, 13; 2, 8, 14; 3, 9, 15; 4, 10, 16; 5, 11, 17;6, 12, 18]';
%loc = [13, 7, 1; 14, 8, 2; 15, 9, 3; 16, 10, 4; 17, 11, 5; 18, 12, 6]';
j = 1;
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
	data_ST  = analyzeST(params);

	% RM to get spont
	params_RM = data{2, 2};
	data_RM = analyzeRM(params_RM);
	spont = data_RM.spont;

	num_ds = size(data_ST, 2);
	if num_ds == 4
		%data_colors = {'#1b9e77', '#d95f02', '#7570b3', '#e7298a'};
		%data_colors = {'#ffffcc', '#a1dab4', '#41b6c4', '#225ea8'};
		% Flipped
		data_colors = {'#034E1C', '#03882F', '#3F985C', '#82BB95'};
	else
		%data_colors = {'#1b9e77', '#d95f02', '#e7298a'};
		%data_colors = {'#a1dab4', '#41b6c4', '#225ea8'};
		% Flipped
		data_colors = {'#034E1C', '#03882F', '#82BB95'};
	end

	% Sort
	spls = cell2mat(cellfun(@(p) p.spl, data_ST, 'UniformOutput',false));
	[~, order] = sort(spls);
	order = fliplr(order);
	max_rate = max(cellfun(@(d) max(d.rate), data_ST));

	label_ind = 1;
	for ind = 1:3 %num_ds
		h(loc(ind, ineuron)) = subplot(3, 6, loc(ind, ineuron));
		j = j+1;
		hold on

		rate = data_ST{order(ind)}.rate;
		rate_std = data_ST{order(ind)}.rate_std;
		rlb = data_ST{order(ind)}.rlb;
		rub = data_ST{order(ind)}.rub;
		fpeaks = data_ST{order(ind)}.fpeaks;
		spl = data_ST{order(ind)}.spl;
		rate_sm = data_ST{order(ind)}.rates_sm;

		% Plot
		rates_sm = smooth_rates(rate, rlb, rub);
		errorbar(fpeaks./1000, rate, rate_std/sqrt(params{order(ind)}.nrep), 'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{ind})
		plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{:,ind})
		%plot(fpeaks, rate_sm, 'linewidth', linewidth, 'color', data_colors{ind})
		label{1} = [num2str(params{order(ind)}.spl) ' dB SPL'];
		%label_ind = label_ind+1;


		plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
		xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
		%label{2} = 'Estimated CF';
		yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)
		yline(0.1, 'k', LineWidth=linewidth)
		%label(label_ind+1) = {'Spont'};

		set(gca, 'Fontsize', 14, 'XTick', plot_range(1)+0.200:0.400:plot_range(2)-0.200);
		xlim(plot_range);
		grid on
		ylim([0 max_rate+5])
		if ind == 1 && ineuron==2
			xticklabels([])
			title('BS', 'FontSize',20)
		elseif ind == 1 && ineuron==5
			title('BE', 'FontSize',20)
			xticklabels([])
		elseif ind == 1
			xticklabels([])
		elseif ind == 2 && (ineuron==1 || ineuron ==4)
			ylabel('Avg. rate (sp/s)')
			xticklabels([])
		elseif ind == 2
			xticklabels([])
		elseif ind == 3 && (ineuron==2 || ineuron ==5)
			xlabel('Spectral Peak  Frequency (Hz)')
		end
		if ineuron == 1 || ineuron == 4
			legend(label, 'Box','off', 'Location','southwest')
		end	
	end
end

%% Move 

%left = linspace(0.05, 0.85, 6);
left = [0.07 0.22 0.37 0.55 0.70 0.85]+0.01;
bottom = fliplr(linspace(0.1, 0.64, 3));
width = 0.13;
height = 0.27;

col = repmat(left, 1, 3);
row = reshape(repmat(bottom, 6, 1), 18, 1);

for ii = 1:18
	set(h(ii), 'position', [col(ii) row(ii) width height])
end

annotation('textbox',...
	[0.00571210967250573 0.725041666666666 0.0455064737242955 0.151041666666667],...
	'String',{'83','dB','SPL'},...
	'HorizontalAlignment','center',...
	'FontSize',18,...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.00952018278750954 0.166708333333333 0.0455064737242955 0.151041666666667],...
	'String',{'43','dB','SPL'},...
	'HorizontalAlignment','center',...
	'FontSize',18,...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.00647372429550647 0.441708333333333 0.0455064737242953 0.151041666666667],...
	'String',{'63','dB','SPL'},...
	'HorizontalAlignment','center',...
	'FontSize',18,...
	'FitBoxToText','off',...
	'EdgeColor','none');

%% Export figure

%exportgraphics(gcf, fullfile(savepath, 'manuscript', 'examples-differentMTFs.png'), 'Resolution', 600)
