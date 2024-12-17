%% Population Analysis
% J. Fritzinger, updated 1/9/23
clear

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find number of neurons in each category 
% Find sessions for target MTF type
%MTF_target = 'BE';
%isMTF = strcmp(sessions.MTF, MTF_target);

% Load in spreadsheet with peak information
spreadsheet_name = 'peak_picking.xlsx';
table = readtable(fullfile(datapath, spreadsheet_name));

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
bin200_MTF = bin200; % & isMTF;

%% Overlayed 

figure('Position',[27,634,1053,204])
tiledlayout(1, 4, 'TileSpacing','compact')
spls = [43, 63, 73, 83];
for ispl = 1:4

	indices = find(bin200_MTF(:,ispl));
	num_index = length(indices);

	red = [187, 249, 186]/255;
	pink = [23, 64, 86]/255;
	color_gradient = [linspace(red(1),pink(1),num_index)', linspace(red(2),pink(2),num_index)',...
		linspace(red(3),pink(3),num_index)'];
	CFs = sessions.CF(indices);
	[~, order] = sort(CFs);

	for isesh = 1:num_index

		% Load in session
		putative = sessions.Putative_Units{indices(isesh)};
		CF = sessions.CF(indices(isesh));
		load(fullfile(datapath, 'neural_data', [putative '.mat']))
		params_ST = data(5+ispl, 2);
		params_RM = data{2, 2};
		data_RM = analyzeRM(params_RM);

		peak_ind = find(strcmp(putative, table.Putative) & spls(ispl)==table.SPL);
		if ~isempty(peak_ind)
		type = table.Type{peak_ind};
		if strcmp(type, 'Di')||strcmp(type, 'Flat')
			continue
		end
		end
		
		
		% Analysis
		data_ST = analyzeST(params_ST, CF);
		data_ST = data_ST{1};
		rates = data_ST.rates_sm;
		fpeaks = data_ST.fpeaks;

		% Get peak
		[~, peak_ind] = max(rates);
		peak_freq = fpeaks(peak_ind);
		

		% Align
		fpeaks_re_CF = fpeaks-CF;
		%fpeaks_re_max = fpeaks-peak_freq;
		log_fpeaks = log2(fpeaks/CF);
		log_peak_freq = log_fpeaks(peak_ind);
		log_fpeaks_re_max = log_fpeaks-log_peak_freq;

		% Record difference between the peak frequency and the CF for
		% plotting
		diff_CF(ispl, isesh) = log_peak_freq;

		% Normalize
		rates = zscore(rates);
		%rates = (rates - data_RM.spont)./ std(rates); % z-score w.r.t. spont
		%rates = rates/max(rates);

		% Plot
		nexttile(ispl)
		hold on
		%plot(fpeaks_re_CF, rates, 'linewidth', 1.5, 'Color',color_gradient(order(isesh), :));
		patch([log_fpeaks',NaN],[rates',NaN],'w','EdgeColor','k','LineWidth',1.5,'EdgeAlpha',0.2);
	end
end


spls = [43, 63, 73, 83];
for ispl = 1:4
	% Plotting Params
	nexttile(ispl)
	xlabel('Octaves w.r.t. CF (oct)')
	xline(0, 'k')
	yline(0, 'k')
	box on
	%title([MTF_target ', ' num2str(spls(ispl)) ' dB SPL, F0=200, Bin'])
	xlim([-2.5 2.5])
	ylim([-2 2.5])
	if ispl == 1
		ylabel('Z-score')
	else
		yticklabels([])
	end
	set(gca, 'fontsize', 14)
end

%%
figure
figure('Position',[184,632,1497,287])
tiledlayout(1, 4, 'TileSpacing','compact')
edges = linspace(-3, 2, 31);
for ispl = 1:4
	nexttile
	histogram(diff_CF(ispl, :), edges)
	xline(median(diff_CF(ispl,:)), 'r', 'LineWidth',2)
	xlabel('Peak Freq - CF (Octaves w.r.t. CF)')
	title([MTF_target ', ' num2str(spls(ispl)) ' dB SPL, F0=200, Bin'])
	set(gca, 'fontsize', 14)
end

figure
boxplot(diff_CF')
ylabel('Peak Freq - CF (Octaves w.r.t. CF)')
xlabel('Levels (dB SPL)')
xticklabels({'43', '63', '73', '83'});
title('Difference between peak rate frequency and CF')
set(gca, 'fontsize', 14)