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
MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
bin200_MTF = bin200 & isMTF;

%% Scatter plot of CF vs HHBW (all MTF)

has_data = bin200_MTF(:,1) | bin200_MTF(:,2) | bin200_MTF(:,3) | bin200_MTF(:,4);
indices = find(has_data);
num_index = length(indices);

width_CF = NaN(num_index, 4);
width = NaN(num_index, 4);
%CFs = sessions.CF(1:334);
CFs = sessions.CF(indices);
for isesh = 1:num_index

	% Load in session
	putative = sessions.Putative_Units{indices(isesh)};
	CF = sessions.CF(indices(isesh));
	load(fullfile(datapath, [putative '.mat']))

	for ispl = 1:4

		% Analysis
		params_ST = data(5+ispl, 2);
		if ~isempty(params_ST{1})
			data_ST = analyzeST(params_ST);
			data_ST = data_ST{1};

			% Get CF rate
			[~, CF_ind] = min(abs(CF-data_ST.fpeaks));
			CF_rate = data_ST.rates_sm(CF_ind);

			% Get width
			width(isesh, ispl) = data_ST.width;
			if data_ST.width ~= 2400
				width_CF(isesh, ispl) = data_ST.width/CF;
			end
		end
	end
end

% Plot
figure('position', [519,533,1150,305])
tiledlayout(1, 4, 'TileSpacing','tight', 'Padding','compact')
spls = [43, 63, 73, 83];
for ispl = 1:4
	
	nexttile(ispl)
	scatter(CFs/1000, width_CF(:,ispl), 'filled')
	box on
	number = width_CF(:,ispl);
	number(isnan(number)) = [];
	title([num2str(spls(ispl)) ' dB SPL, n=' num2str(length(number))])
	xlabel('CF (kHz)')
	if ispl == 1
		ylabel('1/2 Height BW / CF')
	end
	ylim([0.02 5])
	xlim([0.3 10])
	set(gca, 'fontsize', 17)
	hold on
	set(gca, 'XScale', 'log')
	set(gca, 'YScale', 'log')
	xticks([0 200 500 1000 2000 5000 10000]/1000)
	yticks([0 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000])
	yline(1, '--k', 'linewidth', 2)
	grid on

end
