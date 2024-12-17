	%% Population Analysis
% J. Fritzinger, updated 1/9/23
clear

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find number of neurons in each category 
% Find sessions for target MTF type
MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);
%isMTF = contains(sessions.MTF, 'H');


% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
%bin200_MTF = bin200 & isMTF;
bin200_MTF = bin200;

%% Plot all data on one axis for each level 

figure('Position',[1707,-8,940,955])
tiledlayout(1, 4,"TileSpacing","compact", 'Padding','compact')
levels = [43, 63, 73, 83];
plots = [1:3; 5:7; 9:11; 13:15];
for ispl = 1:4

	nexttile
	%h(ispl) = subplot(4, 4, plots(ispl, :));
	hold on
	indices = find(bin200_MTF(:,ispl));
	num_index = length(indices);
	width = NaN(num_index,1);

	% Sort by CF 
	CF_list = sessions.CF(indices);
	[~, order] = sort(CF_list);

	for isesh = 1:num_index
		ind = indices(order(isesh));

		% Load in session
		putative = sessions.Putative_Units{ind};
		CF = sessions.CF(ind);
		load(fullfile(datapath, [putative '.mat']))
		params_ST = data(5+ispl, 2);
		
		% Analysis
		data_ST = analyzeST(params_ST);
		data_ST = data_ST{1};

		% Get CF rate
		[~, CF_ind] = min(abs(CF-data_ST.fpeaks));
		CF_rate = data_ST.rate(CF_ind);
		CF_sm = data_ST.rates_sm(CF_ind);
		
		% Z-score of rates 
		% z_rate = zscore(data_ST.rates_sm); % Z-score
		% CF_z_rate = z_rate(CF_ind); % Z-score

		% Normalize peak rate to 1 
		z_rate = data_ST.rates_sm/max(data_ST.rates_sm);
		CF_z_rate = z_rate(CF_ind);

		% Calculate HHBW
		width(isesh) = data_ST.width;

		% Plot
		%plot(h(ispl), data_ST.fpeaks,data_ST.rates_sm, 'linewidth', 1.5);
		%scatter(h(ispl), CF, CF_sm, 'filled', 'MarkerEdgeColor','r', 'MarkerFaceColor','r')
		plot(data_ST.fpeaks,z_rate+(isesh/2),  'k','linewidth', 1.5);
		scatter(CF, CF_z_rate+(isesh/2), 15, 'filled', 'MarkerEdgeColor',...
			'r', 'MarkerFaceColor','r')
	end

	% Plotting Params
	xlabel('Frequency (Hz)')
	ylabel('Avg. Rate (sp/s)')
	box on
	title(sprintf('%s, %d dB SPL, Bin',...
		MTF_target, levels(ispl)))
	set(gca, 'fontsize', 14)
	xlim([0 9000]);
	set(gca, 'XScale', 'log')

	% Width figure
	% h(4+ispl) = subplot(4, 4, 4*ispl);
	% histogram(h(4+ispl), width, 10)
	% xlabel(h(4+ispl), 'Half-Height BW (Hz)')
	% ylabel(h(4+ispl), '# Neurons')
	% set(h(4+ispl), 'fontsize', 14)
end