%% plot_RM_predictions.m
clear 

%% Load in spreadsheet

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% 

bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
contra200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_con);
contra200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con);
contra200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con);
contra200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_con);
has_RM = cellfun(@(s) contains(s, 'R'), sessions{:,13});

%% Run model for all units 

binmodes = {'Contra', 'Bin'};
SPLs = {'43', '63', '73', '83'};
for ispl = 1:4
	for ibin = 2

		% Find sessions of interest
		if ibin == 1
			has_ST = contra200(:,ispl);
		else
			has_ST = bin200(:,ispl);
		end
		has_data = has_RM & has_ST;
		index = find(has_data);

		for isesh = 1:length(index)

			% Load in data
			s_ind = index(isesh);
			putative = sessions.Putative_Units{s_ind};
			CF = sessions.CF(s_ind);
			load(fullfile(datapath,'neural_data' , [putative '.mat']), 'data');
			param_ST = data(5+ispl,ibin); 
			data_ST = analyzeST(param_ST, CF);
			data_ST = data_ST{1};
			param_ST = param_ST{1};
			
			% General analysis
			params_RM = data{2,2};
			data_RM = analyzeRM(params_RM);

			% Load in RM R2
			filename = sprintf('%s_RM_%s_%s', putative, SPLs{ispl}, binmodes{ibin});
			load(fullfile(datapath, 'RM_predictions', [filename '.mat']), "RM_R2")

			% Add to matrix
			R2(ispl, isesh) = RM_R2.R2;
			CFs(ispl, isesh) = CF;
			x = 1;

			%% Plot 
			%figure('Position',[3,632,623,273])
			% % Plot real response
			% nexttile
			% hold on
			% errorbar(data_ST.fpeaks,data_ST.rate,data_ST.rate_std/(sqrt(param_ST.nrep)), 'LineWidth',1.5);
			% %plot(data_ST.fpeaks,(avModel.*ratio), 'LineWidth',1.5);
			% plot(data_ST.fpeaks, rate_RM, 'LineWidth',1.5)
			% xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			% xlim([param_ST.fpeaks(1) param_ST.fpeaks(end)]);
			% grid on
			% xlabel('Tone Frequency (Hz)')
			% ylabel('Avg Rate (sp/s)')
			% xticklabels(xticks/1000)
			% title(sprintf('R^2 = %0.2f\n', R2))
		
		end
	end
end

%% Plot results 

R2(R2==0) = NaN;
figure('Position',[560,602,983,246])
tiledlayout(1, 4, 'TileSpacing','compact')
SPLs = [43, 63, 73, 83];
for ispl = 1:4
	nexttile

	edges = linspace(0, 1, 20);
	histogram(R2(ispl,:), edges)
	hold on
	xlim([0 1])
	xline(median(R2(ispl, :), 'omitnan'), '--r', 'LineWidth',2)
	xline(mean(R2(ispl, :), 'omitnan'), 'r', 'LineWidth',2)

	% Percentage of R^2 values > 0.5
	num_above = sum(R2(ispl,:)>0.5);
	num_all = sum(R2(ispl,:)>0);

	title(sprintf('%d dB SPL, %d percent > 0.5', SPLs(ispl), round(num_above/num_all*100)))
end

%% 

R2(R2==0) = NaN;
figure()
SPLs = [43, 63, 73, 83];

boxplot(R2')
hold on
num = size(R2,2);
x = [ones(1,num); 2*ones(1,num); 3*ones(1,num); 4*ones(1,num)];
swarmchart(x', R2', 'filled')
ylim([0 1])
xticklabels(SPLs)
xlabel('Level (dB SPL)')
ylabel('R^2')
title('RM R^2 Values')
set(gca, 'Fontsize', 16)
box off 

%% 

R2(R2==0) = NaN;
figure('Position',[560,602,983,246])
tiledlayout(1, 4, 'TileSpacing','compact')

SPLs = [43, 63, 73, 83];
for ispl = 1:4
	nexttile

	edges = linspace(0, 1, 20);
	scatter(CFs(ispl, :), R2(ispl,:), 'filled')	

	% Percentage of R^2 values > 0.5
	num_above = sum(R2(ispl,:)>0.5);
	num_all = sum(R2(ispl,:)>0);

	title(sprintf('%d dB SPL, %d percent > 0.5', SPLs(ispl), round(num_above/num_all*100)))
end








