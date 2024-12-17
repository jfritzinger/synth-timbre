%% plot_STRF_and_RM_R2
%
%
%
clear

%% Load in data

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);
savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Synth-Timbre/data/manuscript/STRF_Models';
%savepath = 'C:\Users\jfritzinger\Box\02 - Code\Synth-Timbre\data\manuscript\STRF_Models';


%% Analyze and plot data 

binmodes = {'Contra', 'Bin'};
SPLs = {'63', '73'};
for ispl = 1:2
	for ibin = 2

		if ispl == 1 && ibin == 1
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con) & cellfun(@(s) contains(s, 'R'), sessions.STRF_con);
		elseif ispl == 1 && ibin == 2
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB) & cellfun(@(s) contains(s, 'R'), sessions.STRF);
		elseif ispl == 2 && ibin == 1
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con) & cellfun(@(s) contains(s, 'R'), sessions.STRF_con);
		elseif ispl == 2 && ibin == 2
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB) & cellfun(@(s) contains(s, 'R'), sessions.STRF);
		end
		index = find(has_data);

		for isesh = 1:length(index)

			% Load in data
			s_ind = index(isesh);
			putative = sessions.Putative_Units{s_ind};
			CF = sessions.CF(s_ind);
			load(fullfile(datapath, 'neural_data', [putative '.mat']), 'data');
			params_STRF = data{4,ibin};
			param_ST = data(6+ispl,ibin); 
			param_RM = data{2, 2};

			% Load in STRF model
			if ispl == 1 && ibin == 1
				filename = [putative '_STRF_63_Contra'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel_con")
			elseif ispl == 1 && ibin == 2
				filename = [putative '_STRF_63_Bin'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel")
			elseif ispl == 2 && ibin == 1
				filename = [putative '_STRF_73_Contra'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel_con")
			elseif ispl == 2 && ibin == 2
				filename = [putative '_STRF_73_Bin'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel")
			end

			% Load in RM analysis 
			filename = sprintf('%s_RM_%s_%s', putative, SPLs{ispl}, binmodes{ibin});
			load(fullfile(datapath, 'RM_predictions', [filename '.mat']), "RM_R2")

			% General analysis
			% data_STRF = analyzeSTRF(params_STRF);
			data_ST = analyzeST(param_ST, CF);
			param_ST = param_ST{1};
			data_ST = data_ST{1};

			%% Plot 
			
			%tiledlayout(1, 2)

			% Plot STRF
			% nexttile
			% STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
			% imagesc(data_STRF.t, data_STRF.f./1000, STRF_mat, data_STRF.clims_strf);
			% set(gca,'Ydir','normal','XLim',data_STRF.tlims,'YLim',[param_ST.fpeaks(2) param_ST.fpeaks(end)]./1000)
			% hold on
			% yline(CF/1000, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			% colormap(redblue)
			% grid on
			% xlabel('Time (s)');
			% ylabel('Frequency (kHz)')

			% Plot real response
			%nexttile
			% figure('Position',[3,632,623,273])
			% hold on
			% errorbar(data_ST.fpeaks,data_ST.rate,data_ST.rate_std/(sqrt(param_ST.nrep)), 'LineWidth',1.5);
			% plot(data_ST.fpeaks,(STRFmodel.avModel.*STRFmodel.ratio), 'LineWidth',1.5);
			% plot(data_ST.fpeaks, RM_R2.rate_RM, 'LineWidth',1.5);
			% xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			% xlim([param_ST.fpeaks(1) param_ST.fpeaks(end)]);
			% grid on
			% xlabel('Tone Frequency (Hz)')
			% ylabel('Avg Rate (sp/s)')
			% ylim([0 STRFmodel.max_all+7])
			% xticklabels(xticks/1000)
			% title('STRF and RM Predictions')
			% hLeg = legend({'Data', sprintf('STRF, R^2=%0.2f', STRFmodel.R2),...
			% 	sprintf('RM, R^2=%0.2f',RM_R2.R2)}, 'Location','eastoutside');
			% hLegend.ItemTokenSize = [6,6];


			%% Create matrix of values 
			R2_STRF(ispl, isesh) = STRFmodel.R2;
			R2_RM(ispl, isesh) = RM_R2.R2;

		end
	end
end

%% Plot results 

figure('Position',[560,602,613,246])
tiledlayout(1, 2, 'TileSpacing','compact')
SPLs = [63, 73];
for ispl = 1:2
	nexttile

	edges = linspace(0, 1, 20);
	histogram(R2_STRF(ispl,:), edges)
	xlim([0 1])
	

	% Percentage of R^2 values > 0.5
	num_above = sum(R2_STRF(ispl,:)>0.5);
	num_all = sum(R2_STRF(ispl,:)>0);

	title(sprintf('%d dB SPL, %d percent > 0.5', SPLs(ispl), round(num_above/num_all*100)))
end



%% Plot 

figure('Position',[560,524,841,324])
tiledlayout(1, 2)
R2_RM(R2_RM==0) = NaN;
R2_STRF(R2_STRF==0) = NaN;
SPLs = {'63', '73'};
for ispl = 1

	nexttile
	scatter(R2_RM(ispl, :), R2_STRF(ispl, :), 'filled', 'MarkerEdgeColor','k');
	xlabel('R^2 RM')
	ylabel('R^2 STRF')
	title(sprintf('%s dB SPL', SPLs{ispl}))

	RM_temp = R2_RM(ispl, :);
	RM_temp(isnan(RM_temp)) = [];
	STRF_temp = R2_STRF(ispl, :);
	STRF_temp(isnan(STRF_temp)) = [];

	% Add regression line 
	[p, S] = polyfit(RM_temp,STRF_temp,1);
	b = p(1); a = p(2);
	hold on
	x1 = linspace(0,1, 100);
	set(gca, 'fontsize', 16)
	%f1 = polyval(p,x1);
	%plot(x1, f1, 'k')
	%ßlegend({'Data Points',sprintf('Linear Fit  y = %.2f + %.2f*x',a,b)}, 'Location','best')

end


