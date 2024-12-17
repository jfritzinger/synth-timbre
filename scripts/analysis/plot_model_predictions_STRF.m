%% plot_STRF_model_predictions 
% 
%
% J. Fritzinger
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

binmodes = {'Contra', 'Binaural'};
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
			putative_neuron = sessions.Putative_Units{s_ind};
			CF = sessions.CF(s_ind);
			load(fullfile(datapath, 'neural_data', [putative_neuron '.mat']), 'data');
			params_STRF = data{4,ibin};
			param_ST = data(6+ispl,ibin); 

			% Load in model
			if ispl == 1 && ibin == 1
				filename = [putative_neuron '_STRF_63_Contra'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel_con")
			elseif ispl == 1 && ibin == 2
				filename = [putative_neuron '_STRF_63_Bin'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel")
			elseif ispl == 2 && ibin == 1
				filename = [putative_neuron '_STRF_73_Contra'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel_con")
			elseif ispl == 2 && ibin == 2
				filename = [putative_neuron '_STRF_73_Bin'];
				load(fullfile(savepath, [filename '.mat']), "STRFmodel")
			end

			% General analysis
			% data_STRF = analyzeSTRF(params_STRF);
			data_ST = analyzeST(param_ST);
			param_ST = param_ST{1};
			data_ST = data_ST{1};

			% %% Plot 
			% figure('Position',[3,632,623,273])
			% tiledlayout(1, 2)
			% 
			% % Plot STRF
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
			% 
			% % Plot real response
			% nexttile
			% hold on
			% errorbar(data_ST.fpeaks,data_ST.rate,data_ST.rate_std/(sqrt(param_ST.nrep)), 'LineWidth',1.5);
			% plot(data_ST.fpeaks,(STRFmodel.avModel.*STRFmodel.ratio), 'LineWidth',1.5);
			% xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			% xlim([param_ST.fpeaks(1) param_ST.fpeaks(end)]);
			% grid on
			% xlabel('Tone Frequency (Hz)')
			% ylabel('Avg Rate (sp/s)')
			% ylim([0 STRFmodel.max_all+7])
			% xticklabels(xticks/1000)
			% title(sprintf('STRF Prediction, R^2 = %0.2f\n', STRFmodel.R2))
			% legend({'Data', 'STRF Model'}, 'Location','best')

			%% Create matrix of values 
			R2(ispl, isesh) = STRFmodel.R2;

		end
	end
end

%% Plot results 

R2(R2==0) = NaN;
figure('Position',[932,467,648,246])
tiledlayout(1, 2, 'TileSpacing','compact')

SPLs = [63, 73];
for ispl = 1:2
	nexttile

	edges = linspace(0, 1, 20);
	histogram(R2(ispl,:), edges)
	xlim([0 1])
	

	% Percentage of R^2 values > 0.5
	num_above = sum(R2(ispl,:)>0.5);
	num_all = sum(R2(ispl,:)>0);

	title(sprintf('%d dB SPL, %d percent > 0.5', SPLs(ispl), round(num_above/num_all*100)))
end
