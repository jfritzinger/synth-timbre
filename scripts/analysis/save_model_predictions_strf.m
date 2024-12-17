%% strf_model_predictions
%
%
% J. Fritzinger, created 10/6/24

%% Load in spreadsheet

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);
%savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Synth-Timbre/data/manuscript/STRF_Models';
savepath = 'C:\Users\jfritzinger\Box\02 - Code\Synth-Timbre\data\manuscript\STRF_Models';

%% Run model for all units 

binmodes = {'Contra', 'Binaural'};
SPLs = {'63', '73'};
for ispl = 1
	for ibin = 1

		% Find sessions of interest
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
			
			% General analysis
			data_STRF = analyzeSTRF(params_STRF);
			data_ST = analyzeST(param_ST);
			param_ST = param_ST{1};
			data_ST = data_ST{1};

			% Calculate model response
			param_ST.Fs = 48000;
			param_ST.mnrep = param_ST.nrep;
			param_ST.dur = 0.3;
			[R2, avModel, stdModel, ratio, max_all] = modelTimbreSTRF(param_ST, data_STRF, data_ST);

			% Display progress
			fprintf('%s Done, %.2f through %s \n', putative_neuron, isesh/length(index), binmodes{ibin})

			%% Plot 
			figure('Position',[3,632,623,273])
			tiledlayout(1, 2)

			% Plot STRF
			nexttile
			STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
			imagesc(data_STRF.t, data_STRF.f./1000, STRF_mat, data_STRF.clims_strf);
			set(gca,'Ydir','normal','XLim',data_STRF.tlims,'YLim',[param_ST.fpeaks(2) param_ST.fpeaks(end)]./1000)
			hold on
			yline(CF/1000, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			colormap(redblue)
			grid on
			xlabel('Time (s)');
			ylabel('Frequency (kHz)')

			% Plot real response
			nexttile
			hold on
			errorbar(data_ST.fpeaks,data_ST.rate,data_ST.rate_std/(sqrt(param_ST.nrep)), 'LineWidth',1.5);
			plot(data_ST.fpeaks,(avModel.*ratio), 'LineWidth',1.5);
			xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			xlim([param_ST.fpeaks(1) param_ST.fpeaks(end)]);
			grid on
			xlabel('Tone Frequency (Hz)')
			ylabel('Avg Rate (sp/s)')
			ylim([0 max_all+7])
			xticklabels(xticks/1000)
			title(sprintf('R^2 = %0.2f\n', R2))

			%% Add to struct
			temp(isesh).putative = putative_neuron;
			temp(isesh).R2 = R2;
			temp(isesh).avModel = avModel;
			temp(isesh).stdModel = stdModel;
			temp(isesh).ratio = ratio;
			temp(isesh).max_all = max_all;
			temp(isesh).SPL = param_ST.spl;

			if ispl == 1 && ibin == 1
				STRFmodel_con = temp(isesh);
				filename = [putative_neuron '_STRF_63_Contra'];
				save(fullfile(savepath, [filename '.mat']), "STRFmodel_con")
			elseif ispl == 1 && ibin == 2
				STRFmodel = temp(isesh);
				filename = [putative_neuron '_STRF_63_Bin'];
				save(fullfile(savepath, [filename '.mat']), "STRFmodel")
			elseif ispl == 2 && ibin == 1
				STRFmodel_con = temp(isesh);
				filename = [putative_neuron '_STRF_73_Contra'];
				save(fullfile(savepath, [filename '.mat']), "STRFmodel_con")
			elseif ispl == 2 && ibin == 2
				STRFmodel = temp(isesh);
				filename = [putative_neuron '_STRF_73_Bin'];
				save(fullfile(savepath, [filename '.mat']), "STRFmodel")
			end
		end


		% %% Save Data
		% % Add to matrix for analysis
		% if ibin == 1
		% 	STRF_contra = temp;
		% 	save(fullfile(path,'Fig6_Data' , 'Fig6_STRF_23_contra.mat'), "STRF_contra")
		% else
		% 	STRF = temp;
		% 	save(fullfile(path,'Fig6_Data' , 'Fig6_STRF_23_New.mat'), "STRF")
		% end
		% clear temp
	end
end


%% FUNCTIONS 
