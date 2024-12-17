%% save_RM_R2.m
%
%
% J. Fritzinger, created 10/6/24

%% Load in spreadsheet

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);
%savepath = 'C:\Users\jfritzinger\Box\02 - Code\Synth-Timbre\data\manuscript\RM_predictions';
savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Synth-Timbre/data/manuscript/RM_predictions';

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
			data_ST = analyzeST(param_ST);
			data_ST = data_ST{1};
			param_ST = param_ST{1};
			
			% General analysis
			params_RM = data{2,2};
			data_RM = analyzeRM(params_RM);

			% Calculate RM R2 
			[R2, BF_spl, rate_RM] = calculateRMR2(param_ST, data_ST, data_RM);

			% Display progress
			fprintf('%s Done, %.2f through %s \n', putative, isesh/length(index), binmodes{ibin})

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

			%% Add to struct
			temp(isesh).putative = putative;
			temp(isesh).R2 = R2;
			temp(isesh).BF_spl = BF_spl;
			temp(isesh).rate_RM = rate_RM;

			RM_R2 = temp(isesh);
			filename = sprintf('%s_RM_%s_%s', putative, SPLs{ispl}, binmodes{ibin});
			save(fullfile(savepath, [filename '.mat']), "RM_R2")
		
		end
	end
end
