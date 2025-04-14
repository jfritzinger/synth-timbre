%% plot_VS_analysis
clear 

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name),...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Load in spreadsheet with peak information
spreadsheet_name = 'peak_picking_VS.xlsx';
table = readtable(fullfile(datapath, spreadsheet_name));

%% Histogram overall 
fontsize = 12;

% Plot histogram  
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = table.binmode == 2;
for ispl = 1:4
	isSPL = table.SPL == spl(ispl);


	index = isSPL  & isBin;

	num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), table.Type(index)));
	num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), table.Type(index)));
	num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), table.Type(index)));
	all = sum([num_peak num_dip num_flat]);

	percent_peak(ispl) = num_peak;
	percent_dip(ispl) = num_dip;
	percent_flat(ispl) = num_flat;
	percent_all = [percent_peak; percent_dip; percent_flat]';

	bar(percent_all)
	%title([num2str(spl(ispl)) ' dB SPL'])
	xticklabels(spl)
	legend('Peak', 'Dip', 'Flat', 'Location','south')
	ylabel('# Neurons')
	xlabel('Level (dB SPL)')
	ylim([0 100])
	set(gca, 'fontsize', fontsize)
	title('Peak/Dip Picking in VS Profiles')
end

%%
figure

% Plot histogram  
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = table.binmode == 2;
for ispl = 2
	isSPL = table.SPL == spl(ispl);

	for iMTF = 1:4

		if iMTF == 1
			MTF_target = 'BS';
			isMTF = strcmp(table.MTF, MTF_target);
		elseif iMTF == 2
			MTF_target = 'BE';
			isMTF = strcmp(table.MTF, MTF_target);
		elseif iMTF == 3
			MTF_target = 'Hybrid';
			isMTF = contains(table.MTF, 'H');
		else
			MTF_target = 'F';
			isMTF = strcmp(table.MTF, MTF_target);
		end
		index = isSPL & isMTF & isBin;

		num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), table.Type(index)));
		num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), table.Type(index)));
		num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), table.Type(index)));
		all = sum([num_peak num_dip num_flat]);

		percent_peak(iMTF) = num_peak/all*100;
		percent_dip(iMTF) = num_dip/all*100;
		percent_flat(iMTF) = num_flat/all*100;
	end
	percent_all = [percent_peak; percent_dip; percent_flat]';

	bar(percent_all, 'stacked')
	%title([num2str(spl(ispl)) ' dB SPL'])
	xticklabels({'BS', 'BE', 'Hybrid', 'Flat'})
	legend('Peak', 'Dip', 'Flat', 'Location','south')
	ylabel('% Neurons')
	xlabel('MTF Type')
	ylim([0 100])
	set(gca, 'fontsize', fontsize)

end

%% Average vector strength 

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name),...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Create has_data
data_ind(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
data_ind(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
data_ind(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
data_ind(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);

for ispl = 1:4
	has_data = data_ind(:,ispl);
	data_index = find(has_data);
	num_neurons = sum(has_data);

	% Plot each neuron
	ii = 1;
	for isesh = 1:num_neurons
		ineuron = data_index(isesh);

		% Load in session
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, 'neural_data', [putative '.mat']))

		if data_ind(ineuron,ispl)==1

			% Load in proper dataset for each idata
			param_ST = data(5+ispl, 2);

			% Analyze synthetic timbre
			data_ST = analyzeST(param_ST, CF);
			data_ST = data_ST{1};
			[~, peak_ind] = max(data_ST.rate);
			peak_f = data_ST.fpeaks(peak_ind);
			[~, peaksm_ind] = max(data_ST.rates_sm);
			peak_fsm = data_ST.fpeaks(peaksm_ind);

			temporal = analyzeST_Temporal(param_ST{1}, data_ST);
			VS(isesh, ispl) = max(temporal.VS);

			% Get ind nearest to CF 
			[val, ind] = min(abs(param_ST{1}.fpeaks-CF));
			VS2(isesh, ispl) = temporal.VS(ind);

		end
		fprintf('%s done, %d percent done\n', putative, round(isesh/num_neurons*100))
	end
end

fprintf('43 dB SPL: Mean VS = %0.02f\n', mean(VS(:,1)))
fprintf('63 dB SPL: Mean VS = %0.02f\n', mean(VS(:,2)))
fprintf('73 dB SPL: Mean VS = %0.02f\n', mean(VS(:,3)))
fprintf('83 dB SPL: Mean VS = %0.02f\n', mean(VS(:,4)))

fprintf('43 dB SPL: Mean VS = %0.02f\n', mean(VS2(:,1)))
fprintf('63 dB SPL: Mean VS = %0.02f\n', mean(VS2(:,2)))
fprintf('73 dB SPL: Mean VS = %0.02f\n', mean(VS2(:,3)))
fprintf('83 dB SPL: Mean VS = %0.02f\n', mean(VS2(:,4)))