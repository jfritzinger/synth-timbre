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