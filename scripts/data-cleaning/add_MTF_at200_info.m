%% add_MTF_at200_info.m
clear

%% Load in table

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% Loop through each neuron and calculate the MTF at 200Hz 

for ineuron = 1:num_data-1 

	% Load session
	putative = sessions.Putative_Units{ineuron};
	CF = sessions.CF(ineuron);
	load(fullfile(datapath, [putative '.mat']))

	% Load in MTF session
	param = data{3, 2};
	if isempty(param)
		param = data{3, 1};
	end

	% Analyze MTF 
	data_MTF = analyzeMTF(param);

	% Plot 
	cluster = param.cluster;
	stims = param.stims;
	[~, ~, ~, ~, data_MTF2] = plotPhysMTF([], param, []);
	disp(data_MTF2.at_200)

	% Add 'at 200' to spreadsheet
	sessions.MTF_at200{ineuron} = data_MTF.at_200;

end


%% Save out spreadsheet 

writetable(sessions, 'test.xlsx')