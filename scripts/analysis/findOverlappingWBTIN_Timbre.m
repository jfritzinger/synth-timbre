%% findOverlappingWBTIN_Timbre
clear 

%% Load in spreadsheets 

% Load in timbre 
[base, ~, ~, ~] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
timbre = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(timbre, 1);

% Load in WB-TIN
filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/WB-TIN/data/manuscript';
name = 'Data_Table';
WBTIN = readtable(fullfile(filepath, name), 'PreserveVariableNames',true);
num_WBTIN = size(WBTIN, 1);

% Create empty table of overlaps 
table_size = [20 3];
varNames = ["Session", "Timbre", "WBTIN"];
varTypes = ["string", "string", "string"];
overlap = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames);

%% Find overlapping sessions 

i = 1;
for ii = 1:num_data 

	sessions = unique(timbre{ii, 10:33});
	if isempty(sessions{1})
		sessions = sessions(2:end);
	end
	num_sessions = length(sessions);

	for jj = 1:num_sessions 
		target = sessions(jj);

		for kk = 1:num_WBTIN
			WB_sessions = unique(WBTIN{kk, 11:33});
			if isempty(WB_sessions{1})
				WB_sessions = WB_sessions(2:end);
			end

			yes = ismember(target,WB_sessions);
			if yes == 1
				overlap.Session(i) = target{1};
				overlap.Timbre(i) = timbre.Putative_Units(ii);
				overlap.WBTIN(i) = WBTIN.Putative_Units(kk);

				i = i+1;
			end
		end
	end
end

%% Save spreadsheet 


