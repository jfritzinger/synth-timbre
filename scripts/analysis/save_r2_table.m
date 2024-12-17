%% save_r2_table.m
%
% Script to create an excel spreadsheet with all R2 model results for each
% neuron at each level. Only for binaural stimuli. Used for
% 'plot_model_evaluations.m'. 
%
%
% Author: J. Fritzinger
% Created: 2022-09-13; Last revision: 2024-09-26 
%
% -------------------------------------------------------------------------
clear 

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%%  Synthetic Timbre Data Table 
 
% Initialize spreadsheet columns
%modelpath = '/Volumes/Synth-Timbre/data/manuscript';
modelpath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript';
varNames = ["Putative", "CF", "MTF", "BMF"...
	"SPL", "SFIE_R", "SFIE_R2",...
	"Energy_R", "Energy_R2", ...
	"SFIE_Pop_R", "SFIE_Pop_R2", ...
	"Lat_Inh_R", "Lat_Inh_R2"];
varTypes = ["string", "double", "string", "double"...
	"double", "double", "double", ...
	"double", "double", ...
	"double", "double",...
	"double", "double"];
est_num_rows = 429; % set to number larger than
num_cols = length(varNames);
table_size = [est_num_rows num_cols];
tables = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
isMTF = strcmp(sessions.MTF, 'BE')|strcmp(sessions.MTF, 'BS');
bin200_MTF = bin200 & isMTF;

% Add R^2 values to the spreadsheet 
spls = [43, 63, 73, 83];
has_data = bin200_MTF(:,1) | bin200_MTF(:,2) | bin200_MTF(:,3) | bin200_MTF(:,4);
indices = find(has_data);
num_index = length(indices);
ii = 1;
for isesh = 1:num_index

	% Load in data
	putative = sessions.Putative_Units{indices(isesh)};
	load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE')
	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
	load(fullfile(modelpath,'SFIE_pop_model', [putative '_SFIE_pop.mat']), 'SFIE_pop')
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')

	for ispl = 1:4
		if ~isempty(SFIE{ispl})
			tables.Putative{ii} = sessions.Putative_Units{indices(isesh)};
			tables.CF(ii) = sessions.CF(indices(isesh));
			tables.MTF{ii} = SFIE{ispl}.MTF_shape;
			tables.BMF(ii) = SFIE{ispl}.BMF;
			tables.SPL(ii) = spls(ispl);
			tables.SFIE_R(ii) = SFIE{ispl}.R;
			tables.SFIE_R2(ii) = SFIE{ispl}.R2;
			tables.Energy_R(ii) = energy{ispl}.R;
			tables.Energy_R2(ii) = energy{ispl}.R2;
			tables.SFIE_Pop_R(ii) = SFIE_pop{ispl}.R;
			tables.SFIE_Pop_R2(ii) = SFIE_pop{ispl}.R2;
			tables.Lat_Inh_R(ii) = lat_inh{ispl}.R;
			tables.Lat_Inh_R2(ii) = lat_inh{ispl}.R2;
			ii = ii + 1;
		end
	end
	fprintf('%s done, %d percent done\n', putative, round(isesh/num_index*100))
end

% Save table
writetable(tables,'model_r2_values_ST.xlsx')




