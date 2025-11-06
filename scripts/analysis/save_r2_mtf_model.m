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
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%%  Synthetic Timbre Data Table

% Initialize spreadsheet columns
modelpath ='/Users/johannafritzinger/Library/CloudStorage/Dropbox-UniversityofMichigan/Johanna Fritzinger/02-projects/synth-timbre/data';

%modelpath = '/Volumes/Synth-Timbre/data/manuscript';
%modelpath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript';
varNames = ["Putative", "CF", "MTF", "BMF"...
	"SPL", "SFIE_R", "SFIE_R2", "Lat_Inh_R", "Lat_Inh_R2"];
varTypes = ["string", "double", "string", "double"...
	"double", "double", "double","double", ...
	"double"];
est_num_rows = 128; % set to number larger than
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
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	for ispl = 2
		if ~isempty(SFIE{ispl})


			% Fill out table
			tables.Putative{ii} = sessions.Putative_Units{indices(isesh)};
			tables.CF(ii) = sessions.CF(indices(isesh));
			tables.MTF{ii} = SFIE{ispl}.MTF_shape;
			tables.BMF(ii) = SFIE{ispl}.BMF;
			tables.SPL(ii) = spls(ispl);
			tables.SFIE_R(ii) = SFIE{ispl}.R_MTF;
			tables.SFIE_R2(ii) = SFIE{ispl}.R_MTF^2;
			tables.SFIE_MTFshape{ii} = SFIE{ispl}.MTF_shape_calc;
			tables.Lat_Inh_R(ii) = lat_inh{ispl}.R_MTF;
			tables.Lat_Inh_R2(ii) = lat_inh{ispl}.R_MTF^2;
			tables.Lat_Inh_MTFshape{ii} = lat_inh{ispl}.MTF_shape_calc;
			ii = ii + 1;
		end
	end
	fprintf('%s done, %d percent done\n', putative, round(isesh/num_index*100))
end

% Save table
writetable(tables,'model_r2_values_ST2_MTF.xlsx')

%% Analysis 


isBE = strcmp(tables.MTF, 'BE');
isBS = strcmp(tables.MTF, 'BS');

lat_inh_BS = tables.Lat_Inh_R(isBS);
fprintf('Lateral Inhibition BS mean R^2 = %0.03f\n', median(lat_inh_BS));
lat_inh_BE = tables.Lat_Inh_R(isBE);
fprintf('Lateral Inhibition BE mean R^2 = %0.03f\n', median(lat_inh_BE));

SFIE_BS = tables.SFIE_R(isBS);
fprintf('SFIE BS mean R^2 = %0.03f\n', median(SFIE_BS));
SFIE_BE = tables.SFIE_R(isBE);
fprintf('SFIE BE mean R^2 = %0.03f\n', median(SFIE_BE));

SFIE_shape = strcmp(tables.MTF, tables.SFIE_MTFshape);
fprintf('SFIE: %d/%d, %0.1f%% correctly classified\n', sum(SFIE_shape), ...
    length(SFIE_shape), sum(SFIE_shape)/length(SFIE_shape)*100)
lat_inh_shape = strcmp(tables.MTF, tables.Lat_Inh_MTFshape);
fprintf('Lat Inh: %d/%d, %0.1f%% correctly classified\n', sum(lat_inh_shape), ...
    length(lat_inh_shape), sum(lat_inh_shape)/length(lat_inh_shape)*100)


%% -------------------------------------------------------------------------

function p_value = ftest(data, prediction1, prediction2)
% Assuming you have these variables:
% data: 40x1 vector of observed data
% gaussian_prediction: 40x1 vector of Gaussian model predictions
% dog_prediction: 40x1 vector of DoG model predictions

%mse_A = immse(data, prediction1);
%mse_B = immse(data, prediction2);
squared_errors_A = (data - prediction1).^2;
squared_errors_B = (data - prediction2).^2;
[~, p_value, ~, ~] = ttest(squared_errors_A, squared_errors_B);

end
