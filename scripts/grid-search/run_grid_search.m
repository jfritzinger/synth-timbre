%% run_grid_search
clear

%% Load in spreadsheet and data

% Set up save paths
[~, computer] = system('hostname');
if ismac
	savepath = '/Volumes/Synth-Timbre/data/manuscript/model-fits';
	addpath('/Users/jfritzinger/Projects/shared-models/efferent-model/')
elseif contains(computer, 'I1') % I1
	savepath = '\\NSC-LCARNEY-H2\Synth-Timbre\data\manuscript\model-fits';
else
	savepath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript\model-fits';
end

% Load in spreadsheet
[base, ~, ~, ~] = getPaths();
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, 'scripts/data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);

% Find all data with WB-TIN
has_WB = cellfun(@(s) contains(s, 'R'), sessions.WBTIN_Units);
ind_WB = find(has_WB);
num_neurons = length(ind_WB);

% Load in data
iWB = 1;
putative_timbre = sessions.Putative_Units{ind_WB(iWB)};
putative_WB = sessions.WBTIN_Units{ind_WB(iWB)};
CF = sessions.CF(ind_WB(iWB));
base = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/WB-TIN/data/manuscript/Neural_Data';
load(fullfile(base, [putative_WB '.mat']), 'data');
mkdir(savepath, putative_timbre)

%% Create stimuli
params = cell(4, 1);

% MTFN parameters
params{1} = data{3,2}; % Gets binaural WB-TIN stimuli
params{1}.Fs = 100000;
params{1}.dur = 1;
params{1}.mnrep = 1;
params{1} = generate_MTF(params{1});
params{1}.num_stim = size(params{1}.stim, 1);

% WB-TIN parameters
for ii = 1:3
	data_WB = data{5+ii,2};
	if ~isempty(data_WB)
		params{1+ii} = data_WB; % Gets binaural WB-TIN stimuli
		params{1+ii}.Fs = 100000;
		params{1+ii}.mnrep = 1;
		params{1+ii}.physio = 1;
		params{1+ii} = generate_WBTIN(params{1+ii});
		params{1+ii}.num_stim = size(params{1+ii}.stim, 1);
	end
end
num_stim = size(params, 1);

%% Run models 

% Run AN model
run_AN_model(params, CF)

% Run IC models with grid search parameters
run_IC_model(putative_timbre)


%% Evaluate models 

eval_models

% Create PDF
pdf_grid_search(putative_timbre, CF)


