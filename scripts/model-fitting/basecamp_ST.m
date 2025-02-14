%% basecampe_ST
clear

%% Load in spreadsheet and data

% Set up save paths
[~, computer] = system('hostname');
if ismac
    savepath = '/Volumes/Synth-Timbre/data/manuscript/model-fits';
elseif contains(computer, 'I1') % I1
    savepath = '\\NSC-LCARNEY-H2\Synth-Timbre\data\manuscript\model-fits';
else
    savepath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript\model-fits';
end

% Load in spreadsheet
[~, datapath, ~, ~] = getPaths();
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);

% Find all data with WB-TIN
has_WB = cellfun(@(s) contains(s, 'R'), sessions.WBTIN_Units);
ind_WB = find(has_WB);
num_neurons = length(ind_WB);

% Load in data
for iWB = 3:num_neurons

    putative_timbre = sessions.Putative_Units{ind_WB(iWB)};
    putative = sessions.WBTIN_Units{ind_WB(iWB)};
    CF = sessions.CF(ind_WB(iWB));
    %base = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/WB-TIN/data/2024-manuscript/Neural_Data';
    base = 'C:\Users\jfritzinger\Box\02 - Code\WB-TIN\data\2024-manuscript\Neural_Data';
    load(fullfile(base, [putative '.mat']), 'data');

    if ~isempty(data{7,2})
        mkdir(savepath, putative_timbre)

        %% Create stimuli

        % MTFN parameters
        params{1} = data{3,2}; % Gets binaural WB-TIN stimuli
        params{1}.Fs = 100000;
        params{1}.dur = 0.7;
        params{1}.mnrep = 1;
        params{1} = generate_MTF(params{1});
        params{1}.num_stim = size(params{1}.stim, 1);

        % WB-TIN parameters
        % Fix so that only 40 dB SNR gets run
        params{2} = data{7,2}; % Gets binaural WB-TIN stimuli
        params{2}.Fs = 100000;
        params{2}.mnrep = 1;
        params{2}.physio = 1;
        params{2}.SNR = 40;
        params{2} = generate_WBTIN(params{2});
        params{2}.num_stim = size(params{2}.stim, 1);


        %% Run AN model

        paramCF = 0.125:0.125:1.5;
        run_AN_model(params, CF, paramCF, putative_timbre);

        %% Run IC model

        run_IC_model_fmincon(putative_timbre, data, CF)


        %% Run PDF results and evaluate fits

        get_best_fit_model(putative, CF, putative_timbre)
        pdf_evaluate_fits(putative, CF, putative_timbre)
    end
end