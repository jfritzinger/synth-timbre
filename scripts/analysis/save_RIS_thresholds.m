%% save_peak_finding_results.m
%
% Script that...
%
%
% Author: J. Fritzinger
% Created: 2026-03-11
%
% -------------------------------------------------------------------------
clear

% --- Replacement for setupUmichClusters ---
% Get number of processors from Slurm environment
if ~isempty(getenv('SLURM_CPUS_ON_NODE'))
    NP = str2num(getenv('SLURM_CPUS_ON_NODE'));
elseif ~isempty(getenv('SLURM_NTASKS'))
    NP = str2num(getenv('SLURM_NTASKS'));
else
    NP = 1; % Default to 1 if running locally
end

% Create the worker pool using the 'local' profile
% On a single Slurm node, 'local' uses the cores allocated to your job
if isempty(gcp('nocreate'))
    thePool = parpool('local', NP);
end

%% Add paths 

addpath(genpath('/home/jofritzi/projects/synth-timbre/scripts/SPIKY'))
addpath(genpath('/home/jofritzi/projects/synth-timbre/scripts/helper-functions'))
addpath(genpath('/home/jofritzi/projects/shared-physio'), '-end')

%%

% Load in spreadsheet
[base, datapath, ~, ~] = getPaths();
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Create table
varNames = ["Putative", "CF", "MTF","MTF_at200", "MTF_str", ...
    "SPL", "binmode", "F0", ...
    "Type", "Prom", "Width", "Lim", "Freq", "Q", "Q_log"];
varTypes = ["string", "double", "string", "string", "double", ...
    "double", "double", "double", ...
    "string", "double", "double", "double", "double", "double", "double"];
est_num_rows = 830; % set to number larger than
num_cols = length(varNames);
table_size = [est_num_rows num_cols];
tables = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames);

% Create has_data
data_ind(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
data_ind(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
% data_ind(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
% data_ind(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
% data_ind(:,5) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_con);
% data_ind(:,6) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con);
% data_ind(:,7) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con);
% data_ind(:,8) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_con);
% data_ind(:,9) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_100);
% data_ind(:,10) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_100);
% data_ind(:,11) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_100);
has_data = any(data_ind(:,2),2);
data_index = find(has_data);
num_neurons = sum(has_data);

%% Plot each neuron

row_idx = 1;
for isesh = 1:num_neurons
    ineuron = data_index(isesh);

    % Load in session
    putative = sessions.Putative_Units{ineuron};
    CF = sessions.CF(ineuron);
    MTF_shape = sessions.MTF{ineuron};
    at200 = sessions.MTF_at200{ineuron};
    load(fullfile(datapath, 'neural_data', [putative '.mat']))

    if CF<2000
        CF_Group = 'Low';
    elseif CF>=2000 && CF<4000
        CF_Group = 'Med';
    else
        CF_Group = 'High';
    end

    % MTF analysis
    params_MTF = data{3, 2};
    if ~isempty(params_MTF)
        data_MTF = analyzeMTF(params_MTF);
    end

    for idata = 2
        if data_ind(ineuron,idata)==1

            % Load in proper dataset for each idata
            if ismember(idata, [1, 2, 3, 4])
                param_ST = data(5+idata, 2);
            elseif ismember(idata, [5, 6, 7, 8])
                param_ST = data(1+idata, 1);
            else
                param_ST = data(1+idata, 2);
            end

            % Analyze synthetic timbre
            spl = param_ST{1}.spl;
            data_ST = analyzeST(param_ST, CF);
            data_ST = data_ST{1};


            temporal = analyzeST_Temporal(param_ST{1}, data_ST);

            kk = 1;
            nfpeaks = length(data_ST.fpeaks);
            nrep = param_ST{1}.nrep;
            for iii = 1:nfpeaks
                rep = temporal.y{iii};
                spike_times = temporal.x{iii}/1000;
                for jj = 1:nrep
                    ind_target = rep == jj;
                    spike_times_per_rep{jj, iii} = spike_times(ind_target);
                    rearr_trains{kk} = spike_times(ind_target);
                    kk = kk+1;
                end
            end

            % 1. Pre-calculate total number of comparisons
            num_trains = length(rearr_trains);
            total_pairs = num_trains * (num_trains - 1) / 2;

            % 2. Pre-allocate linear arrays for the results
            % Using a linear array avoids "slicing" errors in parfor
            dist_results = zeros(1, total_pairs);
            idx1_list = zeros(1, total_pairs);
            idx2_list = zeros(1, total_pairs);

            % 3. Create a map of which pair corresponds to which indices
            % This small loop is very fast and happens before the parfor
            count = 1;
            for i = 1:(num_trains-1)
                for j = (i+1):num_trains
                    idx1_list(count) = i;
                    idx2_list(count) = j;
                    count = count + 1;
                end
            end

            % 4. Run the parallel loop
            parfor p = 1:total_pairs
                % Get the specific indices for this worker
                i = idx1_list(p);
                j = idx2_list(p);

                train1 = rearr_trains{i};
                train2 = rearr_trains{j};

                % Prepare local pair
                spikes = zeros(2, max([length(train2), length(train1)]));
                spikes(1, 1:length(train1)) = train1;
                spikes(2, 1:length(train2)) = train2;

                % Setup local parameter structure
                para = struct('tmin', 0, 'tmax', 300, 'dts', 1, 'select_measures', [0 0 1 0 0 0 0 0]);

                % Use the parallel-safe function we created earlier
                % This avoids the 'Transparency violation' error
                [spikes_proc, para_proc, ret] = SPIKY_check_spikes_parallel(spikes, para);

                if ret == 0
                    try
                        SPIKY_loop_results = SPIKY_loop_f_distances(spikes_proc, para_proc);
                        dist_results(p) = SPIKY_loop_results.RI_SPIKE.matrix(1, 2);
                    catch
                        dist_results(p) = 0;
                    end
                else
                    dist_results(p) = 0;
                end
            end

            % 5. Reconstruct the 2D RI_S_dist matrix from the linear results
            RI_S_dist = zeros(num_trains, num_trains);
            for p = 1:total_pairs
                RI_S_dist(idx1_list(p), idx2_list(p)) = dist_results(p);
            end

            % Save out matrix
            savefile = fullfile(datapath, 'RIS_thresholds', [putative '_RIS.mat']);
            save(savefile, 'RI_S_dist', 'putative', 'param_ST', 'CF')

            [threshold_percent, threshold_freq, d_prime] = calculate_SPIKE_Thresholds(param_ST{1}.fpeaks, RI_S_dist, nrep);


            % Add data to table
            % if strcmp(type, 'Flat')
            % 	continue
            % else
            tables.Putative{row_idx} = sessions.Putative_Units{ineuron};
            tables.CF(row_idx) = CF;
            tables.CF_Group{row_idx} = CF_Group;
            tables.MTF{row_idx} = MTF_shape;
            tables.MTF_at200{row_idx} = at200;
            tables.MTF_str(row_idx) = data_MTF.perc_change;
            tables.SPL(row_idx) = spl;
            %tables.Type{ii} = type;
            tables.binmode(row_idx) = param_ST{1}.binmode;
            tables.F0(row_idx) = param_ST{1}.Delta_F;
            % tables.Width(ii) = width;
            % tables.Lim(ii) = lim;
            % tables.Prom(ii) = prom;
            % tables.Freq(ii) = freq;
            % tables.Q(ii) = freq/width;
            % tables.Q_log(ii) = log10(freq/width);
            tables.D_prime(row_idx) = max(d_prime, [], 'all');
            tables.Threshold(row_idx) = threshold_percent;
            tables.Thresh_Freq{row_idx} = threshold_freq;
            %tables.Slope_Rate{row_idx} = slope_rate;
            row_idx = row_idx+1;
            % end
        end
    end
    fprintf('%s done, %d percent done\n', putative, round(isesh/num_neurons*100))
    writetable(tables,fullfile(datapath, 'peak_picking_w_thresholds_RIS.xlsx'))
end

%% Save table

% Save table
writetable(tables,fullfile(datapath, 'peak_picking_w_thresholds_RIS.xlsx'))

delete(thePool)

