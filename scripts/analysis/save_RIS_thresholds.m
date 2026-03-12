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

            % [~, peak_ind] = max(data_ST.rate);
            % peak_f = data_ST.fpeaks(peak_ind);
            % [~, peaksm_ind] = max(data_ST.rates_sm);
            % peak_fsm = data_ST.fpeaks(peaksm_ind);

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

            % reps = 20;
            % fpeaks = 41;
            % RI_S_dist = zeros(reps*fpeaks,reps*fpeaks);
            % RI_S_dist_per = NaN(1, 190);
            % mean_RIS = zeros(1, fpeaks);
            % std_RIS = zeros(1, fpeaks);
            %
            % parfor iii = 1:fpeaks
            %
            %     kk = 1;
            %     spike_times = spike_times_per_rep(:,ii);
            %     for train1i = 1:(length(spike_times)-1)
            %         train1 = spike_times{train1i};
            %
            %         for train2i = (train1i+1):length(spike_times)
            %
            %             train2 = spike_times{train2i};
            %             spikes = zeros(2,max([length(train2) length(train1)])); %initialize matrix to hold both trains
            %             spikes(1,1:length(train1)) = train1;
            %             spikes(2,1:length(train2)) = train2;
            %
            %             %initialize parameter structure
            %             para = struct('tmin',[],'tmax',[],'dts',[],'select_measures',[]);
            %             para.select_measures = [0 0 1 0 0 0 0 0]; %RI_SPIKE
            %             para.tmin = 0;
            %             para.tmax = 300; % ms
            %             para.dts = 1;
            %             d_para = para;
            %
            %             % Call SPIKY
            %             SPIKY_check_spikes
            %             para = d_para;
            %             try
            %                 SPIKY_loop_results = SPIKY_loop_f_distances(spikes, para);
            %                 RI_S_dist_per(kk) = SPIKY_loop_results.RI_SPIKE.matrix(1,2);
            %             catch
            %                 RI_S_dist_per(kk) = 0;
            %             end
            %             kk = kk + 1;
            %         end
            %     end
            %     mean_RIS(iii) = mean(RI_S_dist_per, 'omitnan');
            %     std_RIS(iii) = std(RI_S_dist_per, 'omitnan');
            % end


            % % Pre-allocate final result containers
            % mean_RIS = zeros(1, nfpeaks);
            % std_RIS = zeros(1, nfpeaks);
            %
            % parfor iii = 1:nfpeaks
            %     % 1. Local Slicing: Extract the specific data for this worker
            %     % Note: Changed 'ii' to 'iii' to match the loop variable
            %     spike_times = spike_times_per_rep(:, iii);
            %
            %     % 2. Temporary Storage: Create a local version of RI_S_dist_per
            %     % This is initialized INSIDE the parfor so it's "Private" to each worker
            %     local_dist_per = NaN(1, 190);
            %     kk = 1;
            %     for train1i = 1:(nrep - 1)
            %         train1 = spike_times{train1i};
            %
            %         for train2i = (train1i + 1):nrep
            %             train2 = spike_times{train2i};
            %
            %             % Setup local spikes matrix
            %             spikes = zeros(2, max([length(train2), length(train1)]));
            %             spikes(1, 1:length(train1)) = train1;
            %             spikes(2, 1:length(train2)) = train2;
            %
            %             % Setup parameters
            %             para = struct('tmin', 0, 'tmax', 300, 'dts', 1, 'select_measures', [0 0 1 0 0 0 0 0]);
            %             d_para = para;
            %
            %             % Call SPIKY
            %             [spikes, d_para, ret] = SPIKY_check_spikes_parallel(spikes, d_para);
            %             para = d_para;
            %
            %             try
            %                 SPIKY_loop_results = SPIKY_loop_f_distances(spikes, para);
            %                 local_dist_per(kk) = SPIKY_loop_results.RI_SPIKE.matrix(1, 2);
            %             catch
            %                 local_dist_per(kk) = 0;
            %             end
            %             kk = kk + 1;
            %         end
            %     end
            %
            %     % 3. Calculate and store results in sliced output variables
            %     mean_RIS(iii) = mean(local_dist_per, 'omitnan');
            %     std_RIS(iii) = std(local_dist_per, 'omitnan');
            % end


            % Find peaks & prominence values
            % [peaks, dips, type, prom, width, lim, ~, ~, freq] = peakFinding(...
            %     data_ST, CF, 'Temporal', param_ST{1});

            % % Calculate thresholds
            % fpeaks = param_ST{1}.fpeaks;
            % [threshold_percent, threshold_freq, slope_rate, d_prime] = calculateThresholds(...
            %     fpeaks, mean_RIS, std_RIS, CF);
            %
            % figure
            % errorbar(fpeaks, mean_RIS, std_RIS/sqrt(190))
            % title(threshold_percent)
            % ylim([0 1])
            % hold on
            % xline(CF)
            % drawnow

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

