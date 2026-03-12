%% SPIKY_metrics
clear


%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Load in example

putative = 'R27_TT2_P8_N02'; %'R24_TT2_P13_N05'; %
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

%% Analysis

% Synthetic timbre analysis
params = data(7, 2);
params = params(~cellfun(@isempty, params));
data_ST  = analyzeST(params, CF);
data_ST = data_ST{1};
param = params{1};
temporal = analyzeST_Temporal(param, data_ST);

%% Plot PSTH and smoothed PSTH

max_rate = max(temporal.PSTH, [], 'all');
freq_values = round(param.fpeaks);

figure('Position',[200,22,403,1233])
tiledlayout(1, 2, 'Padding','compact')

% Plot PSTH full response
nexttile
num_fpeaks = length(freq_values);
hold on
for j = 1:num_fpeaks

	% Plot PSTHs
	counts = temporal.PSTH(j,:);
	edges = temporal.t;
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');

	% Plot smoothed PSTH
	plot(temporal.t(1:end-1), temporal.PSTH_smooth(j,:)+offset,'k', 'LineWidth',1.5);

end
ylim([0 max_rate*num_fpeaks])
xlabel('Time (ms)')
box on
yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
yticklabels(freq_values)
title('PSTH')
set(gca, 'fontsize', 12)

% Plot period histogram
nexttile
max_rate = max(temporal.p_hist, [], 'all');
for j = 1:num_fpeaks

	% Plot PSTHs
	counts = temporal.p_hist(j,:);
	edges = temporal.t_hist;
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
end
ylim([0 max_rate*num_fpeaks])
xlabel('Time (ms)')
box on
yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
yticklabels(freq_values)
xlim([0 5])
xticks(1:5)
grid on
title('Period Histogram')
set(gca, 'fontsize', 12)


%% Cell 'spikes' with two or more spike trains (each cell array contains the spike times of one spike train)

kk = 1;
for ii = 1:41
    rep = temporal.y{ii};
    spike_times = temporal.x{ii}/1000;
    for jj = 1:20
        ind_target = rep == jj;
        spike_times_per_rep{jj, ii} = spike_times(ind_target);
        rearr_trains{kk} = spike_times(ind_target);
        kk = kk+1;
    end
end



%%

%start timing calc
reps = 20;
fpeaks = 41;
RI_S_dist = zeros(reps*fpeaks,reps*fpeaks);

%do RI_SPIKE distance on all
% for train1i = 1:(length(rearr_trains)-1)
%     train1 = rearr_trains{train1i};
% 
%     for train2i = (train1i+1):length(rearr_trains)
% 
%         train2 = rearr_trains{train2i};
%         spikes = zeros(2,max([length(train2) length(train1)])); %initialize matrix to hold both trains
%         spikes(1,1:length(train1)) = train1;
%         spikes(2,1:length(train2)) = train2;
%         % spikes = cell(1, 2);
%         % spikes{1, 1} = train1;
%         % spikes{1, 2} = train2;
% 
%         %initialize parameter structure
%         para = struct('tmin',[],'tmax',[],'dts',[],'select_measures',[]);
%         para.select_measures = [0 0 1 0 0 0 0 0]; %RI_SPIKE
%         para.tmin = 0;
%         para.tmax = 300; % ms
%         para.dts = 1;
%         d_para = para;
% 
%         % Call SPIKY
%         SPIKY_check_spikes
%         para = d_para;
%         try
%             SPIKY_loop_results = SPIKY_loop_f_distances(spikes, para);
%             RI_S_dist(train1i,train2i) = SPIKY_loop_results.RI_SPIKE.matrix(1,2);
%         catch
%             RI_S_dist(train1i,train2i) = 0;
%         end
%     end
% end

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

% RI_S_dist_per = NaN(41, 20, 20);
% for ii = 1:fpeaks
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
%             % spikes = cell(1, 2);
%             % spikes{1, 1} = train1;
%             % spikes{1, 2} = train2;
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
%                 RI_S_dist_per(ii, train1i,train2i) = SPIKY_loop_results.RI_SPIKE.matrix(1,2);
%             catch
%                 RI_S_dist_per(ii, train1i,train2i) = 0;
%             end
%         end
%     end
%     mean_fpeak(ii) = mean(squeeze(RI_S_dist_per(ii,:,:)), 'all', 'omitnan');
% end

%%

[threshold_percent, threshold, d_prime_results] = calculate_SPIKE_Thresholds(param.fpeaks, RI_S_dist, reps);

%% Plot similarity matrix?

nfpeaks = 41;
figure('Position', [100 100 350 350]);
class_mat_alt = zeros(nfpeaks,nfpeaks);
full_dist_matrix = RI_S_dist;

%we must add the matrix and the transposed matrix to get the
%"full" matrix (not halved down the diagonal)
full_dist_matrix = full_dist_matrix' + full_dist_matrix;

%loop through each spike train
for ti = 1:height(full_dist_matrix)
    %determine the average distance for this train by vowel
    avg_dist = zeros(1,nfpeaks);
    for fpeakii = 1:nfpeaks
        these_dists = reps*(fpeakii-1)+1:reps*(fpeakii-1)+reps;

        %use logical indexing to remove the diagonal
        these_dists = these_dists(these_dists~=ti);
        this_vow_dist = full_dist_matrix(ti,these_dists);
        avg_dist(fpeakii) = mean(this_vow_dist);
    end

    %classify
    [~,assigned_fpeak] = min(avg_dist);
    original_fpeak = ceil(ti/reps);
    class_mat_alt(original_fpeak,assigned_fpeak) = class_mat_alt(original_fpeak,assigned_fpeak) + 1;
end
total_assignments = sum(class_mat_alt(1,:));
class_mat_plot = zeros(nfpeaks+1,nfpeaks+1);
class_mat_plot(1:nfpeaks,1:nfpeaks) = class_mat_alt/total_assignments;

ax = axes('Position',[0.15 0.20 0.65 0.65]);
pcolor(class_mat_plot')
clim([0 1])

%loop through class_mat and print percents
for yi = 1:nfpeaks
    for xi = 1:nfpeaks
        if class_mat_alt(yi,xi) > 0
            this_text = num2str(round(class_mat_alt(yi,xi),2));
            if xi == yi
                text(yi+0.5,xi+0.5,this_text,'HorizontalAlignment','center','FontSize',12,...
                    'FontWeight','bold')
            else
                text(yi+0.5,xi+0.5,this_text,'HorizontalAlignment','center','FontSize',12)
            end
        end
    end
end

ax.XAxis.TickValues = (1:nfpeaks)+0.5;
ax.YAxis.TickValues = (1:nfpeaks)+0.5;
title('Timing Confusion Matrix')
ylabel('Classified Spectral Peak')
xlabel('True Spectral Peak')
set(gca,'FontSize',12)

