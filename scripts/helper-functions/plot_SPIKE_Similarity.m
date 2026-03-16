function plot_SPIKE_Similarity(fpeaks, RI_S_dist, reps)
% 1. Calculate the Mean Distance Matrix (Block Averaging)
num_f = length(fpeaks);
mean_dist_matrix = zeros(num_f, num_f);

for i = 1:num_f
    for j = 1:num_f
        idx1 = (i-1)*reps + (1:reps);
        idx2 = (j-1)*reps + (1:reps);
        % Average all trial-to-trial distances between Freq i and Freq j
        mean_dist_matrix(i,j) = mean(mean(RI_S_dist(idx1, idx2), 'omitnan'), 'omitnan');
    end
end

% 2. Create the Figure
figure('Color', 'w', 'Position', [100, 100, 600, 500]);

% Use imagesc for the heatmap
imagesc(fpeaks, fpeaks, mean_dist_matrix);

% Formatting
axis square;
colormap(flipud(hot)); % Hotter colors = further distance (more different)
colorbar;
set(gca, 'YDir', 'normal'); % Flip axis so low frequencies are at the bottom

% Labels
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('SPIKE-Distance Similarity Matrix', 'FontSize', 14);
hcb = colorbar;
title(hcb, 'Distance');

% Optional: Add a log scale if your frequencies are logarithmic
% set(gca, 'XScale', 'log', 'YScale', 'log');

grid on;
set(gca, 'GridColor', 'w', 'GridAlpha', 0.2);

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
end