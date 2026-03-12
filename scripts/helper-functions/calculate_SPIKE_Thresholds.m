function [threshold_percent, threshold, d_prime_results] = calculate_SPIKE_Thresholds(fpeaks, RI_S_dist, reps)
% fpeaks: vector of your 41 frequencies
% RI_S_dist: the (reps*fpeaks) x (reps*fpeaks) distance matrix
% reps: 20

num_f = length(fpeaks);
% This will hold every possible d' combination
d_prime_matrix = NaN(num_f, num_f);

% 1. Calculate d' for every possible pair of frequencies
for i = 1:num_f
    % Get indices for Frequency i (the "Reference")
    idx1 = (i-1)*reps + (1:reps);
    
    % Get Within-Group Distances (the "Noise" distribution)
    within_mat = RI_S_dist(idx1, idx1);
    % We need the distances between different trials of the same frequency
    % Use the upper triangle to avoid the 0-diagonal and duplicate pairs
    within_vals = within_mat(triu(true(size(within_mat)), 1));
    
    mu_within = mean(within_vals, 'omitnan');
    std_within = std(within_vals, 'omitnan');
    
    for j = i+1:num_f
        % Get indices for Frequency j (the "Comparison")
        idx2 = (j-1)*reps + (1:reps);
        
        % Get Between-Group Distances (the "Signal" distribution)
        % This is the block comparing trials of Freq i to trials of Freq j
        between_vals = RI_S_dist(idx1, idx2);
        mu_between = mean(between_vals(:), 'omitnan');
        
        % Calculate d'
        % Logic: How much did the distance increase relative to internal jitter?
        d_prime_matrix(i,j) = (mu_between - mu_within) / std_within;
    end
end

% Set the output variable
d_prime_results = d_prime_matrix;

% 2. Search for the smallest frequency delta that reaches d' = 1
threshold_percent = NaN;
threshold = NaN;
found = false;

% We check diagonals: delta=1 is adjacent frequencies, delta=2 is skipping one, etc.
for delta = 1:(num_f-1) 
    current_d_primes = diag(d_prime_matrix, delta);
    
    % Find if any pair at this 'distance' apart in frequency hits the threshold
    if any(current_d_primes >= 1)
        idx_match = find(current_d_primes >= 1, 1);
        f1 = fpeaks(idx_match);
        f2 = fpeaks(idx_match + delta);
        
        threshold = [f1, f2];
        freq_diff = f2 - f1;
        fpeak_mid = (f1 + f2) / 2;
        threshold_percent = (freq_diff / fpeak_mid) * 100;
        
        fprintf('SPIKE d'' = %.2f\n', current_d_primes(idx_match));
        disp(['SPIKE Threshold = ' num2str(threshold_percent) '%'])
        found = true;
        break;
    end
end

if ~found
    disp('No frequency pair reached the d'' = 1 threshold.');
end

end