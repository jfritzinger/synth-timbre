%% Plots thresholds

% Load in thresholds
path = 'C:\Users\johan\Box\03 - Timbre Modeling\7 - Model\12-1 Results\';
filename = 'thresh.mat';
thresh = load([path filename]);

% Set axes 
DL = [1 2 5 10 20 50 100];
ylabels = [0.5 0.6 0.7 0.8 1 2.0 3.0 4.0 5.0 6.0 7.0 8.0 10 20.0 30.0 40.0 50.0];
xlabels = [1 2.0 5.0 10 20.0 50.0 100];

mean_thresh = mean(thresh.thresh, 1);
SD = std(thresh.thresh, 1);

% Plot
figure;
errorbar(DL, mean_thresh(2:end), SD(2:end), '-s','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
set(gca, 'YScale', 'log', 'Ylim', [0.5 50.0], 'YTick', ylabels, 'YTickLabel', ylabels);
set(gca, 'XScale', 'log', 'Xlim', [0.5 105.0], 'XTick', xlabels, 'XTickLabel', xlabels);
ylabel('Threshold (percent)')
xlabel('Variation in non-target dimension (multiples of DL)')
title('Experiment Results')
