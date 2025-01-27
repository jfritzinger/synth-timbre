%% plot_dog_vs_gaussian
clear 

%% Load 

load('R2_DOG.mat', "R2_gauss_all", "R2_dog_all")


%% Plot 

figure
scatter(R2_gauss_all, R2_dog_all, 'filled', 'MarkerEdgeColor','k')
xlim([0 1])
ylim([0 1])
hold on
plot([0 1], [0 1], 'k')
xlabel('Gaussian')
ylabel('DoG')
title('Adjusted R^2 Fits')
set(gca, 'fontSize', 16)

%% Get p-value from the f-test in each neuron to compare goodness of fit 
% Need to go back and save each fit so we can compare for each neuron


