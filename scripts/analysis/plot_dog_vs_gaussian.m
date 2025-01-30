%% plot_dog_vs_gaussian
clear 

%% Load 

[base, datapath] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")


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

num_sesh = length(dog_analysis);
p_all = NaN(num_sesh, 1);
CF_all = NaN(num_sesh, 1);
for ii = 1:num_sesh

	current_dog = dog_analysis(ii).dog_predicted;
	current_gauss = dog_analysis(ii).gaus_predicted;

	[h,p,ci,stats] = vartest2(current_gauss,current_dog);
	p_all(ii) = log(p);
	CF_all(ii) = dog_analysis(ii).CF;

end


figure('Position',[177,534,971,381])
tiledlayout(1, 2)
nexttile
scatter(CF_all, p_all, 'filled', "MarkerEdgeColor",'k')
hold on
yline(log(0.05), '--')
set(gca, 'XScale', 'log')
ylabel('Log p-value from F-test')
xlabel('CF')
set(gca, 'fontsize', 16)
title('Is DoG significantly different than Gaussian?')

nexttile
edges = linspace(-60, 0, 120);
histogram(p_all, edges)
xlabel('Log p-value from F-test')
ylabel('# Neurons')
set(gca, 'fontsize', 16)
title('Histogram of P-values')