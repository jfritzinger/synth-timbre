%% dog_analysis_plots 
clear 

%% Load in fits  

[base, datapath] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")


%% Plot example fit 

figure('Position',[53,540,1299,377])
tiledlayout(1, 3)

% Example: R24_TT2_P13_N02, CF = 1150Hz, BS
putative = 'R24_TT2_P13_N02';
ind = cellfun(@(d) strcmp(d, putative), {dog_analysis.putative});

% Load in to get spont rate
load(fullfile(datapath, 'neural_data', [putative '.mat']))
params_RM = data{2, 2};
data_RM = analyzeRM(params_RM);
spont = data_RM.spont;

% Plot 
nexttile
hold on
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).rate, 'LineWidth',2);
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).dog_predicted, 'LineWidth',2);
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).gaus_predicted, 'LineWidth',2);
xline(dog_analysis(ind).CF, '--')
yline(spont, 'k')
ylabel('Avg. Rate (sp/s)')
xlabel('Spectral Peak Freq. (Hz)')
title('Example Fits')
set(gca, 'fontSize', 16)
legend('Data', 'DoG', 'Gaussian')

%% Plot adjusted R^2 values 

nexttile
scatter(R2_gauss_all, R2_dog_all, 'filled', 'MarkerEdgeColor','k')
xlim([0 1])
ylim([0 1])
hold on
plot([0 1], [0 1], 'k')
xlabel('Gaussian')
ylabel('DoG')
title('Adjusted R^2 Fits')
set(gca, 'fontSize', 16)


%% Plot F-test for MSE of fits 

num_sesh = length(dog_analysis);
p_all = NaN(num_sesh, 1);
CF_all = NaN(num_sesh, 1);
for ii = 1:num_sesh

	current_dog = dog_analysis(ii).dog_predicted;
	current_gauss = dog_analysis(ii).gaus_predicted;
	rate = dog_analysis(ii).rate;
	CF_all(ii) = dog_analysis(ii).CF;
	
	% Comparing based on how close the curves are to data 
	p_value = ftest(rate, current_gauss, current_dog);
	p_all(ii) = log(p_value);
	
end

sig = p_all<log(0.05);
not_sig = p_all>log(0.05);

nexttile
scatter(CF_all(sig), p_all(sig), 'filled', "MarkerEdgeColor",'k')
hold on
scatter(CF_all(not_sig), p_all(not_sig), 'filled', "MarkerEdgeColor",'k')
yline(log(0.05), '--')
set(gca, 'XScale', 'log')
ylabel('Log p-value from F-test')
xlabel('CF')
set(gca, 'fontsize', 16)
title('Is DoG significantly different than Gaussian?')

% Annotations 
msg = sprintf('%d/%d significant', sum(sig), length(sig));
msg2 = sprintf('%d/%d not significant', sum(not_sig), length(not_sig));
text(0.65, 0.15, msg, 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',16)
text(0.65, 0.09, msg2, 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',16)


%% 

function p_value = ftest(data, prediction1, prediction2)
% Assuming you have these variables:
% data: 40x1 vector of observed data
% gaussian_prediction: 40x1 vector of Gaussian model predictions
% dog_prediction: 40x1 vector of DoG model predictions

% Calculate RSS
RSS_gaussian = sum((data - prediction1).^2);
RSS_dog = sum((data - prediction2).^2);

% Degrees of freedom
df_gaussian = length(data) - 3;
df_dog = length(data) - 6;

% Calculate F-statistic
F = ((RSS_gaussian - RSS_dog) / (df_gaussian - df_dog)) / (RSS_dog / df_dog);

% Compute p-value
p_value = 1 - fcdf(F, df_gaussian - df_dog, df_dog);

% Display results
% fprintf('F-statistic: %.4f\n', F);
% fprintf('p-value: %.4e\n', p_value);

end
