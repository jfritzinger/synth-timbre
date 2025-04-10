%% ST_dog%% dog_analysis_plots 
clear 

%% Load in fits  

[base, datapath, savepath, ppi] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")

%% Plot example fit 

fig = figure('Position',[53,540,8*ppi,7*ppi]);
tiledlayout(2, 2)
fontsize = 12;
legsize = 10;
labelsize = 20;

%% Plot example


% Example: R24_TT2_P13_N02, CF = 1150Hz, BS
%putative = 'R29_TT4_P5_N02'; 
%putative = 'R24_TT2_P13_N02';
%putative = 'R25_TT3_P9_N01';
putative = 'R27_TT3_P8_N01';
ind = cellfun(@(d) strcmp(d, putative), {dog_analysis.putative});

% Load in to get spont rate
load(fullfile(datapath, 'neural_data', [putative '.mat']))
params_RM = data{2, 2};
data_RM = analyzeRM(params_RM);
spont = data_RM.spont;

% Plot 
nexttile
hold on
plot(dog_analysis(ind).fpeaks/1000, dog_analysis(ind).rate, 'k','LineWidth',1);
% plot(dog_analysis(ind).fpeaks, dog_analysis(ind).dog_predicted2, ...
% 	'LineWidth',2, 'color', 'b');
plot(dog_analysis(ind).fpeaks/1000, dog_analysis(ind).dog_predicted, ...
	'LineWidth',1, 'color', 'b');
plot(dog_analysis(ind).fpeaks/1000, dog_analysis(ind).gaus_predicted,...
	'LineWidth',1, 'color', '#1b9e77');
xline(dog_analysis(ind).CF/1000, '--', 'LineWidth',1.5)
yline(spont, 'k')
ylabel('Avg. Rate (sp/s)')
xlabel('Spectral Peak Freq. (kHz)')
title('Example Neuron Predictions')
set(gca, 'fontSize', fontsize)
hleg = legend('Data', 'DoG', 'Gaussian', 'CF', 'fontsize', legsize);
hleg.ItemTokenSize = [16, 6];
grid on
xlim([dog_analysis(ind).fpeaks(1) dog_analysis(ind).fpeaks(end)]/1000)


%% Plot filters 

% Plot DoG Parameters
nexttile
DOGparams = dog_analysis(ind).dog_params;
W = dog_model(dog_analysis(ind).fpeaks, DOGparams);
hold on
plot(dog_analysis(ind).fpeaks/1000,W, 'color', 'b')
xline(dog_analysis(ind).CF/1000, '--', 'linewidth', 1.5)

Fs = 100000;
Gparams = dog_analysis(ind).gauss_params;
f = linspace(0, Fs/2, 100000);
fc = 10^Gparams(1);
sigma = 10^Gparams(2);
g = Gparams(3);
W = gaussian_model(f, fc, sigma, g);
hold on
plot(f/1000,W, 'color', '#1b9e77')

% Plot labels
title('DoG and Gaussian Kernels')
set(gca, 'fontsize', fontsize)
ylabel('Amplitude')
xlabel('Frequency (kHz)')
xlim([dog_analysis(ind).fpeaks(1) dog_analysis(ind).fpeaks(end)]/1000)
set(gca, 'xscale', 'log')
grid on

%% Plot adjusted R^2 values 

sig = [dog_analysis.p_value]<0.05;
notsig = [dog_analysis.p_value]>0.05;

nexttile
scatter(R2_gauss_all(sig), R2_dog_all(sig),30, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor','b')
hold on
scatter(R2_gauss_all(notsig), R2_dog_all(notsig), 30, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',[0.4 0.4 0.4])
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
yticks(0:0.2:1)
grid on
plot([0 1], [0 1], 'k')
xlabel('Gaussian Adjusted R^2')
ylabel('DoG Adjusted R^2')
title('Fit Comparison')
set(gca, 'fontSize', fontsize)
msg = sprintf('%d sig.', sum(sig));
msg2 = sprintf('%d not sig.', sum(notsig));
legend(msg, msg2, 'Location','southeast', 'fontsize', legsize)
set(gca, 'fontsize', fontsize)


%% Plot DoG parameter values 

% Get good fits
good_fit = [dog_analysis.R2_dog]>0.5;

% Plot 
all_dog_params = [dog_analysis(good_fit).dog_params];
all_dog_params = reshape(all_dog_params, 6,[])'; % 6 for OG
CFs = [dog_analysis(good_fit).CF];

% Un-log CF_exc, CF_inh
CF_exc = 10.^all_dog_params(:,5);
CF_inh = 10.^all_dog_params(:,6);
s_exc = 10.^all_dog_params(:,3);
s_inh = 10.^all_dog_params(:,4);
g_exc = all_dog_params(:,1);
g_inh = all_dog_params(:,2);

% Scatter plot of ratio of inhibition to excitation sigma and strengths 
ratio_sigma = log10(s_inh./s_exc);
ratio_g = log10(g_inh./g_exc);
nexttile
hold on
scatter(ratio_sigma, ratio_g, 30, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
xline(0)
yline(0)
% xlim([-12 6])
% ylim([-6 6.1])
xlabel('Log Bandwidth Ratio (\sigma_i_n_h/\sigma_e_x_c)')
ylabel('Log Strength Ratio (g_i_n_h/g_e_x_c)')
title('Fit Parameter Ratios')
grid on

%% Annotate 

% Set annotations
left = [0.03 0.51 0.03 0.51];
bottom = [0.95 0.95 0.47 0.47];
label = {'A', 'B', 'C', 'D'};
for ii = 1:4
	annotation('textbox',[left(ii) bottom(ii) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end

%% Save figure 

filename = 'Fig11_dog_analysis_plots';
saveFigure(filename)
