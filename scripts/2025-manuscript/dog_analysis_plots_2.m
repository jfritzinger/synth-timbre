%% ST_dog%% dog_analysis_plots 
clear 

%% Load in fits  

[base, datapath] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")


%% Plot example fit 

figure('Position',[53,540,450,500])
tiledlayout(2, 1)
fontsize = 16;

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
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).rate, 'k','LineWidth',2);
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).dog_predicted, ...
	'LineWidth',2, 'color', 'b');
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).gaus_predicted,...
	'LineWidth',2, 'color', '#1b9e77');
xline(dog_analysis(ind).CF, '--', 'LineWidth',2)
yline(spont, 'k')
ylabel('Avg. Rate (sp/s)')
xlabel('Spectral Peak Freq. (Hz)')
title('Example Fit')
set(gca, 'fontSize', fontsize)
hleg = legend('Data', 'DoG', 'Gaussian', 'CF', 'fontsize', 16);
hleg.ItemTokenSize = [16, 6];
grid on
xlim([dog_analysis(ind).fpeaks(1) dog_analysis(ind).fpeaks(end)])

%% Plot adjusted R^2 values 

sig = [dog_analysis.p_value]<0.05;
notsig = [dog_analysis.p_value]>0.05;

nexttile
scatter(R2_gauss_all(sig), R2_dog_all(sig),60, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor','b')
hold on
scatter(R2_gauss_all(notsig), R2_dog_all(notsig), 60, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',[0.4 0.4 0.4])
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
yticks(0:0.2:1)
grid on
plot([0 1], [0 1], 'k')
xlabel('Gaussian Adjusted R^2')
ylabel('DoG Adjusted R^2')
title('DoG vs Gaussian Comparison')
set(gca, 'fontSize', fontsize)
msg = sprintf('%d sig.', sum(sig));
msg2 = sprintf('%d not sig.', sum(notsig));
legend(msg, msg2, 'Location','southeast')


%% Plot DoG parameter values 
% DoG params: g_exc, g_inh, s_exc, s_inh,  CF_exc, CF_inh

% Get good fits
good_fit = [dog_analysis.R2_dog]>0.5;


% Plot 
all_dog_params = [dog_analysis(good_fit).dog_params];
all_dog_params = reshape(all_dog_params, 6,[])';
CFs = [dog_analysis(good_fit).CF];

% Un-log CF_exc, CF_inh
CF_exc = 10.^all_dog_params(:,5);
CF_inh = 10.^all_dog_params(:,6);
s_exc = 10.^all_dog_params(:,3);
s_inh = 10.^all_dog_params(:,4);
g_exc = all_dog_params(:,1);
g_inh = all_dog_params(:,2);

figure('Position',[560,12,609,836])
tiledlayout(3, 2)
for ii = 1:6
	nexttile 
	if ii == 5
		scatter(CFs, CF_exc, 'filled')
		ylabel('CF exc (Hz)')
		hold on 
		plot([10,10000], [10 10000], 'k')
		xlim([300 10000])
	elseif ii == 6
		scatter(CFs, CF_inh, 'filled')
		ylabel('CF inh (Hz)')
		hold on 
		plot([10,10000], [10 10000], 'k')
		xlim([300 10000])
	elseif ii == 3
		scatter(CFs, s_exc, 'filled')
		ylabel('sigma exc (Hz)')
		hold on
		plot([10,10000], [10 10000], 'k')
		xlim([300 10000])
		ylim([0 10^5])
	elseif ii == 4
		scatter(CFs, s_inh, 'filled')
		ylabel('sigma inh (Hz)')
		hold on
		plot([10,10000], [10 10000], 'k')
		xlim([300 10000])
		ylim([0 10^5])
	elseif ii == 2
		scatter(CFs, all_dog_params(:,ii), 'filled')
		ylabel('g inh')
		title('Inhibitory Parameters')
	else
		scatter(CFs, all_dog_params(:,ii), 'filled')
		ylabel('g exc')
		title('Excitatory Parameters')
	end
	xlabel('CF')
	set(gca, 'xscale', 'log')
	set(gca, 'yscale', 'log')
	set(gca, 'fontsize', 16)
end

%% Histograms 

figure
tiledlayout(3, 2)

nexttile
histogram(g_exc, 21)
ylabel('# neurons')

nexttile
histogram(g_inh, 21)
ylabel('# neurons')

nexttile
edges = 10.^linspace(log10(10),log10(10000),21);
histogram(s_exc, edges)
set(gca, 'xscale', 'log')
xlabel('sigma exc')
ylabel('# neurons')

nexttile
histogram(s_inh, edges)
set(gca, 'xscale', 'log')
xlabel('sigma inh')
ylabel('# neurons')

nexttile
histogram(CF_exc, edges)
set(gca, 'xscale', 'log')
xlabel('CF exc')
ylabel('# neurons')

nexttile
histogram(CF_inh, edges)
set(gca, 'xscale', 'log')
xlabel('CF inh')
ylabel('# neurons')
%%
figure
histogram(s_inh./s_exc, edges)

%%
% Scatter plot of ratio of inhibition to excitation sigma and strengths 
ratio_sigma = log10(s_inh./s_exc);
ratio_g = log10(g_inh./g_exc);
figure
hold on
scatter(ratio_sigma, ratio_g, 'filled')
xline(0)
yline(0)
% xlim([-12 6])
% ylim([-6 6.1])
xlabel('Log Bandwidth Ratio (\sigma_i_n_h/\sigma_e_x_c)')
ylabel('Log Strength Ratio (g_i_n_h/g_e_x_c)')



