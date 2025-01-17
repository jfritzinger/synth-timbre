%% fit_gaussian_vs_dog

%% Load in spreadsheet 

addpath('/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions')
[base, ~, ~, ~] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Load in data and plot
linewidth = 1.5;

% Load in data
putative = 'R24_TT2_P13_N05';
[base, datapath, savepath, ppi] = getPaths();
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

% RM to get spont
params_RM = data{2, 2};
data_RM = analyzeRM(params_RM);
spont = data_RM.spont;

% Synthetic timbre analysis
params = data(7, 2);
params = params(~cellfun(@isempty, params));
data_ST  = analyzeST(params, CF);
data_ST = data_ST{1};
rate = data_ST.rate;
rate_std = data_ST.rate_std;
rlb = data_ST.rlb;
rub = data_ST.rub;
fpeaks = data_ST.fpeaks;
spl = data_ST.spl;
rate_sm = data_ST.rates_sm;
max_rate = max(rate);

% Plot
figure
hold on
rates_sm = smooth_rates(rate, rlb, rub, CF);
errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
	'linestyle', 'none', 'linewidth', 0.8, 'color', 'b')
plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color','b')
plot(fpeaks./1000, rates_sm, 'linewidth', linewidth, 'color', 'k')
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)
yline(0.1, 'k', LineWidth=linewidth)


%% Recreate stimulus (1 rep) 

params{1}.Fs = 100000;
params{1}.physio = 1;
params{1}.mnrep = 1;
params{1}.dur = 0.3;
params{1} = generate_ST(params{1});
params{1}.num_stim = size(params{1}.stim, 1);
Fs = 100000;
stim = params{1}.stim;
observed_rate = rate;
r0 = spont;

% % Plot spectrum of stimuli
% for istim = 1:40
% 	figure
% 	y = fft(stim(istim,:));
% 	n = length(stim(istim,:)); % number of samples
% 	fs = 100000;
% 	f = (0:n-1)*(fs/n); % frequency range
% 	power = abs(y).^2/n;
% 	plot(f, power)
% 	xlabel('Frequency')
% 	ylabel('Power')
% 	xlim([0 7000])
% end

%% fmincon

% Fit gaussian
init = [CF, 30, 0]; % Initial guess
lb = [0, 0, 0]; % Lower bounds
ub = [Inf, Inf, Inf]; % Upper bounds
gaussian_params = fmincon(@(p) ...
	objective_function(p, 'gaussian', Fs, stim, observed_rate, r0), ...
    init, [], [], [], [], lb, ub, []);

%% Fit DoG model
dog_init = [CF, 100, 200, 5000, 5000]; % Initial guess
dog_lb = [0, 0, 0, 0, 0]; % Lower bounds
dog_ub = [Inf, Inf, Inf, Inf, Inf]; % Upper bounds
dog_params = fmincon(@(p) objective_function(p, 'dog', Fs, stim, observed_rate, r0), ...
    dog_init, [], [], [], [], dog_lb, dog_ub, []);


%% Plot result

% Plot data 
figure
hold on
errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
	 'linewidth', linewidth, 'color', 'b')
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)
title('Gaussian vs DoG Fits')
ylabel('Avg. Rate (sp/s)')
xlabel('Spectral Peak Freq. (kHz)')

% Plot gaussian
f = linspace(0, Fs/2, 100000);
nstim = size(stim, 1);
gaus_predicted = zeros(nstim, 1);
for i = 1:nstim
	fc = gaussian_params(1);
	sigma = gaussian_params(2);
	g = gaussian_params(3);
	W = gaussian_model(f, fc, sigma, g);
	gaus_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
end
plot(fpeaks./1000, gaus_predicted, 'r', 'linewidth', linewidth)
gaussian_adj_r_squared = calculate_adj_r_squared(observed_rate,...
	gaus_predicted, 3);

% Plot DoG
f = linspace(0, Fs/2, 100000);
nstim = size(stim, 1);
dog_predicted = zeros(nstim, 1);
for i = 1:nstim
	fc = dog_params(1);
	sigma_e = dog_params(2);
	sigma_i = dog_params(3);
	ge = dog_params(4);
	gi = dog_params(5);
	W = dog_model(f, fc, sigma_e, sigma_i, ge, gi);
	dog_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
end
plot(fpeaks./1000, dog_predicted, 'g', 'linewidth', linewidth)
dog_adj_r_squared = calculate_adj_r_squared(observed_rate,...
	dog_predicted, 5);
set(gca, 'FontSize',16)
legend('Data', 'CF', 'Spont', 'Gaussian', 'DoG')


% Annotations
gaus_msg = sprintf('Gaussian adjusted R^2=%0.02f', gaussian_adj_r_squared);
text(0.05, 0.95, gaus_msg, 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',16)
dog_msg = sprintf('DoG adjusted R^2=%0.02f', dog_adj_r_squared);
text(0.05, 0.85, dog_msg, 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',16)

%% FUNCTIONS 

% Minimizing function 
function distance = objective_function(params, model, Fs, stim, observed_rate, r0)
    if strcmp(model, 'gaussian')
        fc = params(1);
        sigma = params(2);
        g = params(3);
        f = linspace(0, Fs/2, 100000);
        W = gaussian_model(f, fc, sigma, g);
    else % DoG model
        fc = params(1);
        sigma_e = params(2);
        sigma_i = params(3);
        ge = params(4);
        gi = params(5);
        f = linspace(0, Fs/2, 100000);
        W = dog_model(f, fc, sigma_e, sigma_i, ge, gi);
    end
    
	nstim = size(stim, 1);
    predicted_rate = zeros(nstim, 1);
    for i = 1:nstim
        predicted_rate(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
    end
    
    distance = sum(abs(predicted_rate - observed_rate)); % City-block distance
	fprintf('Fc = %0.0f, Sigma = %0.2f, g = %0.2f\n', params(1), params(2), params(3))
	fprintf('Dist = %0.2f\n', distance)

	% figure
	% plot(predicted_rate)
end


% Function to compute firing rate 
function r = compute_firing_rate(stim, Fs, W, f, r0)
    N = length(stim);
    X = fft(stim);
    P = abs(X/N).^2;
    P = P(1:N/2+1);
    P(2:end-1) = 2*P(2:end-1);
    
    f_signal = (0:(N/2))*Fs/N;
    P_interp = interp1(f_signal, P, f, 'linear', 0);
    
    r = sum(W .* P_interp) + r0;
    r = max(r, 0); % Half-wave rectification

	% figure; hold on
	% plot(f, P_interp)
	% plot(f, W)
	% ylim([0 0.0006])
	% xlim([0 2000])
	% disp(['Max P_interp: ', num2str(max(P_interp))]);
    % disp(['Max W: ', num2str(max(W))]);
    % disp(['Sum of W * P_interp: ', num2str(sum(W .* P_interp))]);
end



% Create Gaussian function 
function W = gaussian_model(f, fc, sigma, g)
    W = g * exp(-(f - fc).^2 / (2 * sigma^2));
end


% Difference of Gaussians (DoG) model
function W = dog_model(f, fc, sigma_e, sigma_i, ge, gi)
    W = ge * exp(-(f - fc).^2 / (2 * sigma_e^2)) - ...
        gi * exp(-(f - fc).^2 / (2 * sigma_i^2));
end


% Adjusted R^2 function
function adj_r_squared = calculate_adj_r_squared(observed, predicted, n_params)
    n = length(observed);
    SSE = sum((observed - predicted).^2);
    SST = sum((observed - mean(observed)).^2);
    adj_r_squared = 1 - (SSE / (n - n_params - 1)) / (SST / (n - 1));
end

