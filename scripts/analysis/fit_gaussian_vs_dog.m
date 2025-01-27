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

% Generate stimulus 
params{1}.Fs = 100000;
params{1}.physio = 1;
params{1}.mnrep = 1;
params{1}.dur = 0.3;
params{1} = generate_ST(params{1});
params{1}.num_stim = size(params{1}.stim, 1);
Fs = 100000;
observed_rate = rate;
r0 = spont;

% Run AN model 
model_params.type = 'SFIE';
model_params.range = 1; % 1 = population model, 2 = single cell model
model_params.species = 2; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CFs = CF;
model_params.nAN_fibers_per_CF = 10;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 5; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

timerVal = tic;
AN = modelAN(params{1}, model_params); % HSR for IC input
elapsedTime = toc(timerVal)/60;
disp(['Model took ' num2str(elapsedTime) ' minutes'])

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
type = 2; % 1: distance, 2: MSE
stim = params{1}.stim;
%stim = squeeze(AN.an_sout);
timerVal = tic;

% Fit gaussian
init = [CF, 30, 100]; % Initial guess
lb = [CF/4, 0, 0]; % Lower bounds
ub = [CF*4, 5000, Inf]; % Upper bounds

options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
	'MaxFunEvals', 10^12, 'maxiterations', 1000, 'ConstraintTolerance', 1e-12, ...
	'StepTolerance', 1e-16, 'display', 'off');
gaussian_params = fmincon(@(p) ...
	dog_objective_function(p, 'gaussian', Fs, stim, observed_rate, r0, type), ...
    init, [], [], [], [], lb, ub, [], options);


% Fit DoG model
dog_init = [20000, 10000, 100, 500,  CF, CF]; % Initial guess
dog_lb = [100,   100,     10,  10,   CF/4, CF/4]; % Lower bounds
dog_ub = [30000, 30000,   1000,1000, CF*4, CF*4]; % Upper bounds

options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
	'MaxFunEvals', 10^12, 'maxiterations', 1000, 'ConstraintTolerance', 1e-12, ...
	'StepTolerance', 1e-16, 'display', 'off');
dog_params = fmincon(@(p) dog_objective_function(p, 'dog', Fs, stim, observed_rate, r0, type), ...
    dog_init, [], [], [], [], dog_lb, dog_ub, [], options);
disp(['Model took ' num2str(toc(timerVal)) ' seconds'])

% Plot result

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
	W = dog_model(f, dog_params);
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
