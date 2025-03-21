%% plot_model_cell_temporal
clear 

%% Parameters
CF = 1200; 

% Stimulus parameters
params.fpeak_mid = 1200;
params.Delta_F = 200;
params.num_harms = 11;
params.stp_otc = 41;
params.Fs = 100000;
params.mnrep = 1;
params.physio = 0;
params.dur = 0.3;
params.ramp_dur = 0.02;
params.spl = 70;
params.g = 24;
params = generate_ST(params);
params.num_stim = size(params.stim, 1);

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CF_range = 1200;
model_params.num_CFs = 1;
model_params.CFs = 1200;
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 10; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw 
% implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - 
% this is the 'noise' associated with spont. activity of AN fibers - 
% see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to 
% omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium 
% SR, 3 = high SR)
model_params.Fs = 100000;


%% Model 

AN_HSR = modelAN(params, model_params); % HSR for IC input
SFIE = wrapperIC(AN_HSR.an_sout, params, model_params); % SFIE output
%save('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Synth-Timbre/Data/Intro_ModelResponse.mat', 'AN_HSR', 'SFIE')

% addpath '/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions'
% [base, datapath, ~] = getPaths();
% filename = 'Intro_ModelResponse.mat';
% load(fullfile(datapath, filename), 'AN_HSR', 'SFIE')

avAN = AN_HSR.average_AN_sout;
avBE = SFIE.average_ic_sout_BE;
avBS = SFIE.average_ic_sout_BS;
CFs = AN_HSR.CFs;

%% Plot stimulus 
figure('Position',[560,64,560,835])
tiledlayout(6, 1)
font_size = 12;

plot_range = [300 2400];
ticks = [200 500 1000 2000 5000];

% Stimulus Creation
nexttile
hold on
harmonics = params.Delta_F:params.Delta_F:10000;
num_harmonics = length(harmonics);
npts = params.dur * params.Fs; % # pts in stimulus
t = (0:(npts-1))/params.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/params.fpeak_mid)*...
	params.g)/20);
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
    comp_freq = harmonics(iharm);
    component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
    stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(params.spl/20) * (1/rms(stimulus));
component_scales_linear = Level_scale * component_scales_linear;

% Now, make the stimulus for this_fpeak
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
for iharm = 1:num_harmonics
    stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);
    
    % Plot each stimulus
    y = fft(stim);
    mdB = 20*log10(abs(y));
    level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
    stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
        3, 'Color', '#882255');
end

% Plot envelope of stimulus
plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
set(gca,'fontsize',font_size)
ylabel('Level (dB SPL)')
ylim([0 70])
grid on
set(gca, 'XScale', 'log')
hold off
xlim(plot_range)
xticks(ticks)
xlabel('Freq. (Hz)')
title('Stimulus')

% AN Plot
[rate, rate_std] = plotST(params, avAN, 0);
fpeaks = params.fpeaks;

nexttile
hold on
plot(fpeaks,rate , 'linewidth', 2, 'color', '#117733');
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 300])
xticks(ticks)
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
xlabel('CF (Hz)')
ylabel('Avg. Rate (sp/s)')
title('AN Avg. Rate')

% IC BE Plot
nexttile
hold on
[rate, rate_std] = plotST(params, avBE, 0);
plot(fpeaks, rate, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 50])
xticks(ticks)
yticks([0 30 60])
xticklabels([])
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
title('BE Avg. Rate')
ylabel('Avg. Rate (sp/s)')
xlabel('CF (Hz)')

% IC temporal plot 
fs = params.Fs;
spike_hist = squeeze(SFIE.ic_BE);
%VS = calcFFT(spike_hist, params, fs);
VS = calcVS(params, spike_hist, fs);
nexttile
plot(fpeaks, VS, 'linewidth', 2, 'color', '#d95f02')
xlim(plot_range)
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
set(gca, 'XScale', 'log')
xticks(ticks)
ylabel('|fft| @ 200 Hz')
ylim([0 1])
set(gca,'fontsize',font_size)
title('BE: fft @ 200 Hz')
xlabel('CF (Hz)')

% IC BS Plot
nexttile
[rate, rate_std] = plotST(params, avBS, 0);
hold on
plot(fpeaks, rate, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 50])
xticks(ticks)
yticks([0 30 60])
xlabel('CF (Hz)')
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
title('BS Avg. Rate')
ylabel('Avg. Rate (sp/s)')

% IC temporal plot 
fs = params.Fs;
spike_hist = squeeze(SFIE.ic_BS);
%VS = calcFFT(spike_hist, params, fs);
VS = calcVS(params, spike_hist, fs);
nexttile
plot(fpeaks, VS, 'linewidth', 2, 'Color','#d95f02')
xlim(plot_range)
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
set(gca, 'XScale', 'log')
xticks(ticks)
ylabel('|fft| @ 200 Hz')
ylim([0 1])
set(gca,'fontsize',font_size)
title('BS: fft @ 200 Hz')
xlabel('CF (Hz)')

%% Calculate PSTHs and period histograms for model PSTH 

figure('Position',[1700,173,463,1048])
tiledlayout(1, 6, 'TileSpacing','compact', 'Padding','compact')

num_fpeaks = length(fpeaks);
fs = 100000;
for ii = 1:2
	if ii == 1
		spike_hist = squeeze(SFIE.ic_BE);
		name = 'BE';
	else
		spike_hist = squeeze(SFIE.ic_BS);
		name = 'BS';
	end

	nexttile([1, 2])
	hold on
	max_rate = max(spike_hist, [], 'all')/2;
	onsetwin = 0.05;
	for j = 1:num_fpeaks

		% Plot PSTHs
		spike_wo_onset = spike_hist(j, onsetwin*fs:params.dur*fs-1);
		t = linspace(0.05, 0.3, fs*0.25);
		offset = max_rate * (j-1);
		plot(t, spike_wo_onset + offset);
	end
	ylim([0 max_rate*num_fpeaks])
	xlabel('Time (s)')
	box on
	yticks(linspace(max_rate/2, max_rate*100-max_rate/2, num_fpeaks))
	yticklabels(round(fpeaks))
	xlim([0.05 0.3])
	grid on
	title([name ', PSTH (excluding onset)'])
	set(gca, 'fontsize',12)

	% Calculate period histogram
	nexttile
	max_rate = max(spike_hist, [], 'all')/2;
	onsetwin = 0.05;
	hold on
	for j = 1:num_fpeaks

		% Plot PSTHs
		spike_wo_onset = spike_hist(j, onsetwin*fs:params.dur*fs-1);
		t = linspace(0.05, 0.3, fs*0.25);

		freq = 200; % Stimulus frequency in Hz
		period = 1 / freq; % Period in ms
		samples_per_period = fs*period;

		period_hist = reshape(spike_wo_onset, samples_per_period,[]);
		avg = mean(period_hist, 2);
		t_period = linspace(0, period, fs*period);

		offset = max_rate * (j-1);
		patch([t_period flip(t_period)],[avg+offset; ...
			repmat(offset,samples_per_period, 1)], 'k')

	end
	ylim([0 max_rate*num_fpeaks])
	xlabel('Time (s)')
	box on
	yticks(linspace(max_rate/2, max_rate*100-max_rate/2, num_fpeaks))
	yticklabels(round(fpeaks))
	xlim([0 0.005])
	xticks(0:0.001:0.005)
	grid on
	title([name ', Period Histogram'])
	set(gca, 'fontsize',12)
end

%% FUNCTIONS 

function R = calcVS(params, spike_hist, fs)
t = linspace(0, 0.25, fs*0.25);
f = 200;
onsetwin = 0.05; % ms
for ii = 1:length(params.fpeaks)
	r = spike_hist(ii, onsetwin*fs:params.dur*fs-1);
	R(ii) = abs(1/sum(r) * sum(r .* exp(1i * 2*pi * f .* t)));
end
end

% function VS = calcFFT(spike_hist, params, fs)
% 
% freq = 200;
% for ii = 1:length(params.fpeaks)
% 
% 	% Cut to integer # of cycles (50 cycles, excluding onset)
% 	onsetwin = 0.05; % ms
% 	spike_wo_onset = spike_hist(ii, onsetwin*fs:params.dur*fs-1);
% 
% 	% Find frequency index for 200 Hz
% 	L = length(spike_wo_onset);
% 
% 	% Take FFT 
% 	Y = fft(spike_wo_onset);
% 	P2 = abs(Y)/L; % Normalize by signal length
% 	P1 = P2(1:L/2+1);
% 	P1(2:end-1) = 2*P1(2:end-1);
% 	f = fs/L*(0:(L/2));
% 
% 	% Get FFT at 200 Hz 
% 	ind = f==freq;
% 	VS(ii) = P1(ind);
% 
% end
% end