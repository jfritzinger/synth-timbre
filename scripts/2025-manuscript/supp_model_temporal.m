%% supp_model_temporal
clear 

%% Parameters
CF = 1200; 

% Stimulus parameters
params.fpeak_mid = 1200;
params.Delta_F = 200;
params.num_harms = 11;
params.stp_otc = 1;
params.Fs = 100000;
params.mnrep = 1;
params.physio = 0;
params.dur = 0.3;
params.ramp_dur = 0.02;
params.spl = 70;
params.g = 24;
params = generate_ST(params);
fs = params.Fs;

%% Model 

addpath '/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions'
[base, datapath, ~] = getPaths();
filename = 'Intro_ModelResponse.mat';
load(fullfile(datapath, filename), 'AN_HSR', 'SFIE')
avAN = AN_HSR.average_AN_sout;
avBE = SFIE.average_ic_sout_BE;
avBS = SFIE.average_ic_sout_BS;
CFs = AN_HSR.CFs;

%% Create heatmaps for BE and BS 

figure('Position',[560,589,560,259])
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact')
for ii = 1:2
	if ii == 1
		spike_hist = squeeze(SFIE.ic_BE);
		name = 'BE';
	else
		spike_hist = squeeze(SFIE.ic_BS);
		name = 'BS';
	end

	% Calculate period histogram
	nexttile
	max_rate = max(spike_hist, [], 'all')/2;
	onsetwin = 0.05;
	hold on
	for j = 1:100

		% Plot PSTHs
		spike_wo_onset = spike_hist(j, onsetwin*fs:params.dur*fs-1);
		t = linspace(0.05, 0.3, fs*0.25);

		freq = 200; % Stimulus frequency in Hz
		period = 1 / freq; % Period in ms
		samples_per_period = fs*period;

		period_hist = reshape(spike_wo_onset, samples_per_period,[]);
		avg(j,1:500) = mean(period_hist, 2)';
		t_period = linspace(0, period, fs*period);

	end
	grayMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', linspace(0, 1, 256)'];
	grayMap = flipud(grayMap);

	hh = pcolor(t_period, CFs./1000, avg);
	set(hh, 'EdgeColor', 'none');
	hold on
	yline(CF/1000, 'r', 'LineWidth',2)
	colormap(grayMap);
	box on
	yticks(linspace(max_rate/2, max_rate*100-max_rate/2, 100))
	set(gca, 'YScale', 'log')
	ylim([CFs(1) CFs(end)]/1000)
	yticks([0.1 0.2 0.5 1 2 5 10])
	ylabel('CFs (kHz)')
	xlim([0 0.005])
	xticks(0:0.001:0.005)
	xticklabels(0:5)
	xlabel('Period (ms)')
	grid on
	title(name)
	set(gca, 'fontsize',12)
end

%% Save figure 

filename = 'FigS4_model_temporal';
saveFigure(filename)