%% Fig. 4: Comparing model results and recordings
% J. Fritzinger, updated 2/2/2021
clear

%% Run Model

param = cell(4,1);
for istim = 1:4

	if istim == 1
		param{istim}.fpeak_mid = 1200;
	elseif istim == 2
		param{istim}.fpeak_mid = 1100;
	elseif istim == 3
		param{istim}.fpeak_mid = 1000;
	else
		param{istim}.fpeak_mid = 1100;
	end

	param{istim}.Delta_F = 200;
	param{istim}.dur = 0.3;
	param{istim}.ramp_dur = 0.02;
	param{istim}.spl = 70;
	param{istim}.g = 24;
	param{istim}.num_harms = 9;
	param{istim}.stp_otc = 41;
	param{istim}.Fs = 100000;
	param{istim}.mnrep = 1;
	param{istim}.physio = 0;
	param{istim} = generate_ST(param{istim});
	param{istim}.num_stim = size(param{istim}.stim, 1);
end

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.num_CFs = 1;
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 3; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

% Run model
SFIE = cell(4,1);
for istim = 1:4
	model_params.CF_range = param{istim}.fpeak_mid;
	model_params.CFs = param{istim}.fpeak_mid;

	AN_HSR = modelAN(param{istim}, model_params); % HSR for IC input
	SFIE{istim} = wrapperIC(AN_HSR.an_sout, param{istim}, model_params); % SFIE output

	if istim == 1
		[avBS_singlepeak, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BS, 0);
	elseif istim == 2
		[avBS_doublepeak, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BS, 0);
	elseif istim == 3
		[avBE_singlepeak, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BE, 0);
	else
		[avBE_doublepeak, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BE, 0);
	end
end

%% Load in data

base = getPaths();
fpath = 'data/aro-2021/Figure 4';
file_names = dir(fullfile(base, fpath, 'R0*'));

for ii = 1:length(file_names)
	data(ii) = load(fullfile(base, fpath, file_names(ii).name));
end

% Single Peak Stimulus
figure('Position',[140,224,1327,535])
freq_limits = [600 1800];
rate_limits = [0 85];
subplot(4, 3, 2)

% Stimulus parameters
ylabel('Level (dB SPL)')
ylim([30 65])
xlim(freq_limits)
xticklabels([])
grid on
hold on
box on

% Stimulus Creation
harmonics = param{1}.Delta_F:param{1}.Delta_F:10000; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = param{1}.dur * param{1}.Fs; % # pts in stimulus
t = (0:(npts-1))/param{1}.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/param{1}.fpeak_mid)*param{1}.g)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
	comp_freq = harmonics(iharm);
	component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
	stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(param{1}.spl/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

% Now, make the stimulus for this_fpeak
for iharm = 1:num_harmonics
	stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);

	% Plot each stimulus
	y = fft(stim);
	mdB = 20*log10(abs(y));
	level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
	stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
		2, 'Color', [0.6350, 0.0780, 0.1840]);
end
plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
set(gca,'fontsize',14)

%  Double Peak Stimulus
subplot(4, 3, 3)
ylim([30 65])
xlim(freq_limits)
xticklabels([])
yticklabels([])
grid on
hold on
box on

% Stimulus Creation
harmonics = param{2}.Delta_F:param{2}.Delta_F:10000; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = param{2}.dur * param{2}.Fs; % # pts in stimulus
t = (0:(npts-1))/param{2}.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/param{2}.fpeak_mid)*param{2}.g)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
	comp_freq = harmonics(iharm);
	component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
	stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(param{2}.spl/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

% Now, make the stimulus for this_fpeak
for iharm = 1:num_harmonics
	stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);

	% Plot each stimulus
	y = fft(stim);
	mdB = 20*log10(abs(y));
	level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
	stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
		2, 'Color', [0.6350, 0.0780, 0.1840]);
end
plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
set(gca,'fontsize',14)

% RM
[spls,~,si] = unique(double([data(1).param.list.spl]).');
num_spls = length(spls);
[freqs_RM,~,fi] = unique(double([data(1).param.list.freq]).');
num_freqs = length(freqs_RM);
rate_size = [num_freqs,num_spls];
num_spikes = data(1).spike_info.num_spikes_peri;
spike_rates = num_spikes*1000/data(1).param.dur; % spikes/sec
[rates_RM,~,~,~] = accumstats({fi,si},spike_rates, rate_size);
spont = mean(rates_RM(:,1));
max_y = max(rates_RM(:));

% Plot
subplot(4, 3, 4)
hold on
plot(log10(freqs_RM([1 end])),[1 1]*spont,'-','LineWidth',2, 'Color',[0.7 0.7 0.7])
plot(log10(freqs_RM),rates_RM(:,5),'color', '#332288','LineWidth',2)
box on
set(gca,'XTick',log10([250 500 1000 2000 5000 10000]), 'XTickLabel',{'0.25','0.5','1','2','5','10'})
xlims = log10(freqs_RM([1 end]));
xlim(xlims.*[0.98;1.02])
ylim([0 max_y])
ylabel('Rate (sp/s)')
xlabel('Frequency at 70dB SPL')
grid on
hold off
legend('Spontaneous Rate', 'fontsize',10)
set(gca,'fontsize',14)


% MTF
dur = data(2).param.dur/1000; % stimulus duration in seconds.
all_mod_depths = double([data(2).param.list.mdepth]).';
all_mod_depths(all_mod_depths < -80) = max(all_mod_depths);
[~,~,mdi] = unique(all_mod_depths);
[fms,~,fmi] = unique(double([data(2).param.list.fm]).');
if fms(1) == 0
	fms(1) = 1.2;
end
num_mod_freqs = length(fms);
num_depths = length(data(2).param.all_mdepths);
map_size = [num_mod_freqs, num_depths];
spike_rates = data(2).spike_info.num_spikes_delayed/...
	(dur - data(2).spike_info.onsetWin/1000);
[rate_MTF,~,rlb,rub] = accumstats({fmi,mdi},spike_rates, map_size);
rate_sm = zeros(size(rate_MTF));
for j = 1:num_depths
	rate_sm(:,j) = smooth_rates(rate_MTF(:,j),rlb(:,j),rub(:,j));
end

% Plot
subplot(4, 3, [7 10])
hold on
line([1 fms(end)], [1 1]*rate_MTF(1),'Color',[0.7 0.7 0.7], 'LineWidth', 2);
line(fms,rate_MTF,'Color','#332288', 'Marker','.', 'MarkerSize',10, 'LineWidth', 2);
box on
grid on
xtick = [1 2 5 10 20 50 100 200 500];
xlim(xtick([1 end]))
xlabel('Modulation Freq (Hz)')
ylabel('Rate (sp/s)')
set(gca,'XTick',xtick,'XScale', 'log')
hold off
legend('No modulation', 'fontsize',10)
set(gca,'fontsize',14)

% Single Peak (Neuron 1, Model)
for ii = 3:4
	[fpeaks_SC(ii, :),~,fpeaksi] = unique([data(ii).param.list.fpeak].');
	num_fpeaks = length(fpeaks_SC);
	dur = data(ii).param.dur/1000; % stimulus duration in seconds.

	rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
	spike_rates = data(ii).spike_info.num_spikes_delayed/...
		(dur - data(ii).spike_info.onsetWin/1000);
	[rate_SC(ii,:),~, rlb, rub] = accumstats({fpeaksi},spike_rates, rate_size);
	for j = 1:num_depths
		rate_SC_sm(ii,:) = smooth_rates(rate_SC(ii,:),rlb,rub);
	end
end

subplot(4,3, [5 8 11])
hold on
plot(param{1}.fpeaks, avBS_singlepeak, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410])
plot(fpeaks_SC(4,:), rate_SC(4,:), 'LineWidth', 2);
xlim(freq_limits)
ylim(rate_limits)
yticks([0 20 40 60 80])
ylabel('Rate (sp/s)')
xlabel('Frequency (Hz)')
grid on
box on
set(gca,'fontsize',14)

% Double Peak (Neuron 1, Model)
subplot(4, 3, [6, 9, 12])
hold on
plot(param{2}.fpeaks, avBS_doublepeak, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410])
plot(fpeaks_SC(3,:), rate_SC(3,:), 'LineWidth', 2);
xlim(freq_limits)
ylim(rate_limits)
yticks([0 20 40 60 80])
yticklabels([])
xlabel('Frequency (Hz)')
grid on
box on
set(gca,'fontsize',14)

%%-------------------------------------------------------------------------

%% Load in BE Cell

figure('Position',[140,224,1327,535])
freq_limits = [400 1600];

base = getPaths();
fpath = 'data/aro-2021/Figure 4-2';
file_names = dir(fullfile(base, fpath, 'R0*'));

for ii = 1:length(file_names)
	data_BE(ii) = load(fullfile(base, fpath, file_names(ii).name));
end

%% Stimuli

% Stimulus parameters
subplot(4, 3, 2)
ylabel('Level (dB SPL)')
ylim([30 65])
xlim(freq_limits)
xticklabels([])
grid on
hold on
box on

% Stimulus Creation
harmonics = param{3}.Delta_F:param{3}.Delta_F:10000; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = param{3}.dur * param{3}.Fs; % # pts in stimulus
t = (0:(npts-1))/param{3}.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/param{3}.fpeak_mid)*param{3}.g)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
	comp_freq = harmonics(iharm);
	component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
	stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(param{3}.spl/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

% Now, make the stimulus for this_fpeak
for iharm = 1:num_harmonics
	stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);

	% Plot each stimulus
	y = fft(stim);
	mdB = 20*log10(abs(y));
	level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
	stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
		2, 'Color', [0.6350, 0.0780, 0.1840]);
end
plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
set(gca,'fontsize',14)

%%  Double Peak Stimulus
subplot(4, 3, 3)

% Stimulus parameters
ylim([30 65])
xlim(freq_limits)
xticklabels([])
yticklabels([])
grid on
hold on
box on

% Stimulus Creation
harmonics = param{4}.Delta_F:param{4}.Delta_F:10000; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = param{4}.dur * param{4}.Fs; % # pts in stimulus
t = (0:(npts-1))/param{4}.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/param{4}.fpeak_mid)*param{4}.g)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
	comp_freq = harmonics(iharm);
	component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
	stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(param{4}.spl/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

% Now, make the stimulus for this_fpeak
for iharm = 1:num_harmonics
	stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);

	% Plot each stimulus
	y = fft(stim);
	mdB = 20*log10(abs(y));
	level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
	stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
		2, 'Color', [0.6350, 0.0780, 0.1840]);
end
plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
set(gca,'fontsize',14)

% RM
[spls,~,si] = unique(double([data_BE(1).param.list.spl]).');
num_spls = length(spls);
[freqs_RM,~,fi] = unique(double([data_BE(1).param.list.freq]).');
num_freqs = length(freqs_RM);
rate_size = [num_freqs,num_spls];
num_spikes = data_BE(1).spike_info.num_spikes_peri;
spike_rates = num_spikes*1000/data_BE(1).param.dur; % spikes/sec
[rates_RM,~,~,~] = accumstats({fi,si},spike_rates, rate_size);
spont = mean(rates_RM(:,1));
max_y = max(rates_RM(:));

% Plot
subplot(4, 3, 4)
hold on
plot(log10(freqs_RM([1 end])),[1 1]*spont,'-','LineWidth',2, 'Color',[0.7 0.7 0.7])
plot(log10(freqs_RM),rates_RM(:,5),'color','#332288', 'LineWidth',2)

% Label Plot
box on
set(gca,'XTick',log10([250 500 1000 2000 5000 10000]), 'XTickLabel',{'0.25','0.5','1','2','5','10'})
xlims = log10(freqs_RM([1 end]));
xlim(xlims.*[0.98;1.02])
ylim([0 max_y])
ylabel('Rate (sp/s)')
xlabel('Frequency at 70dB SPL (kHz)')
grid on
hold off
legend('Spontaneous Rate', 'fontsize',10)
set(gca,'fontsize',14)

% MTF
dur = data_BE(2).param.dur/1000; % stimulus duration in seconds.
all_mod_depths = double([data_BE(2).param.list.mdepth]).';
all_mod_depths(all_mod_depths < -80) = max(all_mod_depths);
[~,~,mdi] = unique(all_mod_depths);
[fms,~,fmi] = unique(double([data_BE(2).param.list.fm]).');
if fms(1) == 0
	fms(1) = 1.2;
end

num_mod_freqs = length(fms);
num_depths = length(data_BE(2).param.all_mdepths);
map_size = [num_mod_freqs, num_depths];
spike_rates = data_BE(2).spike_info.num_spikes_delayed/...
	(dur - data_BE(2).spike_info.onsetWin/1000);
[rate_MTF,~,rlb,rub] = accumstats({fmi,mdi},spike_rates, map_size);
rate_sm = zeros(size(rate_MTF));
for j = 1:num_depths
	rate_sm(:,j) = smooth_rates(rate_MTF(:,j),rlb(:,j),rub(:,j));
end

% Plot
subplot(4, 3, [7 10])
hold on
line([1 fms(end)], [1 1]*rate_MTF(1),'Color',[0.7 0.7 0.7], 'LineWidth', 2);
line(fms,rate_MTF,'Color','#332288', 'Marker','.', 'MarkerSize',10, 'LineWidth', 2);

% Label the plots
box on
grid on
xtick = [1 2 5 10 20 50 100 200 500];
xlim(xtick([1 end]))
xlabel('Modulation Freq (Hz)')
ylabel('Rate (sp/s)')
set(gca,'XTick',xtick,'XScale', 'log')
hold off
legend('No modulation', 'fontsize',10)
set(gca,'fontsize',14)

% Single Peak (Neuron 2, Model)
for ii = 3:4
	[fpeaks_SC(ii, :),~,fpeaksi] = unique([data_BE(ii).param.list.fpeak].');
	num_fpeaks = length(fpeaks_SC);
	dur = data_BE(ii).param.dur/1000; % stimulus duration in seconds.

	rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
	spike_rates = data_BE(ii).spike_info.num_spikes_delayed/...
		(dur - data_BE(ii).spike_info.onsetWin/1000);
	[rate_SC(ii,:),~, rlb, rub] = accumstats({fpeaksi},spike_rates, rate_size);
	for j = 1:num_depths
		rate_SC_sm(ii,:) = smooth_rates(rate_SC(ii,:),rlb,rub);
	end
end

subplot(4,3, [5 8 11])
hold on
plot(param{3}.fpeaks, avBE_singlepeak, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410])
plot(fpeaks_SC(4,:), rate_SC(4,:), 'LineWidth', 2);
xlim(freq_limits)
ylim([0 160])
yticks([0 30 60 90 120 150])
xticklabels([600 800 1000 1200 1400])
ylabel('Rate (sp/s)')
xlabel('Frequency (Hz)')
grid on
box on
set(gca,'fontsize',14)

% Double Peak (Neuron 2, Model)
subplot(4, 3, [6 9 12])
hold on
plot(param{4}.fpeaks, avBE_doublepeak, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410])
plot(fpeaks_SC(3,:), rate_SC(3,:), 'LineWidth', 2);
xlim(freq_limits)
xticklabels([600 800 1000 1200 1400])
ylim([0 160])
yticks([0 30 60 90 120 150])
yticklabels([])
xlabel('Frequency (Hz)')
grid on
box on
set(gca,'fontsize',14)
set(gcf,'color','w');

