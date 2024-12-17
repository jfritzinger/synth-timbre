%% FutureDirections_NatTimbre
% J. Fritzinger, 1/16/24
clear

%% Make bassoon spectrum

figure('position', [142,513,950,645]); % left bottom width height)
%tiledlayout(4, 1, 'Padding','compact')


%nexttile
h(1) = subplot(4, 1, 1);
hold on
ylimits = [0 70];
xlimits = [150 7000];

% Load
%instrument = 'Oboe.ff.C4.wav'; % original example
instrument = 'Bassoon.ff.Db3.wav';
%instrument = 'Bassoon.ff.E2.wav';
%instrument = 'BbClarinet.ff.A3.wav';
if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/03 - Physiology/07 - Natural Timbre/Cut Waveforms/';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Physiology\07 - Natural Timbre\Cut Waveforms\';
end
[stim, Fs] = audioread([fpath instrument]);

% Calculate spectra
dist = 40;
y2 = fft(stim);
m = abs(y2);
mdB = 20*log10(m);
f = (0:length(y2)-1)*Fs/length(y2);
mdB(mdB<0) = 0;
f(f>Fs/2) = [];
mdB = mdB(1:length(f))';
[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
freqs = f(locs);
plot(freqs, pks, 'LineWidth', 1.5, 'LineStyle', '--', 'Color',...
	[0, 0.4470, 0.7410, 0.5]);

% Plot
for ii = 1:length(pks)
	stem(freqs(ii), pks(ii), 'Marker', 'none', 'LineWidth', ...
		2, 'Color', [0, 0.4470, 0.7410]);
end

% Figure Properties
ylim(ylimits)
%yticks(labels_y)
%yticklabels([])
xlim(xlimits)
xticks([200 400 800 1600 3200 6400])
set(gca, 'XScale', 'log')
hold on
box on
set(gca,'fontsize',18)
title('Bassoon Spectrum Example')
xlabel('Frequency (Hz)')
ylabel('Mag. (dB SPL)')


for ineuron = 1:3
	%% Plot physio & model

	switch ineuron
		case 1 % R025S568, TT2N1, Excellent, CF = 1309Hz, BS - saved 
			session = 'R025S568_TT2_N1';
			CF = 1309;
			MTF_shape = 'BS';
		case 2 % R025S554, TT2N1, Good, CF = 1981Hz, BS	- added, not saved 
			session = 'R024S450_TT2_N1';
			CF = 1718;
			MTF_shape = 'BS';
			ylimits = [0 110];
		case 3 % R024S450, TT2N1, Good, CF = 1718Hz, BS - added, not saved
			session = 'R025S554_TT2_N1';
			CF = 1981;
			MTF_shape = 'BS';
			ylimits = [0 130];
		case 4 % R025S539, TT4N1, Good, CF = 1320Hz, BS - added, not saved
			session = 'R025S539_TT4_N1';
			CF = 1320;
			MTF_shape = 'BS';
			ylimits = [0 70];
		case 5 % R027S015, TT3N1, Good, CF = 1000Hz, BS - added, not saved
			session = 'R027S015_TT3_N1';
			CF = 1000;
			MTF_shape = 'BS';
	end

	%% Load in examples

	base = getPaths();
	fpath = 'data/aro-2024';
	filename = fullfile(base, fpath, [session '.mat']);
	load(filename, 'params', 'data', 'cluster', 'stims')

	% Natural Timbre
	has_nt = cellfun(@(p)strcmp(p.type,'Natural_Timbre'),params);
	data_NT = data(has_nt); % Gets binaural WB-TIN stimuli
	param = params(has_nt);

	[param, ~, data_NT] = plotPhysNT(cluster, param, stims, [], data_NT);
	close

	data_NT = data_NT{1};
	param = param{1};


	%% Run model
	% %Parameters
	% param.Fs = 100000;
	% param.files = param.filename;
	% param.mnrep = 20;
	% 
	% % Generate stimuli
	% [param] = generate_NT(param);
	% param.num_stim = size(param.stim, 1);
	% 
	% % Model parameters
	% model_params.BMF = 100;
	% model_params.species = 1;
	% model_params.CF = CF;
	% model_params.CF_range = CF;
	% model_params.num_CFs = 1;
	% model_params.CFs = CF;
	% model_params.nAN_fibers_per_CF = 10;
	% model_params.cohc = 1; % (0-1 where 1 is normal)
	% model_params.cihc = 1; % (0-1 where 1 is normal)
	% model_params.nrep = 3; % how many times to run the AN model
	% model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
	% model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
	% model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
	% model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
	% model_params.type = MTF_shape;
	% model_params.fiberType = 3;
	% 
	% 
	% % Run model
	% AN = modelAN(param, model_params); % HSR for IC input
	% SFIE = wrapperIC(AN.an_sout, param, model_params); % SFIE output
	% 
	% 
	% % Save 
	% filename = [session '_NT_Model.mat'];
	% save(fullfile(base, fpath, filename), 'AN', 'SFIE', 'param')

	%% Load & Analysis

	filename = [session '_NT_Model.mat'];
	load(fullfile(base, fpath, filename), 'AN', 'SFIE', 'param')

	% Analysis
	[avBS, avBS_std] = plotNT(param, SFIE.avIC, 0);
	R_int = corrcoef(data_NT.rate,avBS).^2;
	R = R_int(1, 2);

	% Energy model 
	Fs = param.Fs;
	stimulus = [param.stim zeros(size(param.stim,1),0.1*Fs)];
	gamma_param.srate = Fs;
	tvals = (1:length(stimulus))/Fs;
	gamma_IF_reg = zeros(length(CF),length(tvals));
	impaired = 0; % 0 = not impaired; 1 = 'impaired'
	pin_gamma = zeros(size(stimulus, 1), Fs*param.dur+0.1*Fs);
	for istim = 1:size(stimulus, 1)
		gamma_param.fc = CF;
		pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired);
	end
	pin_gamma = pin_gamma(:,1:param.dur*Fs);
	energy = sqrt(mean(pin_gamma.^2,2));
	[energy_rate, ~] = plotNT(param, energy, 0);
	R_int = corrcoef(data_NT.rate,energy_rate);
	R_eng = R_int(1, 2).^2;


	%% Plot

	pitch = data_NT.pitch_num;
	h(ineuron+1) = subplot(4, 1, ineuron+1);
	plot(pitch, data_NT.rate,'-*', 'LineWidth',2);
	hold on
	if ineuron == 3
		plot(pitch, avBS*2,'-o', 'LineWidth',2)
	else
		plot(pitch, avBS,'-o', 'LineWidth',2)
	end
	if ineuron == 1
		plot(pitch,energy_rate*2000, 'linewidth', 2);
	else
		plot(pitch,energy_rate*10000, 'linewidth', 2);
	end

	set(gca, 'XScale', 'log')
	box on
	ylim(ylimits)
	xlim([57 575])
	xticks([110 220 440 880])

	if ineuron == 3
		xlabel('Fundamental Frequency (Hz)')
	elseif ineuron == 2
		ylabel('Avg. Rate (sp/s)')
	else
		xticklabels([])
	end
	
	legend(['Data, CF=' num2str(CF) 'Hz'], ['SFIE Model, R^2=' num2str(round(R, 3))], ...
		['Energy, R^2=' num2str(round(R_eng, 3))], ...
		'Location','northwest', 'EdgeColor','none', 'NumColumns',2)
	set(gca,'FontSize',18)
	title('')
end

%% Rearrange 

set(h(1), 'Position',[0.114 0.792 0.85 0.15]);
set(h(2), 'Position',[0.114 0.491 0.85 0.19]);
set(h(3), 'Position',[0.114 0.295 0.85 0.19]);
set(h(4), 'Position',[0.114 0.099 0.85 0.19]);

