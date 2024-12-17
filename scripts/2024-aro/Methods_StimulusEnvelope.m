%% Script that creates synthetic timbre envelope 
% J. Fritzinger, updated 1/17/24 
clear

% Parameters
params.Fs = 100000;
params.F0 = 200; % fundamental frequency 
params.Fc = 1400; % peak harmonic frequency 
params.dur = 0.3; % duration
params.ramp_dur = 0.02; % ramp duration
params.stimdB = 70; % overall level
params.G = 24; % slope, in dB/oct 
%params.fpeaks = [850 1100 1400 1700 2050];
params.fpeaks = [900 1400 2250];
params.mnrep = 1;

% Calculates number of stimuli to be played, just for diplaying the time
% requirement estimation
nstim = length(params.fpeaks);

% Create the Stimulus Gating function
fs = params.Fs;
npts = floor(params.dur*fs);
gate = tukeywin(npts,2*params.ramp_dur/params.dur); %raised cosine ramps

% Generate stimuli for all presentations
params.stim = zeros(nstim*params.mnrep, npts);
presentation = 0; %this value is used as an index for storing a stumulus presentation in the 3rd dimenstion of 'stimuli'

figure('Position',[228,592,670,385]);
hold on
line([params.Fc params.Fc], [0 100], 'LineWidth', 3,'LineStyle', ':', 'Color', [.7 .7 .7]);

% Create stimuli for each rep (irep) and each stimulus (istim)
for istim = 1:nstim
	presentation = presentation + 1;

	% Compute one stimulus waveform.
	this_fpeak = params.fpeaks(istim); % Get peak freq for this presentation

	% Compute fixed set of scalars for central stimulus to obtain spectral envelope & desired stimdB dB SPL
	harmonics = params.F0:params.F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
	num_harmonics = length(harmonics);
	npts = params.dur * fs; % # pts in stimulus
	t = (0:(npts-1))/fs; % time vector
	component_scales_linear = 10.^(-1*abs(log2(harmonics/params.Fc)*params.G)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
	interval = zeros(1,npts);
	for iharm = 1:num_harmonics
		comp_freq = harmonics(iharm);
		component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
		interval = interval + component;          %Add component to interval
	end
	Level_scale = 20e-6*10.^(params.stimdB/20) * (1/rms(interval)); % overall lienar scalar to bring this centered stimulus up to stimdB
	component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

	% Time vectors
	npts = params.dur * fs; % # pts in stimulus
	t = (0:(npts-1))/fs; % time vector
	interval = zeros(1,length(t));
	harmonics = params.F0:params.F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
	num_harmonics = length(harmonics);

	% Make the stimulus for this_fpeak
	shift = this_fpeak - params.Fc; % a negative values for low fpeaks; 0 at center; positive for high fpeaks
	for iharm = 1:num_harmonics
		comp_freq = (harmonics(iharm) + shift);
		if comp_freq > 75 % Hz; make sure we don't include comps outside calibrated range (Note: because we'll lop off components, then scale to, say, 70 dB SPL overall - the comp amps will change whenever one component is eliminated.
			interval = interval + component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
			int_single = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
		else
			int_single = 0;
		end
		shifted_harms(iharm) = comp_freq;

		% Plot each stimulus
		y = fft(int_single);
		mdB = 20*log10(abs(y));
		if sum(y)==0
			level(iharm) = NaN;
		else
			level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
		end

		% Plot harmonics
		if istim == 2
			stem(shifted_harms(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
				2, 'Color', '#882255');
		else
			stem(shifted_harms(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
			2, 'Color', '#DEA6C2');
		end
	end

	if istim == 2
		plot(shifted_harms, level, 'LineWidth', 3, 'Color', '#882255');
	else
		plot(shifted_harms, level, 'LineWidth', 3, 'Color', '#DEA6C2');
	end
	ylim([0 70])
	ylabel('Magnitude (dB SPL)')
	xlabel('Frequency (Hz)')
	
	xlim([180 3500])
	xticks(params.Fc)
	xticklabels({'CF'})
	box on
end
set(gca, 'FontSize', 18)
title('Synthetic Timbre Stimulus', 'FontSize',20)

% Create arrow
annotation('arrow',[0.45 0.60],[0.851 0.851]);

% Create arrow
annotation('arrow',[0.38 0.22],[0.851 0.851]);

%% Log 
% clear
% 
% params.Fs = 100000;
% params.F0 = 200; % fundamental frequency 
% params.Fc = 1400; % peak harmonic frequency 
% params.dur = 0.3; % duration
% params.ramp_dur = 0.02; % ramp duration
% params.stimdB = 70; % overall level
% params.G = 24; % slope, in dB/oct 
% %params.fpeaks = [850 1100 1400 1700 2050];
% params.fpeaks = [900 1400 2300];
% params.mnrep = 1;
% 
% % Calculates number of stimuli to be played, just for diplaying the time
% % requirement estimation
% nstim = length(params.fpeaks);
% fs = params.Fs;
% npts = floor(params.dur*fs);
% 
% % Generate stimuli for all presentations
% params.stim = zeros(nstim*params.mnrep, npts);
% presentation = 0; %this value is used as an index for storing a stumulus presentation in the 3rd dimenstion of 'stimuli'
% 
% figure('Position',[228,771,746,351]);
% hold on
% line([params.Fc params.Fc], [0 100], 'LineWidth', 3,'LineStyle', ':', 'Color', [.7 .7 .7]);
% 
% % Create stimuli for each rep (irep) and each stimulus (istim)
% for istim = 1:nstim
% 	presentation = presentation + 1;
% 
% 	% Compute one stimulus waveform.
% 	this_fpeak = params.fpeaks(istim); % Get peak freq for this presentation
% 
% 	% Compute fixed set of scalars for central stimulus to obtain spectral envelope & desired stimdB dB SPL
% 	harmonics = params.F0:params.F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
% 	num_harmonics = length(harmonics);
% 	npts = params.dur * fs; % # pts in stimulus
% 	t = (0:(npts-1))/fs; % time vector
% 	component_scales_linear = 10.^(-1*abs(log2(harmonics/params.Fc)*params.G)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
% 	interval = zeros(1,npts);
% 	for iharm = 1:num_harmonics
% 		comp_freq = harmonics(iharm);
% 		component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
% 		interval = interval + component;          %Add component to interval
% 	end
% 	Level_scale = 20e-6*10.^(params.stimdB/20) * (1/rms(interval)); % overall lienar scalar to bring this centered stimulus up to stimdB
% 	component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars
% 
% 	% Time vectors
% 	npts = params.dur * fs; % # pts in stimulus
% 	t = (0:(npts-1))/fs; % time vector
% 	interval = zeros(1,length(t));
% 	harmonics = params.F0:params.F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
% 	num_harmonics = length(harmonics);
% 
% 	% Make the stimulus for this_fpeak
% 	shift = this_fpeak - params.Fc; % a negative values for low fpeaks; 0 at center; positive for high fpeaks
% 	for iharm = 1:num_harmonics
% 		comp_freq = (harmonics(iharm) + shift);
% 		if comp_freq > 75 % Hz; make sure we don't include comps outside calibrated range (Note: because we'll lop off components, then scale to, say, 70 dB SPL overall - the comp amps will change whenever one component is eliminated.
% 			interval = interval + component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
% 			int_single = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
% 		else
% 			int_single = 0;
% 		end
% 		shifted_harms(iharm) = comp_freq;
% 
% 		% Plot each harmonic
% 		y = fft(int_single);
% 		mdB = 20*log10(abs(y));
% 		if sum(y)==0
% 			level(iharm) = NaN;
% 		else
% 			level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
% 		end
% 
% 		%%%% CARTOON PLOT: Here I plot the individual harmonics as stem
% 		%%%% plots scaled to a particular level %%%%%%%
% 		% Plot harmonics
% 		if istim == 2
% 			stem(shifted_harms(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
% 				1, 'Color', '#882255');
% 		else
% 			stem(shifted_harms(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
% 			1, 'Color', '#DEA6C2');
% 		end
% 	end
% 
% 	%%%% CARTOON PLOT: Here I plot the envelope of the stimulus %%%%% 
% 	if istim == 2
% 		plot(shifted_harms, level, 'LineWidth', 1.5, 'Color', '#882255'); % center harmonic matched to CF
% 	else
% 		plot(shifted_harms, level, 'LineWidth', 1.5, 'Color', '#DEA6C2'); % slid components
% 	end
% 	ylim([0 70])
% 	ylabel('Magnitude (dB SPL)')
% 	xlabel('Frequency')
% 	title('Synthetic Timbre')
% 	set(gca, 'XScale', 'log')
% 	xlim([180 11000])
% 	xticks(params.Fc)
% 	xticklabels({'CF'})
% 	box on
% end
% 
% set(gca, 'FontSize', 16)
% 
% % Create shifting arrows
% annotation('arrow',[0.537 0.650],[0.851 0.851]);
% annotation('arrow',[0.5 0.39],[0.851 0.851]);
% 
% % Create F0 label
% annotation('line',[0.224 0.224],[0.740 0.685]);
% annotation('line',[0.321, 0.321],[0.737 0.688]);
% annotation('line',[0.224 0.321],[0.713 0.713]);
% annotation('textbox',[0.23 0.735 0.15 0.066],...
% 	'String',{'F0 = 200Hz'},...
% 	'FontSize', 14, ...
% 	'EdgeColor','none');


%%


% % Create textbox
% annotation('textbox',...
% 	[0.664233576642336 0.812235294117647 0.161105318039625 0.0661764705882353],...
% 	'String',{'Spectrum & harmonics shifted'},...
% 	'EdgeColor','none');

% 
% % Set parameters of synthetic timbre 
% Fs = 100000;
% F0 = 200; % fundamental frequency 
% Fc = 1400; % peak harmonic frequency 
% dur = 0.3; % duration
% rampdur = 0.02; % ramp duration
% stimdB = 70; % overall level
% G = 24; % slope, in dB/oct 
% 
% figure;
% hold on
% 
% nstim = 3;
% fpeaks = [900 1400 1900];
% for istim = 1:nstim
% 
% 	% Compute one stimulus waveform.
% 	this_fpeak = fpeaks(istim);
% 
% 	% Generate synthetic timbre stimulus
% 	harmonics = F0:F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
% 	num_harmonics = length(harmonics);
% 	npts = dur * Fs; % # pts in stimulus
% 	t = (0:(npts-1))/Fs; % time vector
% 	component_scales_linear = 10.^(-1*abs(log2(harmonics/this_fpeak)*G)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
% 	stimulus = zeros(1,npts);
% 	for iharm = 1:num_harmonics
% 		comp_freq = harmonics(iharm);
% 		component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
% 		stimulus = stimulus + component;          %Add component to interval
% 	end
% 	Level_scale = 20e-6*10.^(stimdB/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
% 	component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars
% 
% 	% Calculate the level of each harmonic
% 
% 	line([Fc Fc], [0 100], 'LineWidth', 3,'LineStyle', ':', 'Color', [.7 .7 .7]);
% 	for iharm = 1:num_harmonics
% 		stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);
% 
% 		% Plot each stimulus
% 		y = fft(stim);
% 		mdB = 20*log10(abs(y));
% 		level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
% 
% 		% Plot harmonics
% 		if istim == 2
% 			stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
% 				1, 'Color', '#882255');
% 		else
% 			stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
% 			1, 'Color', '#C7779F');
% 		end
% 	end
% 
% 	% Plot envelope
% 	if istim == 2
% 		plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255');
% 	else
% 		plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#C7779F');
% 	end
% 	ylim([0 70])
% 	ylabel('Magnitude (dB SPL)')
% 	xlabel('Frequency')
% 	title('Spectral Envelope for Synthetic Timbre')
% 	set(gca, 'XScale', 'log')
% 	xlim([180 11000])
% 	xticks(Fc)
% 	xticklabels({'CF'})
% 
% end
