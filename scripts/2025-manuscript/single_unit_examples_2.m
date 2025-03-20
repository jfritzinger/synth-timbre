%% single_unit_examples_2.m
%
%
%
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Set up figure

linewidth = 1.5;
figure('Position',[4,295,796,612])
%tiledlayout(3, 3)
data_colors = {'#03882F', '#82BB95'};
legsize = 12;

%% Plot

examples = {'R24_TT2_P13_N05', 'R27_TT2_P8_N02', 'R27_TT2_P8_N05', ...
	'R25_TT2_P9_N01', 'R27_TT3_P1_N08', 'R27_TT2_P7_N01', ...
	'R29_TT4_P5_N15', 'R25_TT2_P8_N02', 'R29_TT1_P3_N05'};

for ineuron = 1:9

	% Load in data
	putative = examples{ineuron};
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
	fpeaks_re_CF = log2(fpeaks/CF);

	h(ineuron) = subplot(3, 3, ineuron);
	yyaxis left
	hold on
	rates_sm = smooth_rates(rate, rlb, rub, CF);
	errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
		'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{1})
	plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{1})
	plot(fpeaks./1000, rates_sm, 'linestyle', '-', 'linewidth', linewidth, 'color', 'k')
	xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	% errorbar(fpeaks_re_CF, rate, rate_std/sqrt(params{1}.nrep), ...
	% 	'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{1})
	% plot(fpeaks_re_CF, rate, 'LineWidth',linewidth, 'Color',data_colors{1})
	% plot(fpeaks_re_CF, rates_sm,'linestyle', '-', 'linewidth', linewidth, 'color', 'k')
	% xline(0, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)

	% Figure parameters
	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
	set(gca, 'Fontsize', 14) %, 'XTick', plot_range(1)+0.200:0.400:plot_range(2)-0.200);
	xlim(plot_range);
	grid on
	ylim([0 max_rate+5])


	if mod(ineuron, 3) == 1
		ylabel('Avg. Rate (sp/s)')
	end

	if ineuron > 6
		xlabel('Spectral Peak Freq. (Hz)')
		%xlabel('Spectral Peak Freq. w.r.t. CF (oct)')
	end

	% Legend
	if ineuron == 9 || ineuron == 3
		% hLeg = legend({'', 'Data', 'Smoothed Data', 'Spont Rate', 'CF'}, ...
		% 	'Location','southwest', 'fontsize', legsize);
		% hLeg.ItemTokenSize = [12, 12];
	end

	%xlim([-1 1])
	if ineuron == 1 || ineuron == 4 || ineuron == 7
		xlim([0.4 2.4])
	elseif ineuron == 2 || ineuron == 5 || ineuron == 8
		xlim([0.8 3.2])
	else
		xlim([2.8 9.2])
	end

	% Labeling MTF Type
	if ineuron == 9
		text(0.05, 0.85, MTF_shape, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',16)
	else
		text(0.05, 0.95, MTF_shape, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',16)
	end
	
	% Parameters
	params = params{1};
	params.Fs = 100000;
	params.dur = 0.3; % duration
	params.mnrep = 1;
	params.stp_otc = 1;
	params.physio = 1;
	params.fpeak_mid = CF;
	params= generate_ST(params);

	% Calculates number of stimuli to be played
	nstim = length(params.fpeaks);

	% Create the Stimulus Gating function
	fs = params.Fs;
	npts = floor(params.dur*fs);
	gate = tukeywin(npts,2*params.ramp_dur/params.dur); %raised cosine ramps

	% Generate stimuli for all presentations
	params.stim = zeros(nstim*params.mnrep, npts);
	presentation = 0; %this value is used as an index for storing a stumulus presentation in the 3rd dimenstion of 'stimuli'
	
	% Compute one stimulus waveform.
	this_fpeak = params.fpeaks(1); % Get peak freq for this presentation

	% Compute fixed set of scalars for central stimulus to obtain spectral envelope & desired stimdB dB SPL
	harmonics = params.Delta_F:params.Delta_F:10000; % component freqs for the central stimulus, when this_fpeak = CF
	num_harmonics = length(harmonics);
	npts = params.dur * fs; % # pts in stimulus
	t = (0:(npts-1))/fs; % time vector
	component_scales_linear = 10.^(-1*abs(log2(harmonics/params.fpeak_mid)*params.g)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
	interval = zeros(1,npts);
	for iharm = 1:num_harmonics
		comp_freq = harmonics(iharm);
		component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
		interval = interval + component;          %Add component to interval
	end
	Level_scale = 20e-6*10.^(params.spl/20) * (1/rms(interval)); % overall lienar scalar to bring this centered stimulus up to stimdB
	component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

	% Time vectors
	npts = params.dur * fs; % # pts in stimulus
	t = (0:(npts-1))/fs; % time vector
	interval = zeros(1,length(t));

	% Make the stimulus for this_fpeak
	shift = this_fpeak - params.fpeak_mid; % a negative values for low fpeaks; 0 at center; positive for high fpeaks
	%figure
	
	for iharm = 1:num_harmonics
		comp_freq = (harmonics(iharm) + shift);
		if comp_freq > 75 % Hz; make sure we don't include comps outside calibrated range (Note: because we'll lop off components, then scale to, say, 70 dB SPL overall - the comp amps will change whenever one component is eliminated.
			interval = interval + component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
			int_single = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
		else
			int_single = 0;
		end
		shifted_harms(iharm) = comp_freq;
		shifted_harms_re_CF(iharm) = log2(comp_freq/CF);

		% Plot each stimulus
		y = fft(int_single);
		mdB = 20*log10(abs(y));
		if sum(y)==0
			level(iharm) = NaN;
		else
			level(iharm) = findpeaks(mdB(1:length(mdB)/2),...
				'MinPeakProminence',200);
		end

		% Plot harmonics
		yyaxis right
		hold on
		stem(shifted_harms(iharm)/1000, level(iharm), '-', 'Marker',...
			'none', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
	end
	% Plot 
	plot(shifted_harms/1000, level, '--', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
	ylim([0 120])
	h(ineuron).YAxis(1).Color = 'k';
	h(ineuron).YAxis(2).Color = 'w';
end



%% Set locations

% Titles
titles_x = {'Low CF', 'Medium CF', 'High CF'};
locs = linspace(0.18, 0.8, 3);
for ii = 1:3
	annotation('textbox',[locs(ii) 0.93 0.4 0.0459],...
		'String',titles_x{ii},...
		'FontSize',18,'EdgeColor','none', ...
		'FontWeight','bold');
end

titles_y = {'Sloping', 'Dip', 'Peak'};
locs = linspace(0.13, 0.81, 3);
for ii = 1:3
	annotation('textbox',[0.06 locs(ii) 0.126 0.0459],...
		'String',titles_y{ii}, ...
		'FontSize',20,'EdgeColor','none', 'Rotation',90, ...
		'FontWeight','bold');
end

left = repmat(linspace(0.12, 0.72, 3), 1, 3);
bottom = repmat(linspace(0.1, 0.73, 3), 3, 1);
bottom = fliplr(reshape(bottom, 1, 9));
width = 0.24;
height = 0.20;

for ii = 1:9
	set(h(ii), 'Position', [left(ii) bottom(ii) width height])
end

% Set annotations
labelsize = 24;
annotation('textbox',[0.01 0.95 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.63 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.33 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

%% Export figure

%exportgraphics(gcf, fullfile(savepath, 'manuscript', 'examples-peaksanddips.png'), 'Resolution', 600)
