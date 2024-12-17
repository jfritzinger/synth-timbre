%% Fig. 7: Human vs cat tuning
% Plots a population model comparing human results vs cat results 
% J. Fritzinger, updated 2/23/2021

%%  Run models 

% Stimulus parameters
param.fpeak_mid = 1200;
param.Delta_F = 200;
param.dur = 0.3;
param.ramp_dur = 0.02;
param.spl = 70;
param.g = 24;
param.num_harms = 9;
param.stp_otc = 1;
param.Fs = 100000;
param.mnrep = 1;
param.physio = 0;
param = generate_ST(param);
param.num_stim = size(param.stim, 1);

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.BMF = 100;
model_params.num_CFs = 30;
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 3; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
model_params.CF_range = linspace(400, 2000, 30);
model_params.CFs = linspace(400, 2000, 30);

% Run model
for istim = 1:2
	model_params.species = istim; % 1 = cat, 2 = human
	AN = modelAN(param, model_params); % HSR for IC input
	SFIE = wrapperIC(AN.an_sout, param, model_params); % SFIE output

	if istim == 1
		av_BScat = SFIE.average_ic_sout_BS;
	else
		av_BShuman = SFIE.average_ic_sout_BS;
	end
end

%% Plot 

% Figure Properties
h = tiledlayout(1,2);
h.TileSpacing = 'compact';
h.Padding = 'compact';
set(gcf,'color','w');
yy_limit = [20 140];
y_limit = [0 50];
xlimits = [400 2000];
left_color = [0 0 0];
right_color = [0 0 0];
set(h,'defaultAxesColorOrder',[left_color; right_color]);

% Stimulus Creation
harmonics = param.Delta_F:param.Delta_F:10000; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = param.dur * param.Fs; % # pts in stimulus
t = (0:(npts-1))/param.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/param.fpeak_mid)*param.g)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
    comp_freq = harmonics(iharm);
    component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
    stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(param.spl/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

%--------------------------------------------------------------------------
% Human model (BS)
nexttile
plot(model_params.CFs, av_BShuman, 'LineWidth', 2, 'Color', '#009E73')
xlabel('IC CF (Hz)')
ylabel('Rate (sp/sec)')
ylim(y_limit)
xlim(xlimits)
xticks([400 600 800 1000 1200 1400 1600 1800 2000])
xticklabels({'400' '' '800' '' '1200' '' '1600' '' '2000'})
grid on
hold on
xlabel('Frequency (Hz)')
title('Human Tuning')
ylabel('Avg. Spikes')

yyaxis right
ylim(yy_limit)
yticklabels([])

level = zeros(num_harmonics,1);
for iharm = 1:num_harmonics
    stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);
    
    % Plot each stimulus
    y = fft(stim);
    mdB = 20*log10(abs(y));
    level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
    stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
        2, 'Color', [0.6350, 0.0780, 0.1840], 'LineStyle', '-');
end
set(gca,'fontsize',16)

%--------------------------------------------------------------------------
% Cat model (BS)

nexttile
plot(model_params.CFs, av_BScat, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410])
xlabel('IC CF (Hz)')
xlim(xlimits)
xticks([400 600 800 1000 1200 1400 1600 1800 2000])
xticklabels({'400' '' '800' '' '1200' '' '1600' '' '2000'})
grid on
hold on
xlabel('Frequency (Hz)')
title('Cat Tuning')

ylim(y_limit)
yticklabels([])

yyaxis right
ylim(yy_limit)
ylabel('Sound Level (dB SPL)')
for iharm = 1:num_harmonics
    stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);
    
    % Plot each stimulus
    y = fft(stim);
    mdB = 20*log10(abs(y));
    level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
    stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
        2, 'Color', [0.6350, 0.0780, 0.1840], 'LineStyle', '-');
end
set(gca,'fontsize',14)
