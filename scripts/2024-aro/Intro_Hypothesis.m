%% Intro_Hypothesis
% J. Fritzinger, updated 1/16/24 
clear 

%% Set up figure
figure('Position',[552,305,693,1005])
font_size = 16;

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

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CF_range = [125 10000];
model_params.num_CFs = 100;
model_params.CFs = logspace(log10(125), log10(10000), 100);
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 10; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
model_params.Fs = 100000;


%% Model 

%AN_HSR = modelAN(params, model_params); % HSR for IC input
%SFIE = wrapperIC(AN_HSR.an_sout, params, model_params); % SFIE output
%save('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Aim 2 - Timbre/Data/Intro_ModelResponse.mat', 'AN_HSR', 'SFIE')

base = getPaths();
fpath = 'data/aro-2024';
filename = 'Intro_ModelResponse.mat';
load(fullfile(base, fpath, filename), 'AN_HSR', 'SFIE')
avAN = AN_HSR.average_AN_sout;
avBE = SFIE.average_ic_sout_BE;
avBS = SFIE.average_ic_sout_BS;
CFs = AN_HSR.CFs;

%% Plot temporal 

% Figure parameters
CF_list = [875 1200 1550];
start_time = 0.1; % s - for plotting
start  = start_time * params.Fs;
stop_time = 0.115; % s
stop = stop_time * params.Fs;
b_LP = fir1(5000,100/(params.Fs/2),'low'); % LP at 100 Hz, for envelope plots

iplot = 0; % step thru CFs within each panel

for CF_plot = CF_list  % for TIN response to 1000 Hz tone
    iplot = iplot + 1;

	% AN model responses
	vihc = model_IHC(params.stim,CF_plot,model_params.nrep,1/params.Fs,...
		params.dur*1.2,model_params.cohc,model_params.cihc,model_params.species);
	[an_sout,~,~] = model_Synapse(vihc,CF_plot,model_params.nrep,1/params.Fs,...
		model_params.fiberType,model_params.noiseType,model_params.implnt); % an_sout is the auditory-nerve synapse output - a rate vs. time function that could be used to drive a spike generator
	an_sout_LP = conv(an_sout,b_LP,'Same'); % smooth AN PSTH to remove much of fine-structure, so it doesn't dominate a simple "enevelope slope" DV
	env_an_sout = envelope(an_sout_LP);
    
    % Plot AN Temporal
    t = (1:length(an_sout))/params.Fs; % time vector for plots
    
    h(iplot) = subplot(5, 5, 11+iplot);
    hold on
    set(gca,'fontsize',14, 'xtick',[100 125 150], 'YDir', 'reverse')
    
    % Duplicate the AN responses for Fluctuation profile plot
    plot(an_sout(start:stop), t(start:stop)*1e3,'Color', '#44AA99','linewidth',1.5)
    [peaks_vals,peak_locations] = findpeaks(an_sout(start:stop),params.Fs,'MinPeakProminence',100);
    plot(peaks_vals, (peak_locations + start_time)*1e3,'Color', '#117733','linewidth',3) % add envelope
    
    ylim([start_time*1e3, stop_time*1e3])
    xlim([0 1000]); % freq
    
    if iplot == 1
        ylabel('Time (ms)')
        yticks([100 115])
        yticklabels([100 115])
    elseif iplot == 2
        yticks([100 115])
        yticklabels([])
        xlabel('Rate (sp/s)')
    else
        yticks([100 115])
        yticklabels([])
    end
    set(gca,'fontsize',font_size)
    xticks([0 1000])
    box on
    grid on
    
end

%% Plot stimulus 

plot_range = [300 6000];
ticks = [200 500 1000 2000 5000];

% Stimulus Creation
h(4) = subplot(5,5,1:5);
hold on
harmonics = params.Delta_F:params.Delta_F:10000; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = params.dur * params.Fs; % # pts in stimulus
t = (0:(npts-1))/params.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/params.fpeak_mid)*params.g)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
    comp_freq = harmonics(iharm);
    component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
    stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(params.spl/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

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
box on
hold off
xlim(plot_range)
xticks(ticks)
%xticklabels([])
xlabel('Freq. (Hz)')

%% AN Plot
h(5) = subplot(5, 5, 6:10);
hold on
plot(CFs, avAN, 'linewidth', 2, 'color', '#117733');
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
box on
xlim(plot_range)
ylim([0 300])
xticks(ticks)
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
xlabel('CF (Hz)')

%xticklabels([])
ylabel('Avg. Rate (sp/s)')

%% IC BE Plot
h(6) = subplot(5, 5, 17:20);
hold on
plot(CFs, avBE, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
box on
xlim(plot_range)
ylim([0 60])
%ylabel('Avg. Rate (sp/s)                         ')
xticks(ticks)
yticks([0 30 60])
xticklabels([])
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line

%% IC BS Plot
h(7) = subplot(5, 5, 22:25);
hold on
plot(CFs, avBS, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
box on
xlim(plot_range)
ylim([0 60])
xticks(ticks)
%ylabel('Avg. Rate (sp/s)')
yticks([0 30 60])
xlabel('CF (Hz)')
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line

%% Labels / Annotations / Positions 

% Create textbox
annotation('textbox',...
	[0.030366492146597 0.834755813953488 0.1 0.0450581395348837],...
	'String',{'Stimulus'},...
	'FontSize',18,...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0 0.622093023255814 0.156020942408377 0.102197674418604],...
	'String',{'AN Rate','Response'},...
	'HorizontalAlignment','right',...
	'FontSize',18,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.0209424083769641 0.136453488372093 0.112041884816754 0.0716744186046504],...
	'String',{'IC BS Rate','Response'},...
	'HorizontalAlignment','right',...
	'FontSize',18,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.0199895287958121 0.32343023255814 0.112041884816754 0.0716744186046505],...
	'String',{'IC BE Rate Response'},...
	'HorizontalAlignment','right',...
	'FontSize',18,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.0147068062827229 0.494011627906977 0.131937172774869 0.0600465116279062],...
	'String',{'AN Temporal Response'},...
	'HorizontalAlignment','right',...
	'FontSize',18,...
	'FitBoxToText','off',...
	'EdgeColor','none');

%%
% Create arrow
annotation('arrow',[0.634 0.634],[0.69 0.579], 'linewidth', 1.5); % top 
annotation('arrow',[0.634 0.634],[0.4103 0.306], 'linewidth', 1.5); % bottom

% Create arrow
annotation('arrow',[0.565 0.5],[0.69 0.579], 'linewidth', 1.5); % top 
annotation('arrow',[0.5 0.565],[0.436 0.33], 'linewidth', 1.5); % bottom

% Create arrow
annotation('arrow',[0.69 0.8],[0.694 0.579], 'linewidth', 1.5);
annotation('arrow',[0.8 0.69],[0.433 0.32], 'linewidth', 1.5);


%% Plot MTFs
clear params

% % Stimulus Parameters and Generation
% params.type = 'typMTFN';
% params.ramp_dur = 0.05;
% params.noise_state = 0;
% params.noise_band = [100, 10000];
% params.dur = 1; %1; % s
% params.reptim = 1.5;
% params.fms = [2, 600, 3]; % fm_lo, fm_hi, steps per octave
% params.mdepths = [0,0,1];
% params.binmode = 2;
% params.No = 30;
% params.spl = 30;
% params.all_mdepths = 0;
% params.Fs = 100000;
% params.nrep = 1;
% params.mnrep = 5;
% params.raised_sine = 1;
% params = generate_MTF(params);
% params.num_stim = size(params.stim, 1);
% 
% % Model Parameters 
% model_params.type = 'SFIE';
% model_params.range = 2; % 1 = population model, 2 = single cell model
% model_params.species = 1; % 1 = cat, 2 = human
% model_params.BMF = 100;
% model_params.CF_range = 1200;
% model_params.num_CFs = 1;
% model_params.CFs = 1200;
% model_params.nAN_fibers_per_CF = 10;
% model_params.cohc = 1; % (0-1 where 1 is normal)
% model_params.cihc = 1; % (0-1 where 1 is normal)
% model_params.nrep = 1; % how many times to run the AN model
% model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
% model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
% model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
% model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
% model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
% model_params.Fs = 100000;
% 
% % Run models 
% AN_HSR = modelAN(params, model_params); % HSR for IC input
% SFIE = wrapperIC(AN_HSR.an_sout, params, model_params); % SFIE output
% save('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Aim 2 - Timbre/Data/Intro_MTFModelResponse.mat', 'AN_HSR', 'SFIE', 'params')

filename = 'Intro_MTFModelResponse.mat';
load(fullfile(base, fpath, filename), 'AN_HSR', 'SFIE', 'params')

% Plot BE
h(8) = subplot(5, 5, 16);
[~, avBE, stdBE] = plotMTF(params, SFIE.average_ic_sout_BE, 0);
hold on
line([1 params.all_fms(end)], [1 1]*avBE(1),'Color',[0.4 0.4 0.4], 'LineWidth', 2);
%errorbar(params.all_fms,avBE_SFIE, stdBE_SFIE,'.');
%line(params.all_fms,avBE_SFIE,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w');
plot(params.all_fms,smooth(avBE),'-k', 'LineWidth', 2) % smoothed MTF
hold off
xtick = [1 2 5 10 20 50 100 200 500];
xlim(xtick([1 end]))
ylabel('Avg. Rate (sp/s)                         ')
set(gca,'XTick',xtick,'XScale', 'log')
hLegend = legend('Unmod.', 'Location','northwest', 'EdgeColor','none');
hLegend.ItemTokenSize = [12,8];
grid on
set(gca,'FontSize',font_size)
xticklabels([])
yticklabels([])
ylim([22 32])
box on

% Plot BS
h(9) = subplot(5, 5, 21);
[~, avBS, stdBS] = plotMTF(params, SFIE.average_ic_sout_BS, 0);
hold on
line([1 params.all_fms(end)], [1 1]*avBS(1),'Color',[0.4 0.4 0.4], 'LineWidth', 2);
%errorbar(params.all_fms,avBS_SFIE, stdBS_SFIE,'.');
%line(params.all_fms,avBS_SFIE,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w');
plot(params.all_fms,smooth(avBS),'-k', 'LineWidth', 2) % smoothed MTF
hold off
xtick = [1 2 5 10 20 50 100 200 500];
xlim(xtick([1 end]))
ylim([0 34])
set(gca,'XTick',xtick,'XScale', 'log')
xticklabels([])
yticklabels([])
grid on
set(gca,'FontSize',font_size)
box on
ylim([19 27])
xlabel('Mod. Freq (Hz)')

%% Plot temporal stimulus 
% 
% figure;
% set(gcf, 'color', 'w')
% npts = params.dur * params.Fs; % # pts in stimulus
% t = (0:(npts-1))/params.Fs; % time vector
% plot(t, stimulus, 'Color', '#882255', 'linewidth', 1.5);
% xlim([0.1 0.120])
% xlabel('Time (ms)')
% xticks([0.1 0.110 0.120])
% xticklabels([100 110 120])
% set(gca,'fontsize',12)
% ylabel('Amplitude')
% grid on


%% Move positions
set(h(1), 'Position',[0.41 0.448 0.13 0.124]);
set(h(2), 'Position',[0.58 0.448 0.13 0.124]);
set(h(3), 'Position',[0.75 0.448 0.13 0.124]);
set(h(4), 'Position',[0.33 0.797 0.655 0.124]);
set(h(5), 'Position',[0.33 0.616 0.655 0.124]);
set(h(6), 'Position',[0.33 0.272 0.655 0.124]);
set(h(7), 'Position',[0.33 0.101 0.655 0.124]);
set(h(8), 'Position',[0.17 0.272 0.123 0.124]);
set(h(9), 'Position',[0.17 0.101 0.123 0.124]);