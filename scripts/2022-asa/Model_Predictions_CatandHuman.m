%% Replica of Maxwell et al., 2021 Average Rate and Temporal Plots
% This script plots the same figure that Maxwell et al., 2021 uses in a
% proceedings paper, but with modified colors and ranges to fit with the
% seminar presentation.
% J. Fritzinger, updated 2/23/2021
clear all 

% Figure Properties
figure;
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(gcf,'color','w');
set(gcf, 'position', [560,112,800,835]); % left bottom width height

%% Plot Model Responsese

% Model parameters
nrep = 30;       % repetitions averaged for profile plots
species = 1:2;   % human & cat
implnt = 0;      % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
noiseType = 1;   % 0 for fixed fGn (1 for variable fGn) - this is the 'noise' associated with spontaneous activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
fiberType = 1;   % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
cohc = 1;
cihc = 1;

% Stimulus parameters
params.Fs = 100000;
params.this_fpeak = 1200;
params.F0 = 200;
params.Fc = 1200;
params.dur = 0.3;
params.ramp_dur = 0.02;
params.stimdB = 70;
params.G = 24;
params.steps = 1;
stimulus = generate_spectralcentroid(params);
BMF = 100;
onset_num = 1;

% Figure parameters
CF_list = [1000 1200 1400];
start_time = 0.1; % s - for plotting
start  = start_time * params.Fs;
stop_time = 0.115; % s
stop = stop_time * params.Fs;
b_LP = fir1(5000,100/(params.Fs/2),'low'); % LP at 100 Hz, for envelope plots
set(gcf,'color','w');

%% Model and Plotting
species_order = [1];
for species_ind = 1
    iplot = 0; % step thru CFs within each panel
    inoise = 1;
    for CF_plot = CF_list  % for TIN response to 1000 Hz tone
        iplot = iplot + 1;
        for irep = 1:nrep % do nrep repetitions - plot the last one, but use all 10 for averages
            impaired = 0; % not impaired!
            excit = gammatone_filter_simple(stimulus,CF_plot,params.Fs,impaired);
            RMS_DV(irep) = rms(excit); %mean(an_sout(start:stop))

            % AN model responses
            vihc = model_IHC(stimulus',CF_plot,nrep,1/params.Fs,params.dur*1.2,cohc,cihc,2);
            [an_sout,~,~] = model_Synapse(vihc,CF_plot,nrep,1/params.Fs,fiberType,noiseType,implnt); % an_sout is the auditory-nerve synapse output - a rate vs. time function that could be used to drive a spike generator
            an_sout_LP = conv(an_sout,b_LP,'Same'); % smooth AN PSTH to remove much of fine-structure, so it doesn't dominate a simple "enevelope slope" DV
            env_an_sout = envelope(an_sout_LP);
            env_slope_DV(irep) = mean(abs(diff(env_an_sout)*params.Fs))/mean(an_sout); %mean(abs(diff(env_an_sout(start:stop))*Fs))/mean(an_sout(start:stop));
            rate_DV(irep) = mean(an_sout); % mean(an_sout(start:stop));
            %                 RMS_DV(irep) = mean(an_sout(start:stop));

        end % last rep will be plotted below, but profiles will be based on means

        %% Plot AN Temporal
        t = (1:length(an_sout))/params.Fs; % time vector for plots

        h(iplot) = subplot(5, 5, 11+iplot);
        hold on
        set(gca,'fontsize',14, 'xtick',[100 125 150], 'YDir', 'reverse')

        % Duplicate the AN responses for Fluctuation profile plot
        if species_ind == 1
            plot(an_sout(start:stop), t(start:stop)*1e3,'Color', '#44AA99','linewidth',1.5) % blue
            [peaks_vals,peak_locations] = findpeaks(an_sout(start:stop),params.Fs,'MinPeakProminence',100);
            plot(peaks_vals, (peak_locations + start_time)*1e3,'Color', '#117733','linewidth',3) % add envelope
        else
            plot(an_sout(start:stop), t(start:stop)*1e3,'Color', '#D07F5C','linewidth',1.5) % orange
            [peaks_vals,peak_locations] = findpeaks(an_sout(start:stop),params.Fs,'MinPeakProminence',100);
            plot(peaks_vals, (peak_locations + start_time)*1e3,'Color', '#D95319','linewidth',3) % add envelope
        end

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
        set(gca,'fontsize',17)
        xticks([0 1000])
        box on
        grid on

    end
end
set(h(1), 'Position',[0.23 0.45 0.12 0.12])
set(h(2), 'Position',[0.43 0.45 0.12 0.12])
set(h(3), 'Position',[0.63 0.45 0.12 0.12])

%% Plot stimulus

% Stimulus parameters
params.Fs = 100000;
params.this_fpeak = 1200;
params.F0 = 200;
params.Fc = 1200;
params.dur = 0.3;
params.ramp_dur = 0.02;
params.stimdB = 80;
params.G = 24;
params.steps = 1;
params.num_harms = 13;
[params.freq_lo, params.freq_hi] = get_freq_limits(params.Fc, params.num_harms, params.F0);
stimulus = generate_spectralcentroid(params);

% Model parameters
model_params.BMF = 100;
model_params.species = 1;
model_params.CF = params.Fc;

for ii = 1:2
    model_params.species = 1;
    [avAN(ii, :), avBE(ii, :), avBS(ii, :), ~, CFs] = modelPopulation(stimulus', model_params, params, 0, 1);
    %model_params.species = 2;
    %[avAN_cat(ii, :), avBE_cat(ii, :), avBS_cat(ii, :), ~, CFs_cat] = modelPopulation(stimulus', model_params, params, 0, 1);
end
avAN = mean(avAN);
avBE = mean(avBE);
avBS = mean(avBS);
%avAN_cat = mean(avAN_cat);
%avBE_cat = mean(avBE_cat);
%avBS_cat = mean(avBS_cat);


plot_range = [300 7000];
set(gcf,'color','w');
ticks = [200 500 1000 2000 5000];

% Stimulus Creation
subplot(5,5,1:5)
set(gca, 'Position',[0.080 0.80 0.775 0.124])
hold on
harmonics = params.F0:params.F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = params.dur * params.Fs; % # pts in stimulus
t = (0:(npts-1))/params.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/params.Fc)*params.G)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
    comp_freq = harmonics(iharm);
    component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
    stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(params.stimdB/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

% Now, make the stimulus for this_fpeak
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
set(gca,'FontSize',17)

%% AN Plot
subplot(5, 5, 6:10)
set(gca, 'Position',[0.08 0.63 0.775 0.124])
plot(CFs, avAN, 'linewidth', 2, 'color', '#117733'); %blue
hold on
%plot(CFs_cat, avAN_cat, 'linewidth', 2, 'color', '#D95319'); % orange
grid on
set(gca, 'XScale', 'log')
box on
xlim(plot_range)
%ylim([0 300])
%ylim([0 70])
xticks(ticks)
%xticklabels([])
ylabel('Avg. Rate (sp/s)')
xlabel('CF. (Hz)')
set(gca,'FontSize',17)

%% IC BE Plot
subplot(5, 5, 16:20)
set(gca, 'Position',[0.080 0.280 .775 0.124])
hold on
plot(CFs, avBE, 'linewidth', 2, 'color', '#D95319'); % blue
%plot(CFs_cat, avBE_cat, 'linewidth', 2, 'color', '#0072BD'); % orange
grid on
set(gca, 'XScale', 'log')
box on
xlim(plot_range)
%ylim([0 60])
%ylim([0 3])
ylabel('Avg. Rate (sp/s)')
xticks(ticks)
%yticks([0 30 60])
%xticklabels([])
xlabel('CF. (Hz)')
set(gca,'FontSize',17)

%% IC BS Plot
subplot(5, 5, 21:25)
set(gca, 'Position',[0.080 0.10 0.775 0.124])
hold on
%plot(CFs_cat, avBS_cat, 'linewidth', 2, 'color', '#0072BD'); % orange
plot(CFs, avBS, 'linewidth', 2, 'color', '#D95319'); % blue
grid on
set(gca, 'XScale', 'log')
box on
xlim(plot_range)
%ylim([0 60])
%ylim([0 3])
xticks(ticks)
ylabel('Avg. Rate (sp/s)')
%yticks([0 30 60])
xlabel('CF. (Hz)')
legend('Human Tuning', 'Cat Tuning', 'location','northeast')
set(gca,'FontSize',17)

%% Plot temporal stimulus
figure;
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(gcf,'color','w');
set(gcf, 'position', [560,469,250,190]); % left bottom width height

npts = params.dur * params.Fs; % # pts in stimulus
t = (0:(npts-1))/params.Fs; % time vector
plot(t, stimulus, 'Color', '#882255', 'linewidth', 1.5);
xlim([0.1 0.120])
xlabel('Time (ms)')
xticks([0.1 0.110 0.120])
xticklabels([100 110 120])
yticklabels([])
ylabel('Amplitude')
grid on
set(gca,'FontSize',20)
