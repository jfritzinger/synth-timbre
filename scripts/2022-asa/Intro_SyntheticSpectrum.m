% Methods Figure

% Figure Properties
figure;
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(gcf,'color','w');
set(gcf, 'position', [560,573,724,374]); % left bottom width height


title('Synthetic Timbre Stimulus')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB SPL)')
ylim([0 70])
xlim([200 7000])
grid on
hold on
xticks([200 600 1200 2400 4800])
yticks([0 20 40 60])
%xticks(1200)
%xticklabels('CF')
set(gca, 'XScale', 'log')
box on

% Stimulus parameters
Fs = 100000;
F0 = 200;
Fc = 1200;
dur = 0.3;
rampdur = 0.02;
stimdB = 70;
G = 24;

% Stimulus Creation
harmonics = [F0:F0:10000]; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);
npts = dur * Fs; % # pts in stimulus
t = (0:(npts-1))/Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/Fc)*G)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
    comp_freq = harmonics(iharm);
    component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
    stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(stimdB/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars

% Now, make the stimulus for this_fpeak
for iharm = 1:num_harmonics
    stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);
    
    % Plot each stimulus
    y = fft(stim);
    mdB = 20*log10(abs(y));
    level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
    stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
        2, 'Color', '#882255');
end

% Plot envelope of stimulus
%plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
set(gca,'fontsize',25)