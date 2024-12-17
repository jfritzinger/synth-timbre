%% 2021-1-11 ARO Figure Code
% Code that creates the figures for 2021 ARO Winter Meeting
% J. Fritzinger, updated 1/14/2021


%% Fig. 1: Static stimulus

% Figure Properties
figure(1)
%set(gcf,'color','#F2F2F2');
set(gcf,'color','w');
title('Spectral Centroid Stimulus')
xlabel('Frequency (Hz)')
ylabel('Sound Level (dB SPL)')
ylim([0 70])
xlim([200 7000])
grid on
hold on
xticks([200 400 600 1200 2400 4800])
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
plot(harmonics, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
set(gca,'fontsize',14)
%close figure(1)

%% Fig. 2: Stimulus movie

% Figure properties
x_pos = 700;   % Screen position
y_pos = 400;   % Screen position
width  = 1000; % Width of figure
height = 250; % Height of figure (by default in pixels)

% Stimulus parameters
freq_lo = 600;
freq_hi = 1800;
Fs = 100000;
F0 = 200;
Fc = 1200;
dur = 0.3;
rampdur = 0.02;
stimdB = 70;
G = 24;
steps = 25;
%fpeaks = linspace(freq_lo, freq_hi, steps);
fpeaks = [1200:-50:600 1200:50:1800];

harmonics = F0:F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
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

% Stimulus Creation
for i = 1:length(fpeaks)
    
    % Figure Properties
    h = figure('Position', [x_pos y_pos width height]);
    %set(gcf,'color','#F2F2F2');
    set(gcf,'color','w');
    xlabel('Frequency (Hz)')
    ylabel('Level (dB SPL)')
    ylim([0 70])
    xlim([200 2200])
    grid on
    hold on
    xticks(200:200:2200)
    xline(1200, '-', 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410])
    box on
    
    % Now, make the stimulus for this_fpeak
    level = zeros(length(harmonics), 1);
    shift = fpeaks(i) - Fc; % a negative values for low fpeaks; 0 at center; positive for high fpeaks
    comp_freq = harmonics + shift;
    for iharm = 1:num_harmonics
        if comp_freq(iharm) > 75
            stim = component_scales_linear(iharm) * sin(2*pi*comp_freq(iharm)*t);
            
            % Plot each stimulus
            y = fft(stim);
            mdB = 20*log10(abs(y));
            level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
            stem(comp_freq(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
                2, 'Color', '#882255');
        end
    end
    
    % Plot envelope of stimulus
    plot(comp_freq, level, 'LineWidth', 1.5, 'Color', '#882255', 'LineStyle', ':');
    set(gca,'fontsize',16);

    % Write to GIF file
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1
%         imwrite(imind,cm,'Figure2.gif','gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'Figure2.gif','gif','WriteMode','append');
%     end
%     
%     close(h)
end




