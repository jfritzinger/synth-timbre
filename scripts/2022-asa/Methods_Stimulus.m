%%

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

figure
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(gcf,'color','w');
set(gcf, 'position', [560,112,650,650]); % left bottom width height

% Stimulus Creation
data = [3, 1, 16];
for i = 1:3

    % Figure Properties
    % Figure Properties
   
    %set(gcf,'color','#F2F2F2');
    h(i) = subplot(3, 1, i);
    if i == 3
        xlabel('Frequency (Hz)')
    elseif i == 2
        ylabel('Level (dB SPL)')
    end
    ylim([0 70])
    yticks([0 20 40 60])
    xlim([200 2200])
    grid on
    hold on
    xticks(200:200:2200)
    xticklabels([])
    xline(1200, '-', 'LineWidth', 4, 'Color', [0, 0.4470, 0.7410, 0.7])
    box on
    if i ~= 3
        xlabel([])
        xticklabels([])
    end

    % Now, make the stimulus for this_fpeak
    level = zeros(length(harmonics), 1);
    shift = fpeaks(data(i)) - Fc; % a negative values for low fpeaks; 0 at center; positive for high fpeaks
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
    plot(comp_freq, level, 'LineWidth', 1.5, 'Color',  [136/256, 34/256, 85/256, 0.5], 'LineStyle', '--');
    set(gca,'fontsize',22);

end

% Create textbox
annotation(gcf,'textbox',...
    [0.49 0.91 0.095 0.057],...
    'Color',[0 0.450980392156863 0.741176470588235],...
    'String',{'CF'},...
    'FontSize',26,...
    'EdgeColor','none');

set(h(1), 'position', [0.15,0.64,0.75,0.27]); % left bottom width height
set(h(2), 'position', [0.15,0.37,0.75,0.27]);
set(h(3), 'position', [0.15,0.1,0.75,0.27]);
