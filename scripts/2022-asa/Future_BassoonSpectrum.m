%% Load in names of music files

figure('position', [375,303,700,250]); % left bottom width height)
hold on
ylimits = [0 70];
xlimits = [150 7000];

% Load 
%instrument = 'Oboe.ff.C4.wav'; % original example
instrument = 'Bassoon.ff.Db3.wav';
%instrument = 'Bassoon.ff.E2.wav';
%instrument = 'BbClarinet.ff.A3.wav';
if ismac
    fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Physiology/07 - Natural Timbre/Cut Waveforms/';
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
set(gca,'fontsize',20)
title('Bassoon Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB SPL)')
