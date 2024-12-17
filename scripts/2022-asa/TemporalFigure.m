%% Create stimulus 

% Stimulus parameters
% params.Fs = 100000;
% params.this_fpeak = 1200;
% params.F0 = 200;
% params.Fc = 1200;
% params.dur = 0.3;
% params.ramp_dur = 0.02;
% params.fpeak_mid = 1200;
% params.num_harms = 13;
% params.mnrep = 1;
% params.physio = 0;
% params.spl = 70;
% params.g = 24;
% params.steps = 1;
% stim = generate_spectralcentroid(params);
% BMF = 100;
% onset_num = 1;
% stimulus = stim.stim;
% dur = 0.3;
% Fs = 100000;

% Natural timbre
instrument = 'Bassoon.ff.Db3.wav';
if ismac
    fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/03 - Physiology/07 - Natural Timbre/Cut Waveforms/';
else
    fpath = 'C:\Users\jfritzinger\Box\03 - Physiology\07 - Natural Timbre\Cut Waveforms\';
end
[stimulus, Fs] = audioread([fpath instrument]);
dur = 0.3;

%% Plot temporal stimulus
figure;
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(gcf, 'position', [560,469,250,190]); % left bottom width height

npts = dur * Fs; % # pts in stimulus
t = (0:(npts-1))./Fs; % time vector
%plot(t, stimulus, 'Color', '#882255', 'linewidth', 1.5);
plot(t, stimulus, 'Color', [0, 0.4470, 0.7410, 0.5], 'linewidth', 1.5); 
xlim([0.1 0.2])
xlabel('Time (ms)')
xticks([0.1 0.150 0.20])
xticklabels([100 150 200])
yticklabels([])
ylabel('Amplitude')
grid on
set(gca,'FontSize',16)