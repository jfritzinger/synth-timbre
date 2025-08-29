function stim = generate_spectral_centroid(dur, rampdur, F0, CF, G, stimdB, Fs)
% Generate_single_formant - Lyzenga & Horst's triangular spectrum (from D. Schwarz code)
%    dur = duration in seconds  (0.5 s)
%    ramdur = 0.02 (s)
%    F0 = fundamental freq. in Hertz  (200 Hz)
%    CF = spectral centroid (1200)
%    stimdB = desired SPL of composite waveform (dB SPL) (70 dB SPL)
%    G = slope of triangle edges in dB/octave (24)
%    Fs = sample rate, samples/sec

F_max = 10000; % max frequency
Pref = 20e-6; % 20 micropascals
F_min = F0;

% Calculate intercepts at which signal is at 10000 Hz.
max_rel_level = -G*log2(F_max/CF);
min_rel_level = -G*log2(F_min/CF);
tri_x = log2([F_min CF F_max]);
tri_y = [stimdB - min_rel_level, stimdB, stimdB + max_rel_level];

% Generate time vector.
t = (0:1/Fs:dur)';
t(end) = [];

% Calculate which harmonics we need.
harm_num = ceil(F_min/F0):floor(F_max/F0);
f_harm = harm_num*F0;
log2_harm_freqs = log2(f_harm);
SPLs = interp1(tri_x,tri_y,log2_harm_freqs);
amp = sqrt(2)*Pref*10.^(SPLs/20); % This converts dB SPL to peak amplitude in Pascals
sin_matrix = sin(t*(2*pi*f_harm)); % zero phase

% Multiply each col of sin_matrix by appropriate amplitude, add up the
% columns, and multiply by ramp envelope.
stim = sin_matrix * amp';
stim = Pref * 10.^(stimdB/20) * stim / rms(stim); % scale overall RMS level into desired pascals
stim = stim .* tukeywin(length(stim),2*rampdur/dur); % apply ramp



