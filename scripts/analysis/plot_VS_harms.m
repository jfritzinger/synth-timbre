%% analyze_phase_locking_to_harms
clear 

%% Load data 

[base, datapath, savepath, ppi] = getPaths();
load(fullfile(datapath, 'VS_harms_data.mat'), 'VS')

%% Plot 

figure;
tiledlayout(2, 2)

nexttile
scatter(VS.CF, VS.max_VS_all);
xlabel('CF (Hz)');
ylabel('Phase Locking Value');
title('Max VS');

nexttile
scatter(VS.CF, VS.max_VS_harm_all);
xlabel('CF (Hz)');
ylabel('Phase Locking Value');
title('Max VS Harmonic');

nexttile
scatter(VS.CF, VS.avg_VS_all);
xlabel('CF (Hz)');
ylabel('Phase Locking Value');
title('Avg VS');

nexttile
scatter(VS.CF, VS.avg_VS_harm_all);
xlabel('CF (Hz)');
ylabel('Phase Locking Value');
title('Avg VS Harmonic');