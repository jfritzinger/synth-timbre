function plotModelResponses(trial, CFs, avAN, avBE, avBS, stimulus, Fs, ...
    Which_IC, cond, Fc, F0, Fc_order, F0_order)
% Plots model responses
% J. Fritzinger, updated 8/20/2020
%
% Inputs: trial: first or second trial
%         avAN: average rate AN response
%         avBE: average rate BE response
%         avBS: average rate BS response
%         stimulus: spectral centroid input
%         Fs: sampling rate
%         Which_IC: SFIE or AN
%         cond: which condition is being played
%         Fc: spectral centroid (Hz)
%         F0: fundamental frequency (Hz)
%         Fc_order: which Fc is being played first
%         F0_order: which F0 is being played first 

if trial == 1
    figure;
end

if cond == 1 % timbre
    order = Fc_order;
else
    order = F0_order;
end

% Parameters
xlabels = [200 500 1000 2000 3000 5000];
ylimit = [0 300];
plot_range = [125 5000];

%% Plot average AN model response
h(1) = subplot(4, 2, order(trial));
hold on
semilogx(CFs, avAN) % log axes
xlim(plot_range);
ylim(ylimit) % fix rate (sp/sec) axis to allow comparisons
set(gca, 'XScale', 'log', 'XTick', xlabels);
grid on
xlabel('AN CF (Hz)','fontsize',10)
ylabel('Ave Rate (sp/sec)','fontsize',10)
title('AN Response Profile','fontsize',10)

%% Plot average IC response
if Which_IC == 1  % SFIE model - plot toggles between BE and BS model responses
    
    % Plots BE response
    h(2) = subplot(4, 2, 2 + order(trial));
    hold on
    semilogx(CFs,avBE)
    title('Midbrain BE Response Profile')
    xlabel('IC BF (kHz)')
    ylabel('Ave Rate (sp/sec)')
    set(gca, 'XScale', 'log', 'XTick', xlabels);
    grid on
    xlim([plot_range(1) plot_range(2)]);
    ylim([0 50]);
    
    % Plots BS response
    h(3) = subplot(4, 2, 4 + order(trial));
    hold on
    semilogx(CFs, avBS)
    
%     y_plt = supsmu(CFs, average_ic_sout_BS, 'Span', 0.5);
%     y_done = average_ic_sout_BS-y_plt;
%     y_done(y_done<0)=0;
%     plot(CFs, y_plt);
%     plot(CFs, y_done);
%     
%     log_CFs = log10(CFs);
%     smoothing = supsmu(log_CFs, average_ic_sout_BS, 'Span', 0.5);
%     y_d = average_ic_sout_BS - smoothing;
%     y_d(y_d<0)=0;
%     weighted_avg2 = 10^((sum(log_CFs.*y_d))/(sum(y_d)));

    log_CFs = log10(CFs);
    mean_sub_BS = max(avBS - 28, 0);
    weighted_avg = 10^((sum(log_CFs.*mean_sub_BS))/(sum(mean_sub_BS)));
    
    plot(CFs, mean_sub_BS);
    plot(weighted_avg, 0, 'X', 'Color', 'k');
    
%     weighted_avg_plot = (sum(CFs.*y_done))/(sum(y_done));
%     plot(weighted_avg_plot, 0, 'X', 'Color', 'k');
%     plot(weighted_avg2, 0, 'X', 'Color', 'r');
    
    %legend('IC Response', 'Smoothed Response', 'Rectified IC-Smoothed Response', 'Calculated Spectral Centroid')
    title(['Midbrain BS Response Profile, Est Fc = ' num2str(round(weighted_avg))])
    xlabel('IC CF (Hz)')
    ylabel('Ave Rate (sp/sec)')
    set(gca, 'XScale', 'log', 'XTick', xlabels);
    grid on
    xlim([plot_range(1) plot_range(2)]);
    ylim([0 50]);
    
elseif  Which_IC == 2  % AN_ModFilt model  (models only Band-Enhanced cells)
    
    % Plots BE from ModFilt model
    h(2) = subplot(4, 2, 2 + order(trial));
    hold on
    semilogx(CFs, avBE)
    xlabel('IC BF (Hz)','fontsize',18)
    ylabel('Ave Rate (sp/sec)','fontsize',18)
    xlim([plot_range(1) plot_range(2)]);
    set(gca, 'XScale', 'log', 'XTick',[200 500 1000 2000 3000 5000]);
    grid on
    ylim([0 50]);
end

%% Plot stimuli
h(4) = subplot(4, 2, 6 + order(trial));
new_fs = (2.^9)/0.005;
stimulus_new = resample(stimulus(trial,:), new_fs, Fs);
m = length(stimulus_new);
nfft = pow2(nextpow2(m));  % Find next power of 2
spectrum_plot = 20*log10(abs(2*fft(stimulus_new,nfft)/m/20e-6)); % see: http://12000.org/my_notes/on_scaling_factor_for_ftt_in_matlab/
fres = new_fs/nfft;
freq = fres*(0:nfft-1);
% specplot_max = -inf; % intialize
% specplot_max = max(specplot_max, max(20*log10(abs(fft(stimulus_new,nfft)/nfft/20e-6))));
semilogx(freq,spectrum_plot);
xlim([plot_range(1) plot_range(2)])
ylim([0, 70])
set(gca, 'XScale', 'log', 'XTick',[200 500 1000 2000 3000 5000]);
grid on
box off
title(['Stimulus, Fc = ' num2str(round(Fc(Fc_order(trial)))) ,', F0 = ' num2str(round(F0(F0_order(trial))))])
ylabel('Magnitude (dB SPL)')
xlabel('Frequency (Hz)')

linkaxes(h,'x');



end