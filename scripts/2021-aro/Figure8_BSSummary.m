%% Figure 8: Summary of many datasets
% J. Fritzinger, updated 2/23/2021
clear all

%% Model

% Fs = 100000;
% F0 = 200;
% Fc = 1000;
% CF = 1000;
% dur = 0.3;
% rampdur = 0.02;
% stimdB = 70;
% G = 24;
% this_fpeak = 1200;
% species = 1;
% num_harms = 9;
% steps = 41;
% 
% [freq_lo, freq_hi] = get_freq_limits(Fc, num_harms, F0);
% [stim, stim_list, fpeaks_model] = specslide_spectralcentroid(Fs, Fc, F0, dur, rampdur, steps, stimdB, G, freq_lo, freq_hi);
% [~, ~, ~, ~, avBS, ~, ~] = modelSingleCell(stim, stim_list, fpeaks_model, Fc, steps, CF, species);
% 
% % Normalize BS response 
% avBS_norm = avBS/max(avBS);

%% Load in datafiles

base = getPaths();
fpath = 'data/aro-2021/Figure 8';
file_names = dir(fullfile(base, fpath, 'R0*'));
%file_namesRM = dir([fpath '\RM\R*']);

data(1) = load(fullfile(base, fpath, 'R024S294_TT2N1_ds7.mat')); 
data(2) = load(fullfile(base, fpath, 'R025S293_TT3N4_ds11.mat')); 
data(3) = load(fullfile(base, fpath, 'R025S275_TT4N1_ds17.mat')); 
data(4) = load(fullfile(base, fpath, 'R025S363_TT2N1_ds10.mat')); 
data(5) = load(fullfile(base, fpath, 'R025S368_TT2N1_ds10.mat')); 
data(6) = load(fullfile(base, fpath, 'R025S372_TT2N2_ds9.mat'));  
data(7) = load(fullfile(base, fpath, 'R024S295_TT2N1_ds9.mat'));
data(8) = load(fullfile(base, fpath, 'R024S298_TT2N1_ds9.mat')); 
data(9) = load(fullfile(base, fpath, 'R024S303_TT2N1_ds8.mat'));
data(10) = load(fullfile(base, fpath, 'R025S375_TT2N1_ds9.mat'));
data(11) = load(fullfile(base, fpath, 'R025S377_TT2N1_ds9.mat')); 
data(12) = load(fullfile(base, fpath, 'R024S325_TT2N1_ds7.mat')); 
data(13) = load(fullfile(base, fpath, 'R024S347_TT2N1_ds8.mat'));

units = {'Multi-unit', 'Single-unit', 'Multi-unit', 'Multi-unit', ...
    'Multi-unit', 'Single-unit', 'Multi-unit', 'Multi-unit', 'Multi-unit', ...
    'Multi-unit', 'Multi-unit', 'Multi-unit', 'Multi-unit'};

%% Analyze

% Get spont rate
% for ii = 1:length(dataRM)
%     [spls,~,si] = unique(double([dataRM(ii).param.list.spl]).');
%     num_spls = length(spls);
%     [freqs,~,fi] = unique(double([dataRM(ii).param.list.freq]).');
%     num_freqs = length(freqs);
%     
%     rate_size = [num_freqs,num_spls];
%     num_spikes = dataRM(ii).spike_info.num_spikes_peri;
%     spike_rates = num_spikes*1000/dataRM(ii).param.dur; % spikes/sec
%     [rates,~,~,~] = accumstats({fi,si},spike_rates, rate_size);
%     
%     spont(ii) = mean(rates(:,1));
% end

for ii = 1:length(data)
    
    %label{ii} = ['R0' num2str(data(ii).param.animal) 'S' num2str(data(ii).param.session) ', ' units{ii} ', CF = ' num2str(data(ii).param.fpeak_mid) 'Hz'];
    label{ii} = [units{ii} ', CF = ' num2str(data(ii).param.fpeak_mid) 'Hz'];

    % Find the different peaks (xaxis)
    [fpeaks{ii},~,fpeaksi] = unique([data(ii).param.list.fpeak].');
    num_fpeaks = length(fpeaks{ii});
    dur = data(ii).param.dur/1000; % stimulus duration in seconds.
    
    % Calculate average spike rate (yaxis)
    rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
    spike_rates = data(ii).spike_info.num_spikes_delayed/...
        (dur - data(ii).spike_info.onsetWin/1000);
    [rate{ii},rate_std] = accumstats({fpeaksi},spike_rates, rate_size);
    
%     % Normalize spike rate
%     rate_norm(ii,:) = rate(ii,:) - spont(ii);
%     max_y(ii) = max(rate_norm(ii,:));
%     rate_norm(ii,:) = rate_norm(ii,:)/max_y(ii);
%     
%     % Center the CF
%     shift = 1000 - data(ii).param.fpeak_mid;
%     fpeaks_norm(ii,:) = fpeaks(ii,:) + shift;  
end


%% Plot raw

% Color array
colors = {'#E2BA01', '#D68102', '#D68102', '#D68102', '#D68102', '#D68102', ...
    '#C35702', '#C35702', '#C35702', '#C35702', '#8A3E01', '#8A3E01', ...
    '#632C01', '#3E1C01'};

% h = tiledlayout(1,2);
% h.TileSpacing = 'compact';
% h.Padding = 'compact';
% nexttile
figure;
%subplot(1,2,1)
hold on
for ii = 1:length(data)
    if ii == 2 || ii == 6
        plot(fpeaks{ii}, rate{ii}, '-.', 'color', colors{ii}, 'linewidth', 2.5);
    else
        plot(fpeaks{ii}, rate{ii}, 'color', colors{ii}, 'linewidth', 2.5);
    end
end

xlabel('Frequency (Hz)')
xlim([400 3200])
grid on
box on
title('Band-Suppressed IC Neuron Responses')
ylabel('Rate (spikes/sec)')
set(gca,'fontsize',14)
hold off
set(gcf,'color','w')
set(gca, 'xscale', 'log')
legend(label, 'fontsize',16, 'location', 'northeastoutside')

%% Plot normalized and plot model 

% % Color array
% colors = {'#ECC262', '#E68940', '#BD5503', '#7D3700', ...
%     '#ECC262', '#E68940', '#BD5503', '#7D3700', ...
%     '#ECC262', '#E68940', '#BD5503', '#7D3700'};
% 
% nexttile
% hold on
% for ii = 1:length(data)
%     plot(fpeaks_norm(ii,:), rate_norm(ii,:), 'color', colors{ii}, 'linewidth', 2);
% end
% 
% % Plot model
% plot(fpeaks_model, avBS_norm, 'linewidth', 2, 'Color', '#2F66A9');
% 
% % Parameters 
% xticks([200 400 600 800 1000 1200 1400 1600 1800])
% xticklabels({'4F0', '3F0', '2F0', 'F0', 'CF', 'F0', '2F0', '3F0', '4F0'})
% xlabel('Frequency (Hz)')
% xlim([600 1800])
% grid on
% box on
% title('Normalized IC Response')
% ylabel('Normalized Rate')
% set(gca,'fontsize',14)
% hold off
% set(gcf,'color','w');
