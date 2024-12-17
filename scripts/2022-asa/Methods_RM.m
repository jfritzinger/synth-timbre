%% Response map 

%% Plot

% Figure
figure('position', [560,573,550,450])
hold on

name = 'R024S485_TT2N1'; % excellent 
rabbit = 'R024';
base = getPaths();
fpath = 'data/asa-2022';


datafile = fullfile(base, fpath, rabbit, [name '.mat']);
data = load(datafile, 'saved_data');
data = data.saved_data;
has_rm = cellfun(@(d)strcmp(d.param.type,'type=RM'), data);
ind = 4;

% Calculate data
[spls,~,si] = unique(double([data{ind}.param.list.spl]).');
num_spls = length(spls);
[freqs_RM,~,fi] = unique(double([data{ind}.param.list.freq]).');
num_freqs = length(freqs_RM);

rate_size = [num_freqs,num_spls];
num_spikes = data{ind}.spike_info.num_spikes_peri;
spike_rates = num_spikes*1000/data{ind}.param.dur; % spikes/sec
[rates_RM,~,~,~] = accumstats({fi,si},spike_rates, rate_size);

spont = mean(rates_RM(:,1));
max_y = max(rates_RM(:));

% Plot RM
CF = 1500;
plot(freqs_RM,rates_RM(:,5),'color', '#20116B','LineWidth',2) % 70 dB
plot(freqs_RM,rates_RM(:,4),'color', '#5E50A9','LineWidth',2) % 50 dB 7B6DC1
plot(freqs_RM,rates_RM(:,3),'color', '#A49BD0','LineWidth',2) % 30 dB
plot(freqs_RM([1 end]),[1 1]*spont,'-','LineWidth',2, 'Color',[0.7 0.7 0.7])
xline(CF, '--', 'Color', [0.4 0.4 0.4],'LineWidth',3);
box on
grid on
hold off
%ylim([0 RM.max_y+10])
%set(gca,'XTick',[])
%xticks(x_label)
%yticks([0 100 200])
%xticklabels([])
set(gca,'fontsize',22)
title('Response Map', 'FontSize', 28)
set(gca,'XTick',[250 500 1000 2000 5000 10000])
xlim([250 10000])
legend('70 dB SPL', '50 dB SPL', '30 dB SPL', 'Spont. Rate', 'CF', 'fontsize',19)
ylabel('Avg. Rate (sp/s)')
xlabel('Tone Frequency (Hz)')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
set(gcf, 'color', 'w')