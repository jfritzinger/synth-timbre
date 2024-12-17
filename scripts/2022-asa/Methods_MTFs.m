%% MTF
clear all

% Figure Properties
figure('position', [560,573,550,450]); % left bottom width height
hybrid = 'R025S512_TT4N1';
BE = 'R025S474_TT1N2'; % contra, excellent 
BS = 'R024S482_TT4N2'; % bin, excellent 
rabbit = 'R024';

base = getPaths();
fpath = 'data/asa-2022';

%% BS MTF 

datafile = fullfile(base, fpath, rabbit, [BS '.mat']);
data = load(datafile, 'saved_data');
data = data.saved_data;
has_mtf = cellfun(@(d)strcmp(d.param.type,'typMTFN'), data);
ind = 5;

% Calculate MTF
dur = data{ind}.param.dur/1000; % stimulus duration in seconds.
all_mod_depths = double([data{ind}.param.list.mdepth]).';
all_mod_depths(all_mod_depths < -80) = max(all_mod_depths);
[~,~,mdi] = unique(all_mod_depths);
[fms,~,fmi] = unique(double([data{ind}.param.list.fm]).');
if fms(1) == 0
    fms(1) = 1.2;
end

num_mod_freqs = length(fms);
num_depths = length(data{ind}.param.all_mdepths);
map_size = [num_mod_freqs, num_depths];

spike_rates = data{ind}.spike_info.num_spikes_delayed/...
    (dur - data{ind}.spike_info.onsetWin/1000);
[rate_MTF,rate_std,rlb,rub] = accumstats({fmi,mdi},spike_rates, map_size);
rate_sm = zeros(size(rate_MTF));
for j = 1:num_depths
    rate_sm(:,j) = smooth_rates(rate_MTF(:,j),rlb(:,j),rub(:,j));
end

% Plot
hold on
line([1 fms(end)], [1 1]*rate_MTF(1),'Color',[0.7 0.7 0.7], 'LineWidth', 2);
errorbar(fms,rate_MTF,rate_std,'Color','k', 'Marker','.', 'MarkerSize',10, 'LineWidth', 2);
plot(fms,rate_sm,'-b', 'LineWidth', 2)

% Label the plots
box on
grid on
xtick = [1 2 5 10 20 50 100 200 500];
xlim(xtick([1 end]))
xlabel('Modulation Freq (Hz)')
ylabel('Spike Rate (sp/s)')
set(gca,'XTick',xtick,'XScale', 'log')
hold off
legend('No modulation', 'fontsize',22, 'location', 'northwest')
set(gca,'fontsize',22)
set(gcf, 'color', 'w')
title('Modulation Transfer Function', 'FontSize', 28)
