%% R024S470, 8/20/2021, TT2N1, Cell: Excellent, CF = 1170Hz, Depth 7.75
close all 

% Load excel file
rabbit = 'R024';
session = 'R024S470_TT2N1.mat';
base = getPaths();
fpath = 'data/asa-2022';


datafile = fullfile(base, fpath, rabbit, session);
data_full = load(datafile, 'saved_data');
data_full = data_full.saved_data;
has_sc = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
    strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
    d.param.fpeak_mid==1000,data_full);

load(fullfile(base, fpath, "Fig3_BSModel.mat"))
load(fullfile(base, fpath, "Fig3_Data.mat"))

%% Figure Properties 

data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};

ylimits = [0 80];
fontname = 'Arial'; 
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
figure;
set(gcf, 'Position', [460,427,1300,420])
set(gcf,'color','w');
tiledlayout(1, 2, 'TileSpacing','compact')

%% Calculate explainable variance explained by the model 

for ind = 1:size(BS.rates, 1)
    R_int = corrcoef(BS.rates(ind,:),avBS(ind,:)).^2;
    R(ind) = R_int(1, 2);
    
    [hat_r2er_BS(ind), r2_BS] = r2er_n2m2(avBS(ind, :)*1.7, BS.rate_matrix{ind});
    
    disp(['DSID = ' num2str(BS.dsid{ind}) ', SC = ' ...
        num2str(BS.fpeak_mid{ind}) ', ' num2str(BS.spl{ind}) ' dB SPL'])
    disp(['Variance explained = ' num2str(R(ind)) ', Explainable variance explained = ' num2str(hat_r2er_BS(ind))]);
end

%% Plot Model
nexttile 
hold on
for ind = 1:4
    plot(BS.fpeaks{1},avBS(ind, :)*1.7, 'linewidth', 1.5, 'Color', model_colors{ind});
end
xlabel('Spectral Peak Frequency (Hz)')
ylim(ylimits)
set(gca, 'XTick', [BS.plot_range(1)+200:400:BS.plot_range(2)-200]);
xlim(BS.plot_range);
plot(BS.fpeaks{1}([1 end]),[1 1]*spontBS,'-','Color',[0.5 0.5 0.5], 'LineWidth', 2)
xline(BS.CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);

for ii = 1:length(BS.spl)
    leg(ii,:) = [num2str(BS.spl{ii}+3) ' dB SPL'];
end
legend(char(leg,['Spont. Rate'], ['CF']),...
    'location', 'northeast', 'FontSize',18)

set(gca,'FontSize',20)
box on
grid on
title('Model', 'FontSize',22)
ylabel('Avg. Rate (sp/s)')

%% Plot data
ylimits = [0 80];

nexttile
hold on
for ind = 1:4
    errorbar(BS.fpeaks{ind},BS.rates(ind,:),BS.rates_std(ind,:)/(sqrt(BS.nrep{ind})),...
        'LineWidth', 1, 'Color', data_colors{ind});
end
plot(BS.fpeaks{1}([1 end]),[1 1]*BS.spont,'-','Color',[0.5 0.5 0.5], 'LineWidth', 2)
xline(BS.CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);

for ii = 1:length(BS.spl)
    leg(ii,:) = [num2str(BS.spl{ii}+3) ' dB SPL'];
end
legend(char(leg,['Spont. Rate'], ['CF']),...
    'location', 'northeast', 'FontSize',18)

xlabel('Spectral Peak Frequency (Hz)')
xlim(BS.plot_range);
set(gca, 'XTick', [BS.plot_range(1)+200:400:BS.plot_range(2)-200]);
ylabel('Avg. Rate (sp/s)')
ylim(ylimits)
set(gca,'FontSize',20)
box on
grid on
title('Single-Unit', 'FontSize',22);

%% Plot MTF 

has_mtf = cellfun(@(d)strcmp(d.param.type,'typMTFN'), data_full);
data = data_full(has_mtf);
ind = 1;

figure;
set(gcf, 'Position', [560,527,600,400])
set(gcf,'color','w');

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
set(gca,'fontsize',20)
legend('No modulation', 'fontsize',19, 'location', 'northwest')
title('Band-Suppressed', 'FontSize', 24);


%% RM 

figure;
set(gcf, 'Position', [560,527,800,400])
set(gcf,'color','w');
hold on

% Plot RM
plot(RM.freqs,RM.rates(:,5),'color', '#20116B','LineWidth',2) % 70 dB
plot(RM.freqs,RM.rates(:,4),'color', '#5E50A9','LineWidth',2) % 50 dB 7B6DC1
plot(RM.freqs,RM.rates(:,3),'color', '#A49BD0','LineWidth',2) % 30 dB
plot(RM.freqs([1 end]),[1 1]*RM.spont,'-','LineWidth',2, 'Color',[0.7 0.7 0.7])
xline(BS.CF, '--', 'Color', [0.4 0.4 0.4],'LineWidth',3);
box on
grid on
hold off
%ylim([0 RM.max_y+10])
%set(gca,'XTick',[])
%xticks(x_label)
%yticks([0 100 200])
%xticklabels([])
set(gca,'fontsize',20)
title('Response Map', 'FontSize', 24)
set(gca,'XTick',[250 500 1000 2000 5000 10000])
xlim([250 10000])
legend('70 dB SPL', '50 dB SPL', '30 dB SPL', 'Spont. Rate', 'CF', 'fontsize',19)
ylabel('Avg. Rate (sp/s)')
xlabel('Tone Frequency (Hz)')
set(gca, 'XScale', 'log');
