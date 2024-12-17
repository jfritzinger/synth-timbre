% Future directions 
clear all
close all

% Figure Properties
figure('position', [375,303,700,300]); % left bottom width height)

%% Plot

base = getPaths();
fpath = 'data/asa-2022';

load(fullfile(base, fpath, 'Fig4_Data.mat'))
load(fullfile(base, fpath, 'Fig4_Model.mat'))

% Plot 
new_pitch = str2num(char(pitch));
plot(new_pitch, rate,'-*', 'LineWidth',2);
hold on
plot(new_pitch, avBS*2.5,'-o', 'LineWidth',2)
set(gca, 'XScale', 'log')
box on
ylim([0 120])
xlim([57 575])
xticks([110 220 440 880])
xlabel('Fundamental Frequency (Hz)')
ylabel('Avg. Rate (sp/s)')
legend('Single-Unit Response', 'Model Response', 'Location','northwest')
set(gca,'FontSize',20)

%% Explainable variance explained calculation

R_int = corrcoef(rate,avBS.*3).^2;
R = R_int(1, 2);
[hat_r2er_BS, r2_BS] = r2er_n2m(avBS.*3, rate_matrix);
disp(['Variance explained = ' num2str(R) ', Explainable variance explained = '...
    num2str(hat_r2er_BS)]);



