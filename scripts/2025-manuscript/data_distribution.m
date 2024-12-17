%% Fig0_DistributionOfData
% J. Fritzinger, updated 12/15/23
%
% This script loads in the putative neurons spreadsheet and plots the MTF
% distribution, BMF distribution, WMF distribution, hybrid BMF/WMF
% distribution, and CF distribution for all neurons
clear

%% Load in spreadsheet

[base, ~, savepath, ppi] = getPaths();
datapath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_units = size(sessions, 1);

%% Set up figure

figure('Position',[194,1045,1175,226])
tiledlayout(1, 5, 'Padding','compact')
fontsize = 12;
titlesize = 14;

%% Only get sessions with synthetic timbre 

synth(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
synth(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
synth(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
synth(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);

synth(:,5) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_100);
synth(:,6) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_100);
synth(:,8) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_100);

synth(:,9) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_con);
synth(:,10) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con);
synth(:,11) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con);
synth(:,12) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_con);

any_synth = any(synth(:,1:12), 2);
table = sessions(any_synth, :);

%% Get CFs for each putative neuron

% WBTIN Diotic
CFs = table.CF;
edges = [0 500 1000 2000 4000 8000 13000];
names = categorical({'<0.5', '0.5-1', '1-2', '2-4', '4-8', '8+'});
names = reordercats(names,{'<0.5', '0.5-1', '1-2', '2-4', '4-8', '8+'});
CF = CFs;
CF(CF==0) = [];
[N, edges1] = histcounts(CF, edges);

% Plot
nexttile
bar(names,N,'FaceColor', 'k', 'EdgeColor','k');
grid on
ylabel('# Neurons')
xlabel('CF (kHz)')
set(gca, 'FontSize', fontsize)
title('CF Distribution', 'fontsize', titlesize)
ylim([0 80])

%% Get MTFs for each putative neuron
% WBTIN Diotic

% figure
% target = cellfun(@(d) contains(d, 'R'), sessions.ST_83dB);
% num_sesh = sum(target);
% MTFs = sessions.MTF(target);

MTFs = table.MTF;
num_sesh = length(MTFs);
MTF_type = zeros(num_sesh,1);
for isesh = 1:num_sesh
	MTF_shape = MTFs{isesh};
	if contains(MTF_shape, 'H')
		MTF_type(isesh) = 3;
	elseif strcmp(MTF_shape, 'BE')
		MTF_type(isesh) = 1;
	elseif strcmp(MTF_shape, 'BS')
		MTF_type(isesh) = 2;
	else % Flat
		MTF_type(isesh) = 4;
	end
end
MTF_names = categorical({'BE','BS','Hybrid','Flat'});
MTF_names = reordercats(MTF_names,{'BE','BS','Hybrid','Flat'});
num_BE = sum(MTF_type==1);
num_BS = sum(MTF_type==2);
num_H = sum(MTF_type==3);
num_n = sum(MTF_type==4);
num_types = [num_BE num_BS num_H num_n];

% Plot
nexttile
bar(MTF_names,num_types, 'black');
set(gca, 'FontSize', fontsize)
title('MTF Type', 'FontSize',titlesize)
grid on
ylabel('# Neurons')
%ylim([0 110])


%% BMFs

% Get BMFs/WMFs
BE_MTFs = strcmp(table.MTF, 'BE');
BMFs = table.BMF(BE_MTFs);
BMFs(isnan(BMFs)) = [];

edges = [0.2 2 4 8 16 32 64 128 254 512 1028];
edges2 = zeros(10, 1);
for iedge = 1:10
	edges2(iedge) = sqrt(edges(iedge)*edges(iedge+1));
end

% BMFs
nexttile
histogram(BMFs, edges2,'FaceColor', '#0072BD', 'EdgeColor','k')
hold on
xline(exp(median(log(BMFs(BMFs~=0)))), 'k', 'LineWidth',1.5)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'FontSize', fontsize)
title('BE BMFs', 'fontsize', titlesize)
set(gca, 'XScale', 'log');
xlabel('BMF (Hz)')
ylabel('# Neurons')

%% WMFs

BS_MTFs = strcmp(table.MTF, 'BS');
WMFs = table.WMF(BS_MTFs);
WMFs(isnan(WMFs)) = [];

nexttile
histogram(WMFs, edges2,'FaceColor', '#D95319', 'EdgeColor','k')
hold on
ylabel('# Neurons')
xline(exp(median(log(WMFs(WMFs~=0)))), 'k', 'LineWidth',1.5)
set(gca, 'FontSize', fontsize)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'XScale', 'log');
xlabel('WMF (Hz)')
title('BS WMFs', 'fontsize', titlesize)

%% Hybrids

H_MTFs = contains(table.MTF, 'H');
WMFs = table.WMF(H_MTFs);
WMFs(isnan(WMFs)) = [];
BMFs = table.BMF(H_MTFs);
BMFs(isnan(BMFs)) = [];

% Plot 
nexttile
histogram(BMFs, edges2)
hold on
histogram(WMFs, edges2)
xlabel('BMF or WMF (Hz)')
hLegend = legend('BMF', 'WMF', 'Location','west');
hLegend.ItemTokenSize = [6,6];
set(gca, 'FontSize', fontsize)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'XScale', 'log');
ylabel('# Neurons')
title('Hybrids', 'fontsize', titlesize)

%% Export 

exportgraphics(gcf, fullfile(savepath, 'manuscript', 'data-distribution.png'), 'Resolution', 600)