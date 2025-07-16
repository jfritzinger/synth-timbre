%% supp_STRF_RM_plots
% This takes 5-10 minutes to run
clear

%% Load in data

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
sheetpath = '/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name), ...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);
savepath = 'STRF_Models';

%% Set up figure

figure('Position',[560,299,7*ppi,10*ppi])
tiledlayout(4, 2, 'Padding','compact')
fontsize = 12;
labelsize = 20;

%% Examples

index1 = [1, 3, 5];
index2 = [2, 4, 6];
ibin = 2;
ispl = 1;
SPLs = {'43', '63', '73'};
binmodes = {'Contra', 'Bin'};
for isesh = 1:3
	switch isesh
		case 1
			putative = 'R24_TT2_P13_N02';
			CF = 1150;
		case 2
			putative = 'R24_TT2_P13_N06';
			CF = 1509;
		case 3
			putative = 'R27_TT4_P8_N05';
			CF = 3983;
	end

	% Load in data
	load(fullfile(datapath, 'neural_data', [putative '.mat']), 'data');
	params_STRF = data{4,ibin};
	param_ST = data(6+ispl,ibin);
	param_RM = data{2, 2};

	% Load in STRF model
	if ispl == 1 && ibin == 1
		filename = [putative '_STRF_63_Contra'];
		load(fullfile(datapath, savepath, [filename '.mat']), "STRFmodel_con")
	elseif ispl == 1 && ibin == 2
		filename = [putative '_STRF_63_Bin'];
		load(fullfile(datapath, savepath, [filename '.mat']), "STRFmodel")
	elseif ispl == 2 && ibin == 1
		filename = [putative '_STRF_73_Contra'];
		load(fullfile(datapath, savepath, [filename '.mat']), "STRFmodel_con")
	elseif ispl == 2 && ibin == 2
		filename = [putative '_STRF_73_Bin'];
		load(fullfile(datapath, savepath, [filename '.mat']), "STRFmodel")
	end

	% Load in RM analysis
	filename = sprintf('%s_RM_%s_%s', putative, SPLs{ispl}, binmodes{ibin});
	load(fullfile(datapath, 'RM_predictions', [filename '.mat']), "RM_R2")

	% General analysis
	data_STRF = analyzeSTRF(params_STRF);
	data_ST = analyzeST(param_ST, CF);
	param_ST = param_ST{1};
	data_ST = data_ST{1};

	% Plot STRF
	h(index1(isesh)) = subplot(4, 2, index1(isesh));
	STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
	imagesc(data_STRF.t, data_STRF.f./1000, STRF_mat, data_STRF.clims_strf);
	set(gca,'Ydir','normal','XLim',data_STRF.tlims,'YLim',[param_ST.fpeaks(2) ...
	param_ST.fpeaks(end)]./1000)
	hold on
	yline(CF/1000, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
	colormap(redblue)
	grid on
	xlabel('Time (s)');
	ylabel('Frequency (kHz)')
	set(gca, 'fontsize', fontsize)


	% Plot real response
	h(index2(isesh)) = subplot(4, 2, index2(isesh));
	hold on
	errorbar(data_ST.fpeaks,data_ST.rate,data_ST.rate_std/(sqrt(param_ST.nrep)),...
		'LineWidth',1.5);
	plot(data_ST.fpeaks,(STRFmodel.avModel.*STRFmodel.ratio), 'LineWidth',1.5);
	plot(data_ST.fpeaks, RM_R2.rate_RM, 'LineWidth',1.5);
	xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
	xlim([param_ST.fpeaks(1) param_ST.fpeaks(end)]);
	grid on
	xlabel('Tone Frequency (Hz)')
	ylabel('Avg Rate (sp/s)')
	ylim([0 STRFmodel.max_all+7])
	xticklabels(xticks/1000)
	if isesh == 1
		hLeg = legend({'Data', ['STRF' newline 'Model'],...
		'RM'}, 'Location','best');
		hLeg.ItemTokenSize = [6,6];
	end
	% hLeg = legend({'Data', sprintf('STRF, R^2=%0.2f', STRFmodel.R2),...
	% 	sprintf('RM, R^2=%0.2f',RM_R2.R2)}, 'Location','eastoutside');
	% hLeg.ItemTokenSize = [6,6];
	set(gca, 'fontsize', fontsize)

end

%% Box plots

load(fullfile(datapath, 'RM_predictions', 'AllUnits_Bin.mat'),...
	"R2", 'CF')
R2(R2==0) = NaN;
SPLs = [43, 63, 73, 83];
h(7) = subplot(4, 2, 7);
boxplot(R2')
hold on
num = size(R2,2);
x = [ones(1,num); 2*ones(1,num); 3*ones(1,num); 4*ones(1,num)];
swarmchart(x', R2')
ylim([0 1])
xticklabels(SPLs)
xlabel('Level (dB SPL)')
ylabel('R^2 RM')
set(gca, 'fontsize', fontsize)
box off

%% Scatter plots


load(fullfile(datapath, 'STRF_Models', 'AllUnits_63_Bin.mat'),"R2_STRF", 'R2_RM')
R2_RM(R2_RM==0) = NaN;
R2_STRF(R2_STRF==0) = NaN;
SPLs = {'63', '73'};
ispl = 1;
h(8) = subplot(4, 2, 8);
scatter(R2_RM(ispl, :), R2_STRF(ispl, :), 'filled', 'MarkerEdgeColor','k');
xlabel('R^2 RM')
ylabel('R^2 STRF')
set(gca, 'fontsize', fontsize)

%% Annotations 
left = [0.01 0.48];
bottom = [0.23 linspace(0.5, 0.96, 3)];

annotation('textbox',[left(1) bottom(4) 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1) bottom(3) 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1) bottom(2) 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1) bottom(1) 0.0826 0.0385],'String',{'D'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2) bottom(1) 0.0826 0.0385],'String',{'E'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

% Move plots 
left = [0.12 0.58];
bottom = [0.06 linspace(0.33, 0.8, 3)];
width = 0.35;
height = 0.17;

set(h(1), 'position', [left(1) bottom(4) width height])
set(h(2), 'position', [left(2) bottom(4) width height])
set(h(3), 'position', [left(1) bottom(3) width height])
set(h(4), 'position', [left(2) bottom(3) width height])
set(h(5), 'position', [left(1) bottom(2) width height])
set(h(6), 'position', [left(2) bottom(2) width height])
set(h(7), 'position', [left(1) bottom(1) width height])
set(h(8), 'position', [left(2) bottom(1) width height])

%% Save figure 

filename = 'FigS3_STRF_RM_plots';
saveFigure(filename)