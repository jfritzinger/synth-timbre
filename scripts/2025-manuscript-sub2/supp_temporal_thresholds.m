%% compare_thresholds
clear 

%% Load in data 

[base, datapath, ~, ppi] = getPaths();
datapath = fullfile(base, 'data', '2025-manuscript-sub2');
tables_VS = readtable(fullfile(datapath, "peak_picking_w_thresholds_VS.xlsx"));

[base, datapath, ~, ppi] = getPaths();
tables = readtable(fullfile(datapath, "peak_picking_w_thresholds.xlsx"));

sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% 

% Set up figure 
figure('Position',[50,50,4*ppi,4*ppi])
data_colors = {'#000000'};
legsize = 6;
fontsize = 7;
titlesize = 8;
linewidth = 1;
labelsize = 13;
scattersize = 10; 
capsize = 2;

%%

% Set up figure 
h(1) = subplot(2, 2, 1);
spls = [43, 63, 73, 83];
is200 = tables.F0==200;
isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');

for ibin = 2
	isbin = tables.binmode == ibin;
	for ispl = 2

		% Get data
		islevel = tables.SPL == spls(ispl);
		index = islevel & isbin & is200; % & isMTF;

		% Data
		CFs = tables_VS.CF(index);
		Qs = tables_VS.Threshold(index);

		% Add in units without thresholds
		Qs(isnan(Qs)) = 100;
		Qs(Qs>100) = 100;

		% Plot
		%gfig = gscatter(CFs/1000, Qs,peaks, 'filled');
		%set(gfig, 'MarkerEdgeColor', 'k');
		scatter(CFs/1000, Qs, scattersize, 'filled', 'MarkerEdgeColor','k', ...
			'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
		hold on

		% Plot minimum thresholds based on stimuli
		%plot(CFs_array/1000, min_threshold, 'r', 'LineWidth',2)

		% Plot human threshold
		scatter(1200/1000, 4,scattersize, 'r', 'filled','markeredgecolor', 'k') %'filled')
		yline(4, 'r')

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		%title('Thresholds vs CF')
		xlabel('CF (kHz)')
		if ispl == 1
			ylabel('Q')
		end
		ylim([0.35 100])
		xlim([0.3 10])
		set(gca, 'fontsize', fontsize)
		set(gca, 'XScale', 'log')
		set(gca, 'YScale', 'log')
		ylabel('Threshold from VS (%)')
		xticks([0 200 500 1000 2000 5000 10000]/1000)
		yticks([0.2 0.5 1 2 5 10 20 50 100])
		yticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '50', '>100'})
		grid on
		box off
		hleg = legend('Neural','Human', '', 'Location','southwest', 'box',...
			'off');
		hleg.ItemTokenSize = [8, 8];
	end
end

%% 

h(2) = subplot(2, 2, 2);

spls = [43, 63, 73, 83];
is200 = tables.F0==200;
isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');

for ibin = 2
	isbin = tables.binmode == ibin;
	for ispl = 2

		% Get data
		islevel = tables.SPL == spls(ispl);
		index = islevel & isbin & is200; % & isMTF;

		% Data
		Qs = tables.Threshold(index);
        Qs_VS = tables_VS.Threshold(index);

		% Add in units without thresholds
		Qs(isnan(Qs)) = 50;
		Qs(Qs>50) = 50;
        Qs_VS(isnan(Qs_VS)) = 50;
        Qs_VS(Qs_VS>50) = 50;

		% Plot
		scatter(Qs, Qs_VS, scattersize, 'filled', 'MarkerEdgeColor','k', ...
			'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
		hold on
        plot([0.35 50], [0.35 50], 'k')

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		%title('Thresholds vs CF')
		xlabel('Threshold from Rate (%)')
		if ispl == 1
			ylabel('Q')
		end
		ylim([0.35 50])
		xlim([0.35 50])
		set(gca, 'fontsize', fontsize)
		set(gca, 'XScale', 'log')
		set(gca, 'YScale', 'log')
		ylabel('Threshold from VS (%)')
		xticks([0.2 0.5 1 2 5 10 20 50 70])
		yticks([0.2 0.5 1 2 5 10 20 50 70])
		yticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '>50'})
        xticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '>50'})
		grid on
		box off
		hleg.ItemTokenSize = [8, 8];
	end
end

% Stats 
above_50 = sum(Qs_VS==50)/length(Qs_VS);
better_than_rate = sum(Qs_VS>Qs);
perc_better = better_than_rate/length(Qs_VS);
near_human =  sum(Qs_VS<6)/length(Qs_VS);

%% 

[base, datapath, ~, ppi] = getPaths();
tables_RIS = readtable(fullfile(datapath, "peak_picking_w_thresholds_RIS.xlsx"));

% Set up figure 
h(3) = subplot(2, 2, 3);
for ibin = 2
	isbin = tables_RIS.binmode == ibin;
	for ispl = 2

		% Get data
		islevel = tables_RIS.SPL == spls(ispl);
		index = islevel & isbin; % & isMTF;

		% Data
		CFs = tables_RIS.CF(index);
		Qs = tables_RIS.Threshold_Real(index);

		% Add in units without thresholds
		Qs(isnan(Qs)) = 100;
		Qs(Qs>100) = 100;

		% Plot
		%gfig = gscatter(CFs/1000, Qs,peaks, 'filled');
		%set(gfig, 'MarkerEdgeColor', 'k');
		scatter(CFs/1000, Qs, scattersize, 'filled', 'MarkerEdgeColor','k', ...
			'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
		hold on

		% Plot minimum thresholds based on stimuli
		%plot(CFs_array/1000, min_threshold, 'r', 'LineWidth',2)

		% Plot human threshold
		scatter(1200/1000, 4,scattersize, 'r', 'filled','markeredgecolor', 'k') %'filled')
		yline(4, 'r')

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		%title('Thresholds vs CF')
		xlabel('CF (kHz)')
		if ispl == 1
			ylabel('Q')
		end
		ylim([0.35 100])
		xlim([0.3 10])
		set(gca, 'fontsize', fontsize)
		set(gca, 'XScale', 'log')
		set(gca, 'YScale', 'log')
		ylabel('Threshold from VS (%)')
		xticks([0 200 500 1000 2000 5000 10000]/1000)
		yticks([0.2 0.5 1 2 5 10 20 50 100])
		yticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '50', '>100'})
		grid on
		box off
		hleg = legend('Neural','Human', '', 'Location','southwest', 'box',...
			'off');
		hleg.ItemTokenSize = [8, 8];
	end
end

%% 

h(4) = subplot(2, 2, 4);

spls = [43, 63, 73, 83];
is200 = tables.F0==200;
isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');

for ibin = 2
	isbin = tables.binmode == ibin;
	for ispl = 2

		% Get data
		islevel = tables.SPL == spls(ispl);
		index = islevel & isbin & is200; % & isMTF;

		% Data
		Qs = tables.Threshold(index);
        Qs_VS = tables_RIS.Threshold_Real;

		% Add in units without thresholds
		Qs(isnan(Qs)) = 50;
		Qs(Qs>50) = 50;
        Qs_VS(isnan(Qs_VS)) = 50;
        Qs_VS(Qs_VS>50) = 50;

		% Plot
		scatter(Qs, Qs_VS, scattersize, 'filled', 'MarkerEdgeColor','k', ...
			'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
		hold on
        plot([0.35 50], [0.35 50], 'k')

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		%title('Thresholds vs CF')
		xlabel('Threshold from Rate (%)')
		if ispl == 1
			ylabel('Q')
		end
		ylim([0.35 50])
		xlim([0.35 50])
		set(gca, 'fontsize', fontsize)
		set(gca, 'XScale', 'log')
		set(gca, 'YScale', 'log')
		ylabel('Threshold from VS (%)')
		xticks([0.2 0.5 1 2 5 10 20 50 70])
		yticks([0.2 0.5 1 2 5 10 20 50 70])
		yticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '>50'})
        xticklabels({'0.2', '0.5', '1', '2', '5', '10', '20', '>50'})
		grid on
		box off
		hleg.ItemTokenSize = [8, 8];
	end
end

% Stats 
above_50 = sum(Qs_VS==50)/length(Qs_VS);
better_than_rate = sum(Qs_VS<Qs);
perc_better = better_than_rate/length(Qs_VS);
near_human =  sum(Qs_VS<6)/length(Qs_VS);

%% Annotate 

annotation('textbox',[0.01 0.97 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.47 0.97 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.48 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.47 0.48 0.0826 0.0385],'String',{'D'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

%% Save figure

save_fig = 1;
if save_fig == 1
	filename = 'supp_temporal_thresholds';
	save_figure(filename)
end