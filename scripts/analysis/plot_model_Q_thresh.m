%% Calculate and plot model Q and thresholds
clear

%% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath, "model_Lat_Inh_Q_thresholds.xlsx"));

sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Plot Q

% Plot
figure('position', [296,665,1194,175])
tiledlayout(1, 6, 'TileSpacing','tight')
spls = [43, 63, 73, 83];


for ispl = 2

	% Get data
	islevel = tables.SPL == spls(ispl);
	index = islevel; % & isMTF;

	% Data
	CFs = tables.CF(index);
	Qs = tables.Q(index);
	MTFs = tables.MTF(index);

	% Plot
	nexttile
	scatter(CFs/1000, Qs, 'filled')
	hold on
	box on

	% Fit linear regression line
	mdl = fitlm(log10(CFs), log10(Qs));
	x = 0.3:0.5:10000;
	p(1) = mdl.Coefficients.Estimate(2,1);
	p(2) = mdl.Coefficients.Estimate(1,1);
	p(3) = mdl.Coefficients.pValue(2);
	p(4) = mdl.Rsquared.Ordinary;
	mdlfit(ispl,:) = 10.^(p(1)*log10(x)+p(2));
	mdlplot = squeeze(mdlfit(ispl, :));
	plot(x/1000, mdlplot, 'k');

	% Plot labels
	number = Qs;
	number(isnan(number)) = [];
	xlabel('CF (kHz)')
	if ispl == 1
		ylabel('Q')
	end
	ylim([0.35 50])
	xlim([0.3 10])
	set(gca, 'fontsize', 12)
	set(gca, 'XScale', 'log')
	set(gca, 'YScale', 'log')
	xticks([0 200 500 1000 2000 5000 10000]/1000)
	yticks([0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000])
	grid on

end

%%

nexttile 


nexttile



%% Plot change in Q over level
fontsize = 12;

% Find sessions for target synthetic timbre response
all_neurons = tables.Putative;
neurons = unique(all_neurons);
num_units = size(neurons, 1);

SPLs = [43, 63, 73, 83];
qs = NaN(num_units, 4);
for isesh = 1:num_units

	putative = neurons{isesh};
	isput = cellfun(@(s) strcmp(s, putative), tables.Putative);

	for ispl = 1:4
		ind = isput & tables.SPL==SPLs(ispl);
		if any(ind)
			qs(isesh, ispl) = tables.Q(ind);
			%qs_log(isesh, ispl) = tables.Q_log(ind);
			%CF_group(isesh) = tables.CF_Group(ind);
		end
	end
end


% Get matrix of units with 43, 63, 83 dB
qs2 = qs(:,[1,2,4]);
rows_with_nan = any(isnan(qs2),2);
qs2(rows_with_nan,:) = [];
x = [43, 63, 83];

% Try 1: Criteria using slope
for ii = 1:length(qs2)
	y = qs2(ii, :)';
	tbl = table(x', y, 'VariableNames', {'X', 'Q'});
	mdl = fitlm(tbl, 'Q ~ X');
	slopes(ii) = mdl.Coefficients.Estimate(2); % slope
end
criteria = 0.03;
same = slopes<criteria & slopes > -1*criteria;
decrease = slopes<-1*criteria;
increase = slopes>criteria;
spls = [43, 63, 83];

indices = [3, 6, 9];
for ii = 1:3
	if ii == 1
		values = increase;
		color = [27,158,119]/256;
	elseif ii == 2
		values = same;
		color = [217,95,2]/256;
	else
		values = decrease;
		color = [117,112,179]/256;
	end

	% Increase
	nexttile
	hold on
	if sum(qs2(values,:))~=0
	plot(spls, qs2(values,:)', 'color',color , 'LineWidth',1)
	xticks(spls)
	ylabel('Q')
	xlim([40 86])
	plot(spls, mean(qs2(values,:), 'omitnan'), 'k', 'LineWidth',2)
	plot(spls, median(qs2(values,:), 'omitnan'), ':k', 'LineWidth',2)
	set(gca, 'fontsize', fontsize)
	xlabel('Level (dB SPL)')
	end
	label = ['n=' num2str(sum(values))];
	text(0.05, 0.95, label, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)
	ylim([0 50])

end

%% Plot thresholds
%
% figure('Position',[1,630,1545,232])
% tiledlayout(1, 6)
% fontsize = 12;
%
% % Plot histogram
% spls = [43, 63, 73, 83];
% islevel = tables.SPL == spls(2);
% index = islevel; % & isMTF;
% thresholds = tables.Threshold(index);
% nexttile
% edges = linspace(0, 30, 20);
% histogram(thresholds, edges)
% xlabel('Thresholds (%)')
% ylabel('# Neurons')
% set(gca, 'fontsize', fontsize)
% %title('Histogram of Thresholds')
% hold on
% xline(4, 'r', 'LineWidth',2)
% xlim([0 30])
% box off
% grid on
% percent = sum(thresholds<=4)/length(thresholds)*100;
% legend([num2str(round(percent)) '% < 4% threshold'], '4% human threshold')
% title('Threshold Histogram')
%
% % Threshold vs CF plot
% spls = [43, 63, 73, 83];
% isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');
% nexttile
% ispl = 2;
% islevel = tables.SPL == spls(ispl);
% index = islevel; % & isMTF;
% CFs = tables.CF(index);
% Qs = tables.Threshold(index);
% MTFs = tables.MTF(index);
% Qs(isnan(Qs)) = 50;
% Qs(Qs>50) = 50;
% scatter(CFs/1000, Qs, 'filled')
% hold on
% scatter(1200/1000, 4, 'k', 'square', 'linewidth', 2) %'filled')
% yline(4, 'k')
% number = Qs;
% number(isnan(number)) = [];
% title('Thresholds vs CF, 63 dB SPL')
% xlabel('CF (kHz)')
% if ispl == 1
% 	ylabel('Q')
% end
% ylim([0.04 50])
% xlim([0.3 10])
% set(gca, 'fontsize', fontsize)
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% ylabel('Thresholds (%)')
% xticks([0 200 500 1000 2000 5000 10000]/1000)
% yticks([0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 70])
% yticklabels({'0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1', '2', '5', '10', '20', '>50'})
% grid on
% box off
%
% % Plot global changes
% nexttile
% spls = [43, 63, 73, 83];
% all_thresholds = NaN(4, 163);
%
% for ispl = 1:4
% 	% Get data
% 	islevel = tables.SPL == spls(ispl);
% 	index = islevel; % & isMTF;
%
% 	% Data
% 	thresh = tables.Threshold(index);
%
% 	% Add in units without thresholds
% 	thresh(isnan(thresh)) = 100;
%
% 	% Put into arrays
% 	all_thresholds(ispl,1:length(thresh)) = thresh;
%
% 	hold on
% 	swarmchart(ones(length(thresh), 1)*ispl, thresh)
% end
%
% boxplot(all_thresholds')
% ylim([0 50])
% ylabel('Threshold (%)')
% xlabel('Level (dB SPL)')
% xticklabels([43, 63, 73, 83])
% title('Population threshold changes over level')
% set(gca, 'fontsize', fontsize)
% set(gca, 'YScale', 'log')
%
% % Plot increasing/decreasing/same
% all_neurons = tables.Putative;
% neurons = unique(all_neurons);
% num_units = size(neurons, 1);
% SPLs = [43, 63, 73, 83];
% qs = NaN(num_units, 4);
% for isesh = 1:num_units
% 	putative = neurons{isesh};
% 	isput = cellfun(@(s) strcmp(s, putative), tables.Putative);
% 	for ispl = 1:4
% 		ind = isput & tables.SPL==SPLs(ispl);
% 		if any(ind)
% 			qs(isesh, ispl) = tables.Threshold(ind);
% 		end
% 	end
% end
% qs2 = qs(:,[1,2,4]);
% rows_with_nan = any(isnan(qs2),2);
% qs2(rows_with_nan,:) = [];
% x = [43, 63, 83];
% for ii = 1:length(qs2) % Try 1: Criteria using slope
% 	y = qs2(ii, :)';
%     tbl = table(x', y, 'VariableNames', {'X', 'Q'});
%     mdl = fitlm(tbl, 'Q ~ X');
%     slopes(ii) = mdl.Coefficients.Estimate(2); % slope
% end
% criteria = 0.03;
% same = slopes<criteria & slopes > -1*criteria;
% decrease = slopes<-1*criteria;
% increase = slopes>criteria;
% spls = [43, 63, 83];
% for ii = 1:3
% 	if ii == 1
% 		values = increase;
% 		color = [27,158,119]/256;
% 	elseif ii == 2
% 		values = same;
% 		color = [217,95,2]/256;
% 	else
% 		values = decrease;
% 		color = [117,112,179]/256;
% 	end
%
% 	nexttile
% 	hold on
% 	if sum(qs2(values,:))~=0
% 	plot(spls, qs2(values,:)', 'color',color , 'LineWidth',1)
% 	xticks(spls)
% 	ylabel('Thresholds (%)')
% 	xlim([40 86])
% 	plot(spls, mean(qs2(values,:), 'omitnan'), 'k', 'LineWidth',2)
% 	plot(spls, median(qs2(values,:), 'omitnan'), ':k', 'LineWidth',2)
% 	set(gca, 'fontsize', fontsize)
% 	xlabel('Level (dB SPL)')
%
% 	ylim([0 40])
% 	end
% 	label = ['n=' num2str(sum(values))];
% 	text(0.05, 0.95, label, 'Units', 'normalized', ...
% 		'VerticalAlignment', 'top', 'FontSize',fontsize)
% 	if ii == 1
% 		title('Increasing threshold')
% 	elseif ii == 2
% 		title('Steady thresholds')
% 	else
% 		title('Decreasing')
% 	end
%
%
% end

