%% plot_pop_thresholds
clear 

%% Load in data 

[base, datapath, ~, ppi] = getPaths();
tables = readtable(fullfile(datapath, "peak_picking_w_thresholds.xlsx"));

sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Set up figure 

figure('Position',[185,556,10*ppi,4*ppi])
tiledlayout(1, 3)
data_colors = {'#03882F', '#82BB95'};
legsize = 10;
linewidth = 2;
fontsize = 12;
labelsize = 20;

%% Example unit 

examples = {'R24_TT2_P13_N05', 'R27_TT2_P8_N02', 'R27_TT2_P8_N05', ...
	'R25_TT2_P9_N01', 'R27_TT3_P1_N08', 'R27_TT2_P7_N01', ...
	'R29_TT4_P5_N15', 'R25_TT2_P8_N02', 'R29_TT1_P3_N05'};
ineuron = 1;
h(1) = subplot(1, 3, 1);

% Load in data
putative = examples{ineuron};
[base, datapath, savepath, ppi] = getPaths();
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

% RM to get spont
params_RM = data{2, 2};
data_RM = analyzeRM(params_RM);
spont = data_RM.spont;

% Synthetic timbre analysis
params = data(7, 2);
params = params(~cellfun(@isempty, params));
data_ST  = analyzeST(params, CF);
data_ST = data_ST{1};
rate = data_ST.rate;
rate_std = data_ST.rate_std;
rlb = data_ST.rlb;
rub = data_ST.rub;
fpeaks = data_ST.fpeaks;
spl = data_ST.spl;
rate_sm = data_ST.rates_sm;
max_rate = max(rate);

% Plot
fpeaks_re_CF = log2(fpeaks/CF);
hold on
rates_sm = smooth_rates(rate, rlb, rub, CF);
errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
	'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{1})
plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{1})
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)

% Figure parameters
plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
set(gca, 'Fontsize', fontsize)
xlim(plot_range);
grid on
ylim([0 max_rate+5])
ylabel('Avg. Rate (sp/s)')
xlabel('Spectral Peak Freq. (Hz)')
xlim([0.4 2.4])
%title('Example Unit with Threshold')

[threshold_percent, threshold_freq, slope_rate] = calculateThresholds(fpeaks, rate, rate_std, CF);
plot(threshold_freq/1000, slope_rate, 'k', 'LineWidth',3)
xline(threshold_freq(1)/1000, 'k')
xline(threshold_freq(2)/1000, 'k')

% Annotate 
msg = sprintf('Threshold = %0.02f%%', threshold_percent);
text(0.6, 0.95, msg, 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',legsize)

%% Plot histogram 

% Get data
spls = [43, 63, 73, 83];
is200 = tables.F0==200;
isbin = tables.binmode == 2;
islevel = tables.SPL == spls(2);
index = islevel & isbin & is200; % & isMTF;

% Data
thresholds = tables.Threshold(index);

% Plot
h(2) = subplot(1, 3, 2);
edges = linspace(0, 30, 20);
histogram(thresholds, edges)
xlabel('Thresholds (%)')
ylabel('# Neurons')
set(gca, 'fontsize', fontsize)
%title('Histogram of Thresholds')
hold on
xline(4, 'r', 'LineWidth',2)
xlim([0 30])
ylim([0 32])
box off
grid on

% Percent <4% threshold
percent = sum(thresholds<=4)/length(thresholds)*100;
legend([num2str(round(percent)) '% < 4% threshold'], '4% human threshold')

%% Plot scatter of thresholds vs CF

% Get minimum thresholds 
% CFs_array = linspace(200, 10000, 100);
% min_threshold = 35./CFs_array*100;


% Set up figure 
h(3) = subplot(1, 3, 3);


%% Thresholds over level  

% Get minimum thresholds 
% CFs_array = linspace(200, 10000, 100);
% min_threshold = 35./CFs_array*100;


% Set up figure 
%h(4) = subplot(1, 4, 4);
spls = [43, 63, 73, 83];
is200 = tables.F0==200;
isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');

figure
all_thresholds = NaN(4, 163);
for ibin = 2
	isbin = tables.binmode == ibin;
	for ispl = 1:4

		% Get data
		islevel = tables.SPL == spls(ispl);
		index = islevel & isbin & is200; % & isMTF;

		% Data
		CFs = tables.CF(index);
		thresh = tables.Threshold(index);
		MTFs = tables.MTF(index);
		peaks = tables.Type(index);

		% Add in units without thresholds
		thresh(isnan(thresh)) = 100;

		% Put into arrays
		all_thresholds(ispl,1:length(thresh)) = thresh;

		hold on
		swarmchart(ones(length(thresh), 1)*ispl, thresh)

	end
end
% 
% figure
% %swarmchart()
% swarmchart(1:4, all_thresholds, 'filled')
boxplot(all_thresholds')
ylim([0 50])
ylabel('Threshold (%)')
xlabel('Level (dB SPL)')
xticklabels([43, 63, 73, 83])
title('Population threshold changes over level')
set(gca, 'fontsize', fontsize)

% [p,tbl,stats] = anova1(all_thresholds');
% results = multcompare(stats);
kruskalwallis(all_thresholds')

[p, tbl, stats] = kruskalwallis(all_thresholds', 1:4);
%multcompare(stats, 'CType', 'dunn-sidak');

%%

% Find sessions for target synthetic timbre response
all_neurons = tables.Putative;
neurons = unique(all_neurons);
num_units = size(neurons, 1);
isbin = tables.binmode == 2;
is200 = tables.F0 == 200;

SPLs = [43, 63, 73, 83];
qs = NaN(num_units, 4);
for isesh = 1:num_units

	putative = neurons{isesh};
	isput = cellfun(@(s) strcmp(s, putative), tables.Putative);

	for ispl = 1:4
		ind = isput & isbin & is200 & tables.SPL==SPLs(ispl);
		if any(ind)
			qs(isesh, ispl) = tables.Threshold(ind);
			CF_group(isesh) = tables.CF_Group(ind);
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

figure('Position',[185,556,10*ppi,4*ppi])
tiledlayout(1, 3)
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
	%h(indices(ii)) = subplot(3, 3, indices(ii));
	nexttile
	hold on
	plot(spls, qs2(values,:)', 'color',color , 'LineWidth',linewidth)
	xticks(spls)
	ylabel('Thresholds (%)')
	xlim([40 86])
	plot(spls, mean(qs2(values,:), 'omitnan'), 'k', 'LineWidth',2)
	plot(spls, median(qs2(values,:), 'omitnan'), ':k', 'LineWidth',2)
	set(gca, 'fontsize', fontsize)
	xlabel('Level (dB SPL)')

	label = ['n=' num2str(sum(values))];
	text(0.05, 0.95, label, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)
	hLeg = legend;
	num_lines = size(hLeg.String,2);
	for iii = 1:num_lines
		if iii==num_lines
			leg{iii} = 'Mean';
		elseif iii == num_lines-1
			leg{iii} = 'Median';
		else
			leg{iii} = '';
		end
	end
	hLeg = legend(leg, 'FontSize',legsize);
	hLeg.ItemTokenSize = [15,6];
	hLeg.Box = 'off';
	ylim([0 40])
	if ii == 1
		title('Increasing threshold')
	elseif ii == 2
		title('Steady thresholds')
	else
		title('Decreasing')
	end


end



%% Arrange and annotate 

left = linspace(0.07, 0.75, 3);
width = 0.23;
height = 0.73;
for ii = 1:3
	set(h(ii), 'Position', [left(ii) 0.15 width height])
end

% Set annotations
left = linspace(0.01, 0.69, 3);
annotation('textbox',[left(1) 0.95 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2) 0.95 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(3) 0.95 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');




