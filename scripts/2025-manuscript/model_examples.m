%% ST_model_examples
clear 

%% Create figure

[base, datapath, savepath, ppi] = getPaths();

figure('Position',[1720,621,4.7*ppi,7*ppi])
h = gobjects(6, 1);

fontsize = 11;
legsize = 9;
labelsize = 20;
linewidth = 1.5;

%% Load in examples and plot

for ii = 1:6
	switch ii
		case 1 % BS Good
			putative = 'R24_TT2_P13_N05';
			CF = 1326;
			MTF_shape = 'BS';
		case 2 % BS Bad
			putative = 'R27_TT4_P8_N10';
			CF = 4652;
			MTF_shape = 'BS';
			% putative = 'R25_TT4_P7_N01';
			% CF = 1516;
		case 3 % Inhibition areas
			putative = 'R29_TT4_P5_N02';
			CF = 758;
			MTF_shape = 'BS';
			% putative = 'R24_TT2_P13_N02';
			% CF = 1150;
		case 4 % Hook
			putative = 'R29_TT3_P5_N07';
			CF = 1320;
			MTF_shape = 'BS';
			% putative = 'R25_TT3_P9_N01';
			% CF = 865;
		case 5 % BE Okay
			putative = 'R27_TT2_P8_N05';
			CF = 5278;
			MTF_shape = 'BE';
		case 6 % BE Bad
			putative = 'R29_TT1_P2_N04';
			CF = 5949;
			MTF_shape = 'BE';
	end
	ispl = 2;

	% Load in data
	[base, datapath, savepath, ppi] = getPaths();
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	% Load in model data
	modelpath = '/Volumes/Synth-Timbre/data/manuscript';
	load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE')
	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')

	% Load in lateral inhibition model
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')
	%load(fullfile(modelpath, 'model-fits', putative, ))

	% Get spont rate
	param_RM = data(2,2);
	data_RM = analyzeRM(param_RM{1});

	% Plot synthetic timbre (raw)
	h(ii) = subplot(3, 2, ii);
	param_ST = data(5+ispl, 2);
	if isempty(param_ST{1})
		continue
	end
	data_ST = analyzeST(param_ST, CF);
	data_ST = data_ST{1};

	% Z-score
	rate = data_ST.rate; %zscore(data_ST.rate);
	rate_sm = data_ST.rates_sm; %zscore(data_ST.rates_sm);
	hold on
	errorbar(data_ST.fpeaks,rate, data_ST.rate_std/sqrt(param_ST{1}.nrep),...
		'linewidth', 0.8, 'Color','k', 'LineStyle','none');
	plot(data_ST.fpeaks,rate_sm, 'linewidth', linewidth,'Color','k');
	yline(data_RM.spont, 'k')

	% Normalize and plot models
	%plot(data_ST.fpeaks, zscore(energy{ispl}.rate), 'LineWidth',linewidth, 'Color','#4634F1')
	%plot(data_ST.fpeaks, zscore(SFIE{ispl}.rate), 'LineWidth',linewidth, 'Color','#009E73')
	%plot(data_ST.fpeaks, zscore(lat_inh{ispl}.rate), 'LineWidth',linewidth, 'Color','#D55E00')

	% Scale models
	spont = data_RM.spont;
	max_rate = max(rate_sm);
	energy_rate = energy{ispl}.rate .* (max_rate/max(energy{ispl}.rate));
	SFIE_rate = SFIE{ispl}.rate .* (max_rate/max(SFIE{ispl}.rate));
	lat_inh_rate = lat_inh{ispl}.rate .* (max_rate/max(lat_inh{ispl}.rate));
	
	% Plot
	plot(data_ST.fpeaks, energy_rate, 'LineWidth',linewidth, 'Color','#4634F1')
	plot(data_ST.fpeaks, SFIE_rate, 'LineWidth',linewidth, 'Color','#009E73')
	plot(data_ST.fpeaks, lat_inh_rate, 'LineWidth',linewidth, 'Color','#D55E00')

	% % Annotate SFIE model R^2
	% message = sprintf('R^2 SFIE = %.02f', SFIE{ispl}.R2);
	% text(0.05, 0.25, message, 'Units', 'normalized', ...
	% 	'VerticalAlignment', 'top', 'FontSize',legsize, 'Color',...
	% 	'#009E73')
	% % Annotate energy model R^2
	% message = sprintf('R^2 Energy = %.02f', energy{ispl}.R2);
	% text(0.05, 0.15, message, 'Units', 'normalized', ...
	% 	'VerticalAlignment', 'top', 'FontSize',legsize, 'Color',...
	% 	'#4634F1')
	% % Annotate lateral inhibition model R^2
	% message = sprintf('R^2 Broad inh = %.02f', lat_inh{ispl}.R2);
	% text(0.05, 0.35, message, 'Units', 'normalized', ...
	% 	'VerticalAlignment', 'top', 'FontSize',legsize, 'Color',...
	% 	'#D55E00')

	% Plot parameters
	plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)];
	xline(CF, '--', 'Color',[0.7 0.7 0.7], 'linewidth', linewidth+1)
	xlim(plot_range)
	if ii == 1
		xticks([3200 3600 4000 4400 4800]);
	end
	ticks = xticks;
	xticklabels(ticks./1000)
	if ii == 5 || ii == 6
		xlabel('Spec. Peak Freq. (kHz)')
	end
	if ii == 1 || ii == 3 || ii == 5
		ylabel('Rate (sp/s)')
	end
	box on
	ylim([0 max_rate+max_rate*0.15])
	%ylim([-5 3])

	% BE/BS labels
	set(gca, 'fontsize', fontsize)
	text(0.05, 0.95, MTF_shape, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',legsize)
	grid on

	if ii == 1
		hLeg = legend('Data', '', 'Energy', 'SFIE', 'Broad Inh.', 'CF', 'Location',...
			'northoutside', 'FontSize',legsize, 'NumColumns', 2);
		hLeg.ItemTokenSize = [12, 12];
	end

end


%% Arrange

left = [0.13 0.6];
bottom = linspace(0.1, 0.68, 3);
height = 0.22;
width = 0.35;

set(h(1), 'Position', [left(1) bottom(3) width height])
set(h(2), 'Position', [left(2) bottom(3) width height])
set(h(3), 'Position', [left(1) bottom(2) width height])
set(h(4), 'Position', [left(2) bottom(2) width height])
set(h(5), 'Position', [left(1) bottom(1) width height])
set(h(6), 'Position', [left(2) bottom(1) width height])

% Annotate 
left = [0.01 0.51];
bottom = linspace(0.32, 0.91, 3);
annotation('textbox',[left(1) 0.97 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2) bottom(3) 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1) bottom(2) 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2) bottom(2) 0.0826 0.0385],'String',{'D'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1) bottom(1) 0.0826 0.0385],'String',{'E'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2) bottom(1) 0.0826 0.0385],'String',{'F'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
