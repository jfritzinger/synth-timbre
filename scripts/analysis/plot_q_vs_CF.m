%% Population Analysis
% J. Fritzinger, updated 1/9/23
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath, "peak_picking2.xlsx"));

%% Scatter plot of CF vs HHBW (all MTF)

% Plot
figure('position', [519,299,1150,539])
tiledlayout(2, 4, 'TileSpacing','tight', 'Padding','compact')
spls = [43, 63, 73, 83];
is200 = tables.F0 == 200;
isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');

for ibin = 1:2
	isbin = tables.binmode == ibin;

	for ispl = 1:4

		% Get data
		islevel = tables.SPL == spls(ispl);
		index = islevel & isbin & is200 & isMTF;

		% Data
		CFs = tables.CF(index);
		Qs = tables.Q(index);
		MTFs = tables.MTF(index);

		% Plot
		nexttile
		gscatter(CFs/1000, Qs, MTFs, 'filled')
		hold on
		box on

		% Find fit line 
		% Fit linear regression line
		mdl = fitlm(log10(CFs), log10(Qs));
		x = 0.3:0.5:10000;
		p(1) = mdl.Coefficients.Estimate(2,1);
		p(2) = mdl.Coefficients.Estimate(1,1);
		p(3) = mdl.Coefficients.pValue(2);
		p(4) = mdl.Rsquared.Ordinary;
		mdlfit(ibin, ispl,:) = 10.^(p(1)*log10(x)+p(2));
		mdlplot = squeeze(mdlfit(ibin, ispl, :));
		plot(x/1000, mdlplot, 'k');

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		title([num2str(spls(ispl)) ' dB SPL, n=' num2str(length(number))])
		xlabel('CF (kHz)')
		if ispl == 1
			ylabel('Q')
		end
		ylim([0.35 50])
		xlim([0.3 10])
		set(gca, 'fontsize', 17)
		set(gca, 'XScale', 'log')
		set(gca, 'YScale', 'log')
		xticks([0 200 500 1000 2000 5000 10000]/1000)
		yticks([0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000])
		grid on

	end
end

%% Fit lines 

figure('Position',[95,458,850,336])
tiledlayout(1, 2)

binmode = {'Contra', 'Binaural'};
for ibin = 1:2
	nexttile
	mdlplot = squeeze(mdlfit(ibin, :, :));
	plot(x/1000, mdlplot)

	ylim([0.9 11])
	xlim([0.3 10])
	set(gca, 'fontsize', 17)
	set(gca, 'XScale', 'log')
	set(gca, 'YScale', 'log')
	xticks([0 200 500 1000 2000 5000 10000]/1000)
	yticks([0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000])
	legend('43 dB SPL', '63 dB SPL', '73 dB SPL', '83 dB SPL', 'Location','northwest')
	grid on
	title([binmode{ibin} ' Linear Regression'])
	xlabel('CF (kHz)')
	ylabel('Q')
end